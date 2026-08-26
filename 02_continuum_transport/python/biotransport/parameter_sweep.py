"""High-intensity parameter sweep and honeycomb optimization landscape."""

from __future__ import annotations
import argparse
from pathlib import Path
import sys
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as opt

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from biotransport.bridge import MdBridgeModel, ProcessSimulator


def run_intense_parameter_sweep(
    params_json: str = "data/sample_md_params.json",
    output_png: str = "data/intense_optimization_landscape.png",
    sample_density: int = 70,
    hex_grid: int = 24,
) -> None:
    print(f"Running high-intensity honeycomb parameter exploration ({sample_density}x{sample_density} dense samples)...")
    model = MdBridgeModel.load_json(params_json)
    sim = ProcessSimulator(model)

    tmps_bar = np.linspace(0.2, 4.0, sample_density)
    shears = np.linspace(1000.0, 15000.0, sample_density)

    TMP_mesh, Shear_mesh = np.meshgrid(tmps_bar, shears)
    TMP_flat = TMP_mesh.flatten()
    Shear_flat = Shear_mesh.flatten()

    bulk_conc = 5.0  # g/L
    c_gel = sim.props.gel_concentration_g_L
    mu_0 = sim.water_viscosity(sim.thermo.temperature_K)
    rm = 1.2e12  # m^-1
    r_c0 = 8.5e14  # m/kg
    n_comp = sim.props.compressibility_exponent_n  # 0.45

    channel_h = 0.0005  # 0.5 mm channel height
    channel_l = 0.10  # 10 cm channel length

    fluxes = []
    energy_efficiencies = []

    for tmp_b, gamma in zip(TMP_flat, Shear_flat):
        tmp_pa = tmp_b * 100_000.0
        km = sim.mass_transfer_coefficient(gamma)

        # Compressible cake resistance scaling with pressure: r_c = r_c0 * (TMP / 1 bar)^n
        rc_spec = r_c0 * ((tmp_pa / 100_000.0) ** n_comp)

        def residual(j_val: float) -> float:
            pe = np.clip(j_val / km, 0.0, 6.0)
            cw = min(c_gel, float(bulk_conc * np.exp(pe)))

            j_crit = 0.25 * km
            excess = max(0.0, j_val - j_crit)
            sticking = 0.003 * ((cw / c_gel) ** 1.5)
            k_shear = (gamma / 4000.0) * 1.2e-3
            m_cake = (sticking * excess * bulk_conc) / max(1e-6, k_shear)
            r_cake = rc_spec * m_cake

            total_r = rm + r_cake
            pi_w = sim.virial_osmotic_pressure(cw)
            phi_rel = np.clip(cw / (c_gel * 1.2), 0.0, 0.95)
            mu_w = mu_0 * ((1.0 - phi_rel) ** (-1.2))
            eff_tmp = max(100.0, tmp_pa - pi_w)
            return j_val - eff_tmp / (mu_w * total_r)

        j_max = tmp_pa / (mu_0 * rm)
        try:
            j_sol = float(opt.brentq(residual, 0.0, j_max, xtol=1e-10))
        except Exception:
            j_sol = km * np.log(c_gel / bulk_conc)
        j_lmh = j_sol * 3.6e6
        fluxes.append(j_lmh)

        # Total energy expenditure: Transmembrane pumping + Crossflow viscous dissipation
        # Crossflow pressure drop: delta_P_cf = mu_0 * L * gamma / h (Pa)
        dp_crossflow_bar = (mu_0 * channel_l * gamma / channel_h) / 100_000.0
        # Total hydraulic power equivalent (bar)
        equivalent_power_bar = tmp_b + 0.85 * (dp_crossflow_bar * 12.0) ** 1.3
        energy_eff = j_lmh / max(0.1, equivalent_power_bar)
        energy_efficiencies.append(energy_eff)

    fluxes = np.array(fluxes)
    energy_efficiencies = np.array(energy_efficiencies)

    fig, axes = plt.subplots(1, 2, figsize=(15, 6), dpi=300)
    fig.suptitle("Coupled Multiscale Continuum Optimization (Honeycomb Tessellation)", fontsize=14, fontweight="bold")

    # Plot 1: Hexbin Honeycomb Permeate Flux
    ax1 = axes[0]
    hb1 = ax1.hexbin(
        Shear_flat,
        TMP_flat,
        C=fluxes,
        gridsize=hex_grid,
        cmap="viridis",
        edgecolors="#222222",
        linewidths=0.35,
        mincnt=1,
    )
    cbar1 = fig.colorbar(hb1, ax=ax1, pad=0.03)
    cbar1.set_label(r"Steady Permeate Flux $(\mathrm{L} \cdot \mathrm{m}^{-2} \cdot \mathrm{h}^{-1})$", fontsize=11)
    ax1.set_xlabel(r"Crossflow Shear Rate $\dot{\gamma}\ (\mathrm{s}^{-1})$", fontsize=11)
    ax1.set_ylabel("Transmembrane Pressure TMP (bar)", fontsize=11)
    ax1.set_title(r"A. Honeycomb Productivity Map $J_{\mathrm{steady}}(\Delta P, \dot{\gamma})$", fontsize=12, fontweight="bold", loc="left")
    ax1.grid(True, linestyle="--", alpha=0.3)

    # Plot 2: Hexbin Honeycomb Energy-Optimal Envelope (with Pareto sweet spot)
    ax2 = axes[1]
    hb2 = ax2.hexbin(
        Shear_flat,
        TMP_flat,
        C=energy_efficiencies,
        gridsize=hex_grid,
        cmap="plasma",
        edgecolors="#222222",
        linewidths=0.35,
        mincnt=1,
    )
    cbar2 = fig.colorbar(hb2, ax=ax2, pad=0.03)
    cbar2.set_label(r"Specific Energy Efficiency $(\mathrm{LMH} \cdot \mathrm{bar}^{-1})$", fontsize=11)
    ax2.set_xlabel(r"Crossflow Shear Rate $\dot{\gamma}\ (\mathrm{s}^{-1})$", fontsize=11)
    ax2.set_ylabel("Transmembrane Pressure TMP (bar)", fontsize=11)
    ax2.set_title(r"B. Energy-Optimal Operating Envelope (Pareto Sweet Spot)", fontsize=12, fontweight="bold", loc="left")
    ax2.grid(True, linestyle="--", alpha=0.3)

    plt.tight_layout()
    out_p = Path(output_png)
    out_p.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(out_p, bbox_inches="tight")
    plt.close()
    print(f"Successfully generated coupled honeycomb landscape: {out_p}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--params", type=str, default="data/sample_md_params.json")
    parser.add_argument("--out", type=str, default="data/intense_optimization_landscape.png")
    parser.add_argument("--density", type=int, default=70)
    args = parser.parse_args()
    run_intense_parameter_sweep(args.params, args.out, sample_density=args.density)
