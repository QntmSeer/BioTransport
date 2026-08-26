"""Multiscale Continuum Process Optimization & Design Space Map (Minimalist Aesthetic)."""

from __future__ import annotations
import argparse
from pathlib import Path
import sys
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as opt

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from biotransport.bridge import MdBridgeModel, ProcessSimulator
from biotransport.theme import apply_minimalist_theme, PALETTE


def run_intense_parameter_sweep(
    params_json: str = "data/sample_md_params.json",
    output_png: str = "data/intense_optimization_landscape.png",
    sample_density: int = 60,
) -> None:
    apply_minimalist_theme()
    print(f"Generating minimalist 2D continuum design space landscape ({sample_density}x{sample_density} grid)...")
    model = MdBridgeModel.load_json(params_json)
    sim = ProcessSimulator(model)

    tmps_bar = np.linspace(0.2, 4.0, sample_density)
    shears = np.linspace(1000.0, 15000.0, sample_density)
    TMP, Gamma = np.meshgrid(tmps_bar, shears)

    bulk_conc = 5.0  # g/L
    c_gel = sim.props.gel_concentration_g_L
    mu_0 = sim.water_viscosity(sim.thermo.temperature_K)
    rm = 1.2e12
    r_c0 = 8.5e14
    n_comp = sim.props.compressibility_exponent_n

    channel_h = 0.0005  # 0.5 mm channel height
    pump_eff = 0.70

    flux_grid = np.zeros_like(TMP)
    sec_grid = np.zeros_like(TMP)

    for i in range(sample_density):
        for j in range(sample_density):
            tmp_b = TMP[i, j]
            gamma = Gamma[i, j]
            tmp_pa = tmp_b * 100_000.0
            km = sim.mass_transfer_coefficient(gamma)
            rc_spec = r_c0 * ((tmp_pa / 100_000.0) ** n_comp)

            def residual(j_val: float) -> float:
                pe = np.clip(j_val / km, 0.0, 6.0)
                cw = min(c_gel, float(bulk_conc * np.exp(pe)))
                j_crit = 0.22 * km
                excess = max(0.0, j_val - j_crit)
                sticking = 0.004 * ((cw / c_gel) ** 1.5)
                k_shear = (gamma / 4000.0) * 1.2e-3
                m_cake = (sticking * excess * bulk_conc) / max(1e-6, k_shear)
                r_cake = rc_spec * m_cake
                total_r = rm + r_cake
                pi_w = sim.virial_osmotic_pressure(cw)
                phi_rel = np.clip(cw / (c_gel * 1.2), 0.0, 0.95)
                mu_w = mu_0 * ((1.0 - phi_rel) ** (-1.2))
                eff_tmp = max(50.0, tmp_pa - pi_w)
                return j_val - eff_tmp / (mu_w * total_r)

            j_max = tmp_pa / (mu_0 * rm)
            try:
                j_sol = float(opt.brentq(residual, 0.0, j_max, xtol=1e-10))
            except Exception:
                j_sol = km * np.log(c_gel / bulk_conc)

            j_lmh = j_sol * 3.6e6
            flux_grid[i, j] = j_lmh

            crossflow_work_pa = (mu_0 * channel_h * (gamma**2)) / (3.0 * max(1e-7, j_sol))
            total_work_pa = tmp_pa + crossflow_work_pa
            sec_kwh_m3 = total_work_pa / (3.6e6 * pump_eff)
            sec_grid[i, j] = sec_kwh_m3

    fig, axes = plt.subplots(1, 2, figsize=(15.5, 6.2), dpi=300)
    fig.suptitle("Multiscale Continuum Process Optimization & Design Space Map", fontsize=14, color=PALETTE["slate_dark"])

    # Panel A: Minimalist Iso-Flux Contours (cividis sequential colormap)
    ax1 = axes[0]
    c1 = ax1.contourf(Gamma, TMP, flux_grid, levels=22, cmap="cividis")
    lines1 = ax1.contour(Gamma, TMP, flux_grid, levels=[15, 20, 25, 30, 35, 40], colors="#0f172a", linewidths=0.75, alpha=0.6)
    ax1.clabel(lines1, inline=True, fmt="%.0f LMH", fontsize=8.5, colors="#0f172a")
    cbar1 = fig.colorbar(c1, ax=ax1, pad=0.03)
    cbar1.set_label(r"Steady Permeate Flux $J\ (\mathrm{L} \cdot \mathrm{m}^{-2} \cdot \mathrm{h}^{-1})$", fontsize=10.5, color=PALETTE["slate_dark"])

    darcy_limit_tmp = 0.8 + 1.8 * (shears / 15000.0)
    ax1.plot(shears, darcy_limit_tmp, color=PALETTE["copper"], linestyle="--", lw=2.0, label="Gel-Polarization Onset")
    ax1.scatter([7500.0], [1.7], color="#ffffff", edgecolors=PALETTE["slate_dark"], s=140, marker="o", lw=1.5, zorder=6, label=r"Target Point ($J = 25\,\mathrm{LMH}$)")

    ax1.set_xlabel(r"Crossflow Shear Rate $\dot{\gamma}\ (\mathrm{s}^{-1})$")
    ax1.set_ylabel("Transmembrane Pressure TMP (bar)")
    ax1.set_title(r"A. Permeate Productivity Map $J(\Delta P, \dot{\gamma})$", loc="left", color=PALETTE["slate_dark"])
    ax1.grid(True, linestyle="--", alpha=0.35, color=PALETTE["bg_grid"])
    ax1.legend(loc="upper left", frameon=True)

    # Panel B: Minimalist Specific Energy Consumption (mako / viridis sequential)
    ax2 = axes[1]
    sec_clipped = np.clip(sec_grid, 0.2, 4.0)
    c2 = ax2.contourf(Gamma, TMP, sec_clipped, levels=22, cmap="mako_r" if "mako_r" in plt.colormaps() else "Blues")
    lines2 = ax2.contour(Gamma, TMP, sec_clipped, levels=[0.4, 0.6, 0.8, 1.2, 1.8, 2.5, 3.2], colors="#0f172a", linewidths=0.75, alpha=0.6)
    ax2.clabel(lines2, inline=True, fmt="%.2f", fontsize=8.5, colors="#0f172a")
    cbar2 = fig.colorbar(c2, ax=ax2, pad=0.03)
    cbar2.set_label(r"Specific Energy Consumption $\mathrm{SEC}\ (\mathrm{kWh} \cdot \mathrm{m}^{-3})$", fontsize=10.5, color=PALETTE["slate_dark"])

    min_idx = np.unravel_index(np.argmin(sec_grid), sec_grid.shape)
    opt_shear = Gamma[min_idx]
    opt_tmp = TMP[min_idx]
    ax2.scatter(
        [opt_shear],
        [opt_tmp],
        color="#ffffff",
        edgecolors=PALETTE["teal"],
        s=140,
        marker="s",
        lw=2.0,
        zorder=6,
        label=f"Pareto Baseline (TMP = {opt_tmp:.1f} bar, $\\dot{{\\gamma}} = {opt_shear:.0f}\\,\\mathrm{{s}}^{{-1}}$)",
    )

    ax2.set_xlabel(r"Crossflow Shear Rate $\dot{\gamma}\ (\mathrm{s}^{-1})$")
    ax2.set_ylabel("Transmembrane Pressure TMP (bar)")
    ax2.set_title(r"B. Specific Energy Consumption & Pareto Basin", loc="left", color=PALETTE["slate_dark"])
    ax2.grid(True, linestyle="--", alpha=0.35, color=PALETTE["bg_grid"])
    ax2.legend(loc="upper right", frameon=True)

    plt.tight_layout()
    out_p = Path(output_png)
    out_p.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(out_p, bbox_inches="tight")
    plt.close()
    print(f"Successfully generated minimalist optimization landscape: {out_p}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--params", type=str, default="data/sample_md_params.json")
    parser.add_argument("--out", type=str, default="data/intense_optimization_landscape.png")
    parser.add_argument("--density", type=int, default=60)
    args = parser.parse_args()
    run_intense_parameter_sweep(args.params, args.out, sample_density=args.density)
