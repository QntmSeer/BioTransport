"""Publication-ready visualization routines with minimalist academic styling."""

from __future__ import annotations
from pathlib import Path
from typing import Any
import numpy as np
import matplotlib.pyplot as plt
from biotransport.theme import apply_minimalist_theme, PALETTE


def plot_multiscale_results(
    results: dict[str, Any],
    output_path: Path | str,
) -> None:
    apply_minimalist_theme()
    t_min = np.asarray(results["time_s"]) / 60.0
    flux = np.asarray(results["flux_lmh"])
    cw = np.asarray(results["wall_conc_g_l"])
    rm = results.get("membrane_resistance_m_inv", 1.0e12)
    rc = np.asarray(results.get("cake_resistance_m_inv", np.zeros_like(flux)))
    cb = results.get("bulk_conc_g_l", 5.0)
    c_gel = results.get("gel_conc_g_l", 420.0)

    fig, axes = plt.subplots(2, 2, figsize=(14, 9), dpi=300)
    fig.suptitle("Continuum Bioparticle Transport & Fouling Dynamics", fontsize=14, color=PALETTE["slate_dark"])

    # Panel A: Flux Decline
    ax1 = axes[0, 0]
    ax1.plot(t_min, flux, color=PALETTE["slate_dark"], lw=2.2, label=r"Permeate Flux $J(t)$")
    ax1.set_xlabel("Filtration Time (min)")
    ax1.set_ylabel(r"Permeate Flux $(\mathrm{L} \cdot \mathrm{m}^{-2} \cdot \mathrm{h}^{-1})$")
    ax1.set_title("A. Permeate Flux Decline", loc="left", color=PALETTE["slate_dark"])
    ax1.legend(loc="upper right", frameon=True)

    # Panel B: Wall Concentration
    ax2 = axes[0, 1]
    ax2.plot(t_min, cw, color=PALETTE["teal"], lw=2.2, label=r"Membrane Wall $C_w(t)$")
    ax2.axhline(y=cb, color=PALETTE["slate_light"], linestyle="--", lw=1.5, label=f"Bulk Feed $C_b = {cb:.1f}\\,\\mathrm{{g/L}}$")
    ax2.axhline(y=c_gel, color=PALETTE["copper"], linestyle=":", lw=1.5, label=f"Gel Boundary $C_{{\\mathrm{{gel}}}} = {c_gel:.0f}\\,\\mathrm{{g/L}}$")
    ax2.set_xlabel("Filtration Time (min)")
    ax2.set_ylabel(r"Concentration $(\mathrm{g} \cdot \mathrm{L}^{-1})$")
    ax2.set_title("B. Polarization Dynamics", loc="left", color=PALETTE["slate_dark"])
    ax2.legend(loc="center right", frameon=True)

    # Panel C: Resistance Breakdown
    ax3 = axes[1, 0]
    ax3.plot(t_min, np.full_like(t_min, rm) / 1e12, color=PALETTE["slate_light"], linestyle="--", lw=1.8, label=r"Membrane Resistance $R_m$")
    ax3.plot(t_min, rc / 1e12, color=PALETTE["copper"], lw=2.0, label=r"Cake Resistance $R_c(t)$")
    ax3.plot(t_min, (rm + rc) / 1e12, color=PALETTE["slate_dark"], lw=2.2, label=r"Total Resistance $R_{\mathrm{total}}(t)$")
    ax3.set_xlabel("Filtration Time (min)")
    ax3.set_ylabel(r"Hydraulic Resistance $(10^{12}\ \mathrm{m}^{-1})$")
    ax3.set_title("C. Resistance Breakdown", loc="left", color=PALETTE["slate_dark"])
    ax3.legend(loc="upper left", frameon=True)

    # Panel D: 2D Space-Time Concentration Boundary Layer Heatmap
    ax4 = axes[1, 1]
    delta_bl_um = 35.0
    y_coords_um = np.linspace(0.0, delta_bl_um, 40)
    c_field_2d = np.zeros((len(y_coords_um), len(cw)))

    for t_idx, cw_val in enumerate(cw):
        c_field_2d[:, t_idx] = cb + (cw_val - cb) * np.exp(-y_coords_um / (delta_bl_um * 0.35))

    im = ax4.imshow(
        c_field_2d,
        aspect="auto",
        origin="lower",
        extent=[t_min[0], t_min[-1], 0.0, delta_bl_um],
        cmap="cividis",
    )
    cbar = fig.colorbar(im, ax=ax4, pad=0.03)
    cbar.set_label(r"Local Concentration $(\mathrm{g} \cdot \mathrm{L}^{-1})$", color=PALETTE["slate_dark"])
    ax4.set_xlabel("Filtration Time (min)")
    ax4.set_ylabel(r"Distance from Membrane $y\ (\mu\mathrm{m})$")
    ax4.set_title("D. Boundary Layer Concentration Field $C(y, t)$", loc="left", color=PALETTE["slate_dark"])

    plt.tight_layout()
    out_p = Path(output_path)
    out_p.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(out_p, bbox_inches="tight")
    plt.close()
    print(f"Saved minimalist 4-panel publication figure to: {out_p}")


def plot_limiting_flux_curves(
    sim: Any,
    output_path: Path | str,
) -> None:
    apply_minimalist_theme()
    import scipy.optimize as opt

    tmps_kpa = np.linspace(20.0, 400.0, 50)
    shears = [2000.0, 4000.0, 8000.0, 12000.0]
    minimal_line_colors = [PALETTE["slate_light"], PALETTE["blue"], PALETTE["teal"], PALETTE["slate_dark"]]

    fig, ax = plt.subplots(figsize=(8.5, 5.5), dpi=300)

    for gamma, col in zip(shears, minimal_line_colors):
        fluxes = []
        for p in tmps_kpa:
            tmp_pa = p * 1000.0
            km = sim.mass_transfer_coefficient(gamma)
            bulk_c = 5.0
            c_gel = sim.props.gel_concentration_g_L
            mu_0 = sim.water_viscosity(sim.thermo.temperature_K)
            rm = 1.0e12
            rc_spec = sim.specific_cake_resistance(tmp_pa)

            def residual(j_val: float) -> float:
                pe = np.clip(j_val / km, 0.0, 6.0)
                cw = min(c_gel, float(bulk_c * np.exp(pe)))
                j_crit = 0.25 * km
                excess = max(0.0, j_val - j_crit)
                sticking = 0.001 * ((cw / c_gel) ** 1.5)
                k_shear = (gamma / 4000.0) * 1.0e-3
                m_cake = (sticking * excess * bulk_c) / max(1e-6, k_shear)
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
                j_sol = km * np.log(c_gel / bulk_c)
            fluxes.append(j_sol * 3.6e6)

        ax.plot(tmps_kpa, fluxes, color=col, lw=2.2, label=rf"$\dot{{\gamma}} = {gamma:.0f}\ \mathrm{{s}}^{{-1}}$")

    ax.set_xlabel(r"Transmembrane Pressure TMP $(\mathrm{kPa})$")
    ax.set_ylabel(r"Steady Permeate Flux $(\mathrm{L} \cdot \mathrm{m}^{-2} \cdot \mathrm{h}^{-1})$")
    ax.set_title("Parametric Limiting Flux Profiles across Shear Rates", loc="left", color=PALETTE["slate_dark"])
    ax.legend(frameon=True, loc="lower right")

    plt.tight_layout()
    out_p = Path(output_path)
    out_p.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(out_p, bbox_inches="tight")
    plt.close()
    print(f"Saved minimalist limiting flux curves to: {out_p}")
