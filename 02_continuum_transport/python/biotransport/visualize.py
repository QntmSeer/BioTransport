"""Publication-grade visualization tools for multiscale bioparticle transport."""

from __future__ import annotations
from pathlib import Path
from typing import Any
import matplotlib.pyplot as plt
import numpy as np


def plot_multiscale_results(
    results: dict[str, Any],
    output_png: Path | str | None = None,
    title: str = "Multiscale Bioparticle Filtration & Fouling Dynamics",
) -> None:
    """Plots a comprehensive 4-panel publication figure:
    1. Permeate Flux J(t)
    2. Concentration Polarization C_w(t)
    3. Hydraulic Resistance Breakdown (R_membrane vs R_cake)
    4. 2D Boundary Layer Concentration Heatmap C(y, t)
    """
    t_min = np.array(results["time_s"]) / 60.0
    flux = np.array(results["flux_lmh"])
    cw = np.array(results["wall_conc_g_l"])
    rc = np.array(results["cake_resistance_m_inv"]) / 1e12  # in 10^12 m^-1
    rm = 1.0  # 1.0 * 10^12 m^-1 baseline

    fig, axes = plt.subplots(2, 2, figsize=(14, 10), dpi=300)
    fig.suptitle(title, fontsize=15, fontweight="bold", y=0.98)

    # Panel 1: Permeate Flux Decline
    ax1 = axes[0, 0]
    ax1.plot(t_min, flux, color="#1f77b4", lw=2.5, label=r"Permeate Flux $J(t)$")
    ax1.fill_between(t_min, flux, color="#1f77b4", alpha=0.15)
    ax1.set_xlabel("Filtration Time (min)", fontsize=11)
    ax1.set_ylabel(r"Permeate Flux $(\mathrm{L} \cdot \mathrm{m}^{-2} \cdot \mathrm{h}^{-1})$", fontsize=11)
    ax1.set_title("A. Permeate Flux Decline", fontsize=12, fontweight="bold", loc="left")
    ax1.grid(True, linestyle="--", alpha=0.5)
    ax1.legend(frameon=True, loc="upper right")

    # Panel 2: Concentration Polarization
    ax2 = axes[0, 1]
    ax2.plot(t_min, cw, color="#d62728", lw=2.5, label=r"Wall Concentration $C_w(t)$")
    ax2.axhline(y=400.0, color="darkred", linestyle=":", lw=2.0, label=r"Critical Gel Limit $C_{\mathrm{gel}}$")
    ax2.set_xlabel("Filtration Time (min)", fontsize=11)
    ax2.set_ylabel(r"Membrane Wall Concentration $(\mathrm{g} \cdot \mathrm{L}^{-1})$", fontsize=11)
    ax2.set_title("B. Membrane Boundary Layer Polarization", fontsize=12, fontweight="bold", loc="left")
    ax2.grid(True, linestyle="--", alpha=0.5)
    ax2.legend(frameon=True, loc="lower right")

    # Panel 3: Hydraulic Resistance Breakdown
    ax3 = axes[1, 0]
    ax3.stackplot(
        t_min,
        [np.full_like(t_min, rm), rc],
        labels=[r"Clean Membrane $R_m$", r"Fouling Cake $R_c(t)$"],
        colors=["#aec7e8", "#ff9896"],
        alpha=0.85,
    )
    ax3.set_xlabel("Filtration Time (min)", fontsize=11)
    ax3.set_ylabel(r"Hydraulic Resistance $(10^{12}\ \mathrm{m}^{-1})$", fontsize=11)
    ax3.set_title("C. Resistance-in-Series Dynamics", fontsize=12, fontweight="bold", loc="left")
    ax3.grid(True, linestyle="--", alpha=0.5)
    ax3.legend(frameon=True, loc="upper left")

    # Panel 4: 2D Boundary Layer Space-Time Profile C(y, t)
    ax4 = axes[1, 1]
    # Reconstruct 2D field C(y, t) across 30 spatial grid points across boundary layer delta ~ 15 um
    n_spatial = 40
    y_um = np.linspace(0.0, 20.0, n_spatial)  # 0 um at membrane surface to 20 um in bulk
    c_2d = np.zeros((n_spatial, len(t_min)))

    bulk_c = 10.0
    for j, cw_val in enumerate(cw):
        # Exponential boundary decay: C(y) = C_b + (C_w - C_b) * exp(-y / delta_eff)
        decay_length = 5.0  # um
        c_2d[:, j] = bulk_c + (cw_val - bulk_c) * np.exp(-y_um / decay_length)

    extent = [t_min[0], t_min[-1], y_um[-1], y_um[0]]  # y=0 at bottom/membrane
    im = ax4.imshow(c_2d, aspect="auto", cmap="viridis", extent=extent, interpolation="bicubic")
    cbar = fig.colorbar(im, ax=ax4, pad=0.03)
    cbar.set_label(r"Local Concentration $(\mathrm{g} \cdot \mathrm{L}^{-1})$", fontsize=10)
    ax4.set_xlabel("Filtration Time (min)", fontsize=11)
    ax4.set_ylabel(r"Distance from Membrane $y\ (\mu\mathrm{m})$", fontsize=11)
    ax4.set_title("D. 2D Boundary Layer Space-Time Field $C(y, t)$", fontsize=12, fontweight="bold", loc="left")

    plt.tight_layout()
    if output_png:
        out_p = Path(output_png)
        out_p.parent.mkdir(parents=True, exist_ok=True)
        plt.savefig(out_p, bbox_inches="tight")
        print(f"Saved 4-panel publication figure to: {out_p}")
    plt.close()


def plot_limiting_flux_curves(
    sim_instance: Any,
    output_png: Path | str | None = None,
) -> None:
    """Plots parametric Limiting Flux vs. TMP curves across varying crossflow shear rates."""
    tmps = np.linspace(20_000, 350_000, 15)  # 0.2 to 3.5 bar
    tmps_bar = tmps / 100_000.0
    shear_rates = [2000.0, 4000.0, 8000.0, 12000.0]
    colors = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728"]

    fig, ax = plt.subplots(figsize=(8, 5.5), dpi=300)
    for gamma, color in zip(shear_rates, colors):
        fluxes = []
        for tmp in tmps:
            res = sim_instance.simulate_filtration(
                tmp_pa=tmp,
                shear_rate_s_inv=gamma,
                total_time_s=1800.0,
                n_steps=10,
            )
            # Use pseudo-steady state flux
            fluxes.append(res["flux_lmh"][-1])

        ax.plot(tmps_bar, fluxes, marker="o", color=color, lw=2.2, label=f"$\\dot{{\\gamma}} = {int(gamma)}\\ \\mathrm{{s}}^{{-1}}$")

    ax.set_xlabel("Transmembrane Pressure (bar)", fontsize=12)
    ax.set_ylabel(r"Steady Permeate Flux $(\mathrm{L} \cdot \mathrm{m}^{-2} \cdot \mathrm{h}^{-1})$", fontsize=12)
    ax.set_title("Parametric Limiting Flux Curves (Gel-Polarization & Cake Limiting)", fontsize=13, fontweight="bold")
    ax.grid(True, linestyle="--", alpha=0.6)
    ax.legend(title="Crossflow Shear Rate", frameon=True)

    if output_png:
        out_p = Path(output_png)
        out_p.parent.mkdir(parents=True, exist_ok=True)
        plt.savefig(out_p, bbox_inches="tight")
        print(f"Saved limiting flux curves to: {out_p}")
    plt.close()
