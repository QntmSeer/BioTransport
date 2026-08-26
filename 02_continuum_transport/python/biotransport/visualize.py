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
    """Plots 3-panel figure: Permeate Flux J(t), Wall Concentration C_w(t), and Cake Resistance R_c(t)."""
    t_min = np.array(results["time_s"]) / 60.0  # Convert seconds to minutes

    fig, axes = plt.subplots(1, 3, figsize=(15, 4.5), dpi=300)
    fig.suptitle(title, fontsize=14, fontweight="bold", y=1.02)

    # Panel 1: Permeate Flux Decline
    ax1 = axes[0]
    ax1.plot(t_min, results["flux_lmh"], color="#1f77b4", lw=2.5, label="Permeate Flux $J(t)$")
    ax1.set_xlabel("Filtration Time (min)", fontsize=11)
    ax1.set_ylabel("Permeate Flux (LMH)", fontsize=11)
    ax1.set_title("Permeate Flux Decline", fontsize=12)
    ax1.grid(True, linestyle="--", alpha=0.6)
    ax1.legend(frameon=True)

    # Panel 2: Concentration Polarization (Wall Concentration)
    ax2 = axes[1]
    ax2.plot(t_min, results["wall_conc_g_l"], color="#d62728", lw=2.5, label="Wall Conc $C_w(t)$")
    ax2.set_xlabel("Filtration Time (min)", fontsize=11)
    ax2.set_ylabel("Membrane Wall Concentration (g/L)", fontsize=11)
    ax2.set_title("Concentration Polarization", fontsize=12)
    ax2.grid(True, linestyle="--", alpha=0.6)
    ax2.legend(frameon=True)

    # Panel 3: Cake Resistance Growth
    ax3 = axes[2]
    ax3.plot(t_min, np.array(results["cake_resistance_m_inv"]) / 1e12, color="#2ca02c", lw=2.5, label="Cake Resistance $R_c(t)$")
    ax3.set_xlabel("Filtration Time (min)", fontsize=11)
    ax3.set_ylabel("Cake Resistance ($10^{12} \\mathrm{m}^{-1}$)", fontsize=11)
    ax3.set_title("Fouling Cake Layer Growth", fontsize=12)
    ax3.grid(True, linestyle="--", alpha=0.6)
    ax3.legend(frameon=True)

    plt.tight_layout()
    if output_png:
        out_p = Path(output_png)
        out_p.parent.mkdir(parents=True, exist_ok=True)
        plt.savefig(out_p, bbox_inches="tight")
        print(f"Saved publication figure to: {out_p}")
    plt.close()
