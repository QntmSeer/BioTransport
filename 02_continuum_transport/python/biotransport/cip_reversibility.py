"""Cleaning-in-Place (CIP) and Fouling Reversibility Diagnostics (Minimalist Styling)."""

from __future__ import annotations
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from biotransport.theme import apply_minimalist_theme, PALETTE


def simulate_cip_reversibility(
    r_m: float = 1.2e12,
    r_rev_cake: float = 4.8e12,
    r_irr_adsorp: float = 8.5e11,
    output_png: Path | str = "data/cip_fouling_reversibility.png",
) -> dict[str, float]:
    apply_minimalist_theme()
    r_total = r_m + r_rev_cake + r_irr_adsorp

    pct_rm = (r_m / r_total) * 100.0
    pct_rev = (r_rev_cake / r_total) * 100.0
    pct_irr = (r_irr_adsorp / r_total) * 100.0

    cycles = np.arange(1, 6)
    flux_recovery_pct = [100.0, 99.2, 98.5, 97.9, 97.3]
    fouled_end_flux_pct = [22.4, 21.8, 21.2, 20.6, 20.1]

    fig, axes = plt.subplots(1, 2, figsize=(13.5, 5.2), dpi=300)
    fig.suptitle("Membrane Cleaning-in-Place (CIP) & Fouling Reversibility Analysis", fontsize=14, color=PALETTE["slate_dark"])

    # Panel 1: Hydraulic Resistance Breakdown
    ax1 = axes[0]
    labels = [
        f"Membrane $R_m$\n({pct_rm:.1f}%)",
        f"Reversible Cake $R_{{rev}}$\n({pct_rev:.1f}%)",
        f"Pore Adsorption $R_{{irr}}$\n({pct_irr:.1f}%)",
    ]
    colors = [PALETTE["slate_light"], PALETTE["teal"], PALETTE["copper"]]
    explode = (0.02, 0.03, 0.03)
    ax1.pie(
        [r_m, r_rev_cake, r_irr_adsorp],
        labels=labels,
        colors=colors,
        explode=explode,
        autopct="%1.1f%%",
        startangle=140,
        textprops={"fontsize": 9.5, "color": PALETTE["slate_dark"], "fontweight": "500"},
        wedgeprops={"edgecolor": "#ffffff", "linewidth": 1.2},
    )
    ax1.set_title("A. Total Hydraulic Resistance Distribution", color=PALETTE["slate_dark"])

    # Panel 2: Multi-Cycle CIP Flux Recovery
    ax2 = axes[1]
    width = 0.32
    ax2.bar(cycles - width / 2, flux_recovery_pct, width, label="Post-CIP Recovery (0.1M NaOH)", color=PALETTE["teal"], edgecolor=PALETTE["slate_dark"], lw=0.8)
    ax2.bar(cycles + width / 2, fouled_end_flux_pct, width, label="End-of-Run Fouled Flux", color=PALETTE["slate_light"], edgecolor=PALETTE["slate_dark"], lw=0.8)

    for i, fr in enumerate(flux_recovery_pct):
        ax2.text(cycles[i] - width / 2, fr + 1.5, f"{fr:.1f}%", ha="center", fontsize=8.5, color=PALETTE["slate_dark"], fontweight="600")

    ax2.set_xlabel("Reusability Batch Cycle")
    ax2.set_ylabel("Normalized Permeate Flux (%)")
    ax2.set_ylim([0, 115])
    ax2.set_title(r"B. 5-Cycle CIP Flux Recovery ($FR \geq 97.3\%$)", loc="left", color=PALETTE["slate_dark"])
    ax2.legend(loc="lower left", frameon=True)

    plt.tight_layout()
    out_p = Path(output_png)
    out_p.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(out_p, bbox_inches="tight")
    plt.close()
    print(f"Successfully generated minimalist CIP reversibility summary: {out_p}")

    return {
        "reversible_fouling_pct": float(pct_rev),
        "irreversible_fouling_pct": float(pct_irr),
        "cycle_5_flux_recovery_pct": float(flux_recovery_pct[-1]),
    }


if __name__ == "__main__":
    simulate_cip_reversibility()
