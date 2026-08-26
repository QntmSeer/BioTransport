"""Hermia's 4-Mechanism Fouling Diagnostic Engine (Minimalist Styling)."""

from __future__ import annotations
from pathlib import Path
from typing import Any
import numpy as np
import matplotlib.pyplot as plt
from biotransport.theme import apply_minimalist_theme, PALETTE


def analyze_hermia_fouling(
    results: dict[str, Any],
    output_png: Path | str = "data/hermia_fouling_analysis.png",
) -> dict[str, float]:
    apply_minimalist_theme()
    t_s = np.array(results["time_s"])
    j_lmh = np.array(results["flux_lmh"])
    j_norm = j_lmh / max(1e-6, j_lmh[0])

    t_s = t_s[1:]
    j_norm = j_norm[1:]

    y_complete = np.log(np.clip(j_norm, 1e-6, 1.0))
    r_complete = np.corrcoef(t_s, y_complete)[0, 1] ** 2

    y_standard = (j_norm) ** (-0.5)
    r_standard = np.corrcoef(t_s, y_standard)[0, 1] ** 2

    y_intermediate = (j_norm) ** (-1.0)
    r_intermediate = np.corrcoef(t_s, y_intermediate)[0, 1] ** 2

    y_cake = (j_norm) ** (-2.0)
    r_cake = np.corrcoef(t_s, y_cake)[0, 1] ** 2

    scores = {
        "Complete Pore Blocking (n=2.0)": float(r_complete),
        "Standard Pore Constriction (n=1.5)": float(r_standard),
        "Intermediate Blocking (n=1.0)": float(r_intermediate),
        "Cake Filtration (n=0.0)": float(r_cake),
    }

    fig, axes = plt.subplots(2, 2, figsize=(11.5, 8.5), dpi=300)
    fig.suptitle("Hermia's 4-Mechanism Fouling Diagnostic Breakdown", fontsize=14, color=PALETTE["slate_dark"])

    t_min = t_s / 60.0

    # Panel 1: Complete Blocking
    axes[0, 0].plot(t_min, y_complete, color=PALETTE["slate_med"], lw=1.8, marker="o", markersize=2.5)
    axes[0, 0].set_title(f"A. Complete Blocking ($R^2 = {r_complete:.4f}$)", loc="left", color=PALETTE["slate_dark"])
    axes[0, 0].set_xlabel("Filtration Time (min)")
    axes[0, 0].set_ylabel(r"$\ln(J / J_0)$")

    # Panel 2: Standard Constriction
    axes[0, 1].plot(t_min, y_standard, color=PALETTE["teal"], lw=1.8, marker="s", markersize=2.5)
    axes[0, 1].set_title(f"B. Standard Pore Constriction ($R^2 = {r_standard:.4f}$)", loc="left", color=PALETTE["slate_dark"])
    axes[0, 1].set_xlabel("Filtration Time (min)")
    axes[0, 1].set_ylabel(r"$(J / J_0)^{-1/2}$")

    # Panel 3: Intermediate Blocking
    axes[1, 0].plot(t_min, y_intermediate, color=PALETTE["copper"], lw=1.8, marker="^", markersize=2.5)
    axes[1, 0].set_title(f"C. Intermediate Blocking ($R^2 = {r_intermediate:.4f}$)", loc="left", color=PALETTE["slate_dark"])
    axes[1, 0].set_xlabel("Filtration Time (min)")
    axes[1, 0].set_ylabel(r"$(J / J_0)^{-1}$")

    # Panel 4: Cake Formation
    axes[1, 1].plot(t_min, y_cake, color=PALETTE["blue"], lw=1.8, marker="d", markersize=2.5)
    axes[1, 1].set_title(f"D. Cake Filtration ($R^2 = {r_cake:.4f}$)", loc="left", color=PALETTE["slate_dark"])
    axes[1, 1].set_xlabel("Filtration Time (min)")
    axes[1, 1].set_ylabel(r"$(J / J_0)^{-2}$")

    plt.tight_layout()
    out_p = Path(output_png)
    out_p.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(out_p, bbox_inches="tight")
    plt.close()
    print(f"Successfully generated minimalist Hermia diagnostic: {out_p}")

    return scores


if __name__ == "__main__":
    from biotransport.bridge import MdBridgeModel, ProcessSimulator

    m = MdBridgeModel.load_json("data/sample_md_params.json")
    sim = ProcessSimulator(m)
    res = sim.simulate_filtration(tmp_pa=200_000.0, total_time_s=3600.0)
    analyze_hermia_fouling(res)
