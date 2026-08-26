"""Hermia's 4-Mechanism Fouling Diagnostic Engine."""

from __future__ import annotations
from pathlib import Path
from typing import Any
import numpy as np
import matplotlib.pyplot as plt


def analyze_hermia_fouling(
    results: dict[str, Any],
    output_png: Path | str = "data/hermia_fouling_analysis.png",
) -> dict[str, float]:
    """Fits the filtration time-series data against Hermia's 4 fouling models:
    1. Complete Blocking (n = 2.0)
    2. Standard Pore Constriction (n = 1.5)
    3. Intermediate Blocking (n = 1.0)
    4. Cake Filtration (n = 0.0)
    """
    t_s = np.array(results["time_s"])
    j_lmh = np.array(results["flux_lmh"])
    j_norm = j_lmh / max(1e-6, j_lmh[0])

    # Avoid zero division
    t_s = t_s[1:]
    j_norm = j_norm[1:]

    # 1. Complete Blocking: ln(J/J0) vs t
    y_complete = np.log(np.clip(j_norm, 1e-6, 1.0))
    r_complete = np.corrcoef(t_s, y_complete)[0, 1] ** 2

    # 2. Standard Pore Constriction: (J/J0)^(-0.5) vs t
    y_standard = (j_norm) ** (-0.5)
    r_standard = np.corrcoef(t_s, y_standard)[0, 1] ** 2

    # 3. Intermediate Blocking: (J/J0)^(-1) vs t
    y_intermediate = (j_norm) ** (-1.0)
    r_intermediate = np.corrcoef(t_s, y_intermediate)[0, 1] ** 2

    # 4. Cake Filtration: (J/J0)^(-2) vs t
    y_cake = (j_norm) ** (-2.0)
    r_cake = np.corrcoef(t_s, y_cake)[0, 1] ** 2

    scores = {
        "Complete Pore Blocking (n=2.0)": float(r_complete),
        "Standard Pore Constriction (n=1.5)": float(r_standard),
        "Intermediate Blocking (n=1.0)": float(r_intermediate),
        "Cake Filtration (n=0.0)": float(r_cake),
    }

    fig, axes = plt.subplots(2, 2, figsize=(12, 9), dpi=300)
    fig.suptitle("Hermia's 4-Mechanism Fouling Diagnostic Breakdown", fontsize=14, fontweight="bold")

    t_min = t_s / 60.0

    # Panel 1: Complete Blocking
    axes[0, 0].plot(t_min, y_complete, "o-", color="#1f77b4", markersize=3, lw=1.8)
    axes[0, 0].set_title(f"A. Complete Blocking ($R^2 = {r_complete:.4f}$)", fontsize=11, fontweight="bold")
    axes[0, 0].set_xlabel("Filtration Time (min)")
    axes[0, 0].set_ylabel(r"$\ln(J / J_0)$")
    axes[0, 0].grid(True, linestyle="--", alpha=0.5)

    # Panel 2: Standard Constriction
    axes[0, 1].plot(t_min, y_standard, "s-", color="#ff7f0e", markersize=3, lw=1.8)
    axes[0, 1].set_title(f"B. Standard Pore Constriction ($R^2 = {r_standard:.4f}$)", fontsize=11, fontweight="bold")
    axes[0, 1].set_xlabel("Filtration Time (min)")
    axes[0, 1].set_ylabel(r"$(J / J_0)^{-1/2}$")
    axes[0, 1].grid(True, linestyle="--", alpha=0.5)

    # Panel 3: Intermediate Blocking
    axes[1, 0].plot(t_min, y_intermediate, "^-", color="#2ca02c", markersize=3, lw=1.8)
    axes[1, 0].set_title(f"C. Intermediate Blocking ($R^2 = {r_intermediate:.4f}$)", fontsize=11, fontweight="bold")
    axes[1, 0].set_xlabel("Filtration Time (min)")
    axes[1, 0].set_ylabel(r"$(J / J_0)^{-1}$")
    axes[1, 0].grid(True, linestyle="--", alpha=0.5)

    # Panel 4: Cake Formation
    axes[1, 1].plot(t_min, y_cake, "d-", color="#d62728", markersize=3, lw=1.8)
    axes[1, 1].set_title(f"D. Cake Filtration ($R^2 = {r_cake:.4f}$)", fontsize=11, fontweight="bold")
    axes[1, 1].set_xlabel("Filtration Time (min)")
    axes[1, 1].set_ylabel(r"$(J / J_0)^{-2}$")
    axes[1, 1].grid(True, linestyle="--", alpha=0.5)

    plt.tight_layout()
    out_p = Path(output_png)
    out_p.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(out_p, bbox_inches="tight")
    plt.close()
    print(f"Successfully generated Hermia diagnostic: {out_p}")

    return scores


if __name__ == "__main__":
    from biotransport.bridge import MdBridgeModel, ProcessSimulator

    m = MdBridgeModel.load_json("data/sample_md_params.json")
    sim = ProcessSimulator(m)
    res = sim.simulate_filtration(tmp_pa=200_000.0, total_time_s=3600.0)
    analyze_hermia_fouling(res)
