"""Constant-Volume Diafiltration (CVD) and Buffer Exchange Module (Minimalist Styling)."""

from __future__ import annotations
from pathlib import Path
from typing import Any
import numpy as np
import matplotlib.pyplot as plt
from biotransport.theme import apply_minimalist_theme, PALETTE


def simulate_diafiltration(
    v_retentate_L: float = 10.0,
    membrane_area_m2: float = 0.5,
    avg_flux_lmh: float = 40.0,
    r_obs_product: float = 0.998,
    r_obs_impurity: float = 0.05,
    target_dv: float = 8.0,
    output_png: Path | str = "data/diafiltration_summary.png",
) -> dict[str, Any]:
    apply_minimalist_theme()
    dvs = np.linspace(0.0, target_dv, 100)

    q_permeate_l_h = avg_flux_lmh * membrane_area_m2
    time_hours = dvs * (v_retentate_L / q_permeate_l_h)

    c_imp_norm = np.exp(-(1.0 - r_obs_impurity) * dvs) * 100.0
    product_yield = np.exp(-(1.0 - r_obs_product) * dvs) * 100.0
    buffer_consumed_L = dvs * v_retentate_L

    fig, axes = plt.subplots(1, 3, figsize=(15.5, 4.8), dpi=300)
    fig.suptitle("Constant-Volume Diafiltration & Buffer Exchange Kinetics", fontsize=14, color=PALETTE["slate_dark"])

    # Panel 1: Impurity Clearance & Product Yield
    ax1 = axes[0]
    ax1.plot(dvs, c_imp_norm, color=PALETTE["copper"], lw=2.2, label=r"Impurity $C_{\mathrm{imp}}(DV)$")
    ax1.plot(dvs, product_yield, color=PALETTE["teal"], lw=2.2, label=r"Bioparticle Yield $Y(DV)$")
    ax1.axvline(x=5.0, color=PALETTE["slate_light"], linestyle=":", lw=1.5, label=r"$5\times\mathrm{DV}$ Target ($99.3\%$ Cleared)")
    ax1.set_xlabel("Diavolumes (DV)")
    ax1.set_ylabel("Normalized Percentage (%)")
    ax1.set_title("A. Buffer Exchange & Clearance Kinetics", loc="left", color=PALETTE["slate_dark"])
    ax1.legend(loc="center right", frameon=True)

    # Panel 2: Process Time Evolution
    ax2 = axes[1]
    ax2.plot(dvs, time_hours * 60.0, color=PALETTE["slate_dark"], lw=2.2)
    ax2.set_xlabel("Diavolumes (DV)")
    ax2.set_ylabel("Cumulative Process Time (min)")
    ax2.set_title("B. Diafiltration Time Scaling", loc="left", color=PALETTE["slate_dark"])

    # Panel 3: Buffer Volume Requirement
    ax3 = axes[2]
    ax3.plot(dvs, buffer_consumed_L, color=PALETTE["blue"], lw=2.2)
    ax3.set_xlabel("Diavolumes (DV)")
    ax3.set_ylabel("Buffer Consumed (L)")
    ax3.set_title(r"C. Buffer Consumption ($V_{\mathrm{buf}} = \mathrm{DV} \times V_{\mathrm{ret}}$)", loc="left", color=PALETTE["slate_dark"])

    plt.tight_layout()
    out_p = Path(output_png)
    out_p.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(out_p, bbox_inches="tight")
    plt.close()
    print(f"Successfully generated minimalist diafiltration dashboard: {out_p}")

    return {
        "final_dv": target_dv,
        "final_product_yield_pct": float(product_yield[-1]),
        "final_impurity_remaining_pct": float(c_imp_norm[-1]),
        "total_process_time_hours": float(time_hours[-1]),
        "total_buffer_consumed_L": float(buffer_consumed_L[-1]),
    }


if __name__ == "__main__":
    simulate_diafiltration()
