"""Constant-Volume Diafiltration (CVD) and Buffer Exchange Module."""

from __future__ import annotations
from pathlib import Path
from typing import Any
import numpy as np
import matplotlib.pyplot as plt


def simulate_diafiltration(
    v_retentate_L: float = 10.0,
    membrane_area_m2: float = 0.5,
    avg_flux_lmh: float = 40.0,
    r_obs_product: float = 0.998,
    r_obs_impurity: float = 0.05,
    target_dv: float = 8.0,
    output_png: Path | str = "data/diafiltration_summary.png",
) -> dict[str, Any]:
    """Simulates multi-diavolume constant-volume diafiltration for bioparticle purification."""
    dvs = np.linspace(0.0, target_dv, 100)

    # Volumetric permeate flow rate (L/h)
    q_permeate_l_h = avg_flux_lmh * membrane_area_m2
    # Time per DV = V_retentate / Q_permeate (hours)
    time_hours = dvs * (v_retentate_L / q_permeate_l_h)

    # Impurity clearance & Product recovery
    c_imp_norm = np.exp(-(1.0 - r_obs_impurity) * dvs) * 100.0  # % of initial
    product_yield = np.exp(-(1.0 - r_obs_product) * dvs) * 100.0  # % recovery

    buffer_consumed_L = dvs * v_retentate_L

    fig, axes = plt.subplots(1, 3, figsize=(16, 5), dpi=300)
    fig.suptitle("Constant-Volume Diafiltration & Buffer Exchange Kinetics", fontsize=14, fontweight="bold")

    # Panel 1: Impurity Clearance & Product Yield
    ax1 = axes[0]
    ax1.plot(dvs, c_imp_norm, "r-", lw=2.5, label=r"Impurity Concentration $C_{\mathrm{imp}}(DV)$")
    ax1.plot(dvs, product_yield, "b-", lw=2.5, label=r"Bioparticle Yield $Y_{\mathrm{product}}(DV)$")
    ax1.axvline(x=5.0, color="gray", linestyle="--", alpha=0.7, label=r"$5\times\mathrm{DV}$ Target ($99.3\%$ Cleared)")
    ax1.set_xlabel("Diavolumes (DV)", fontsize=11)
    ax1.set_ylabel("Normalized Percentage (%)", fontsize=11)
    ax1.set_title("A. Buffer Exchange & Clearance Kinetics", fontsize=12, fontweight="bold", loc="left")
    ax1.grid(True, linestyle="--", alpha=0.5)
    ax1.legend(loc="center right", frameon=True)

    # Panel 2: Process Time Evolution
    ax2 = axes[1]
    ax2.plot(dvs, time_hours * 60.0, color="#2ca02c", lw=2.5)
    ax2.set_xlabel("Diavolumes (DV)", fontsize=11)
    ax2.set_ylabel("Cumulative Process Time (min)", fontsize=11)
    ax2.set_title("B. Diafiltration Time Scaling", fontsize=12, fontweight="bold", loc="left")
    ax2.grid(True, linestyle="--", alpha=0.5)

    # Panel 3: Buffer Volume Requirement
    ax3 = axes[2]
    ax3.plot(dvs, buffer_consumed_L, color="#9467bd", lw=2.5)
    ax3.set_xlabel("Diavolumes (DV)", fontsize=11)
    ax3.set_ylabel("Buffer Consumed (L)", fontsize=11)
    ax3.set_title("C. Buffer Consumption ($V_{\mathrm{buffer}} = \mathrm{DV} \times V_{\mathrm{ret}}$)", fontsize=12, fontweight="bold", loc="left")
    ax3.grid(True, linestyle="--", alpha=0.5)

    plt.tight_layout()
    out_p = Path(output_png)
    out_p.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(out_p, bbox_inches="tight")
    plt.close()
    print(f"Successfully generated diafiltration dashboard: {out_p}")

    return {
        "final_dv": target_dv,
        "final_product_yield_pct": float(product_yield[-1]),
        "final_impurity_remaining_pct": float(c_imp_norm[-1]),
        "total_process_time_hours": float(time_hours[-1]),
        "total_buffer_consumed_L": float(buffer_consumed_L[-1]),
    }


if __name__ == "__main__":
    simulate_diafiltration()
