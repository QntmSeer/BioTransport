"""Techno-Economic Analysis (TEA) and Membrane Sizing Module for Bioparticle Purification."""

from __future__ import annotations
from pathlib import Path
from typing import Any
import numpy as np
import matplotlib.pyplot as plt


def calculate_techno_economics(
    batch_volume_L: float = 100.0,
    concentration_factor: float = 10.0,
    target_time_hours: float = 3.0,
    avg_flux_lmh: float = 45.0,
    membrane_cost_per_m2: float = 850.0,  # USD / m^2 PES membrane
    electricity_cost_per_kwh: float = 0.15,
    cassette_area_m2: float = 0.5,
    output_png: Path | str = "data/techno_economics_summary.png",
) -> dict[str, Any]:
    """Calculates industrial scale-up membrane sizing, cassette requirements, and operational costs."""
    v_permeate_L = batch_volume_L * (1.0 - 1.0 / concentration_factor)
    req_membrane_area_m2 = v_permeate_L / (avg_flux_lmh * target_time_hours)
    num_cassettes = int(np.ceil(req_membrane_area_m2 / cassette_area_m2))
    installed_area_m2 = num_cassettes * cassette_area_m2

    # Pumping power: Q_feed * TMP / pump_eff
    # Crossflow flux typically 5-10x permeate rate
    q_feed_m3_s = (5.0 * avg_flux_lmh * installed_area_m2 / 1000.0) / 3600.0
    tmp_pa = 200_000.0
    pump_eff = 0.65
    pump_power_kw = (q_feed_m3_s * tmp_pa / pump_eff) / 1000.0
    energy_kwh = pump_power_kw * target_time_hours

    capex_membrane_usd = installed_area_m2 * membrane_cost_per_m2
    opex_energy_usd = energy_kwh * electricity_cost_per_kwh
    buffer_cost_usd = (v_permeate_L * 0.05)  # $0.05 / L buffer

    # Batch scaling curves (10 L to 1000 L)
    batch_sizes = np.linspace(10.0, 1000.0, 50)
    area_curve = (batch_sizes * (1.0 - 1.0 / concentration_factor)) / (avg_flux_lmh * target_time_hours)
    energy_curve = ((5.0 * avg_flux_lmh * area_curve / 1000.0) / 3600.0 * tmp_pa / pump_eff / 1000.0) * target_time_hours

    fig, axes = plt.subplots(1, 2, figsize=(14, 5.5), dpi=300)
    fig.suptitle("Bioprocess Techno-Economics & Membrane Scale-Up Analysis", fontsize=14, fontweight="bold")

    # Panel 1: Membrane Area Scaling vs Batch Size
    ax1 = axes[0]
    ax1.plot(batch_sizes, area_curve, color="#1f77b4", lw=2.5, label="Required Membrane Area")
    ax1.scatter([batch_volume_L], [req_membrane_area_m2], color="crimson", s=100, zorder=5, label=f"Baseline {batch_volume_L:.0f}L Batch ($A = {req_membrane_area_m2:.2f}\\,\\mathrm{{m}}^2$)")
    ax1.set_xlabel("Batch Volume (L)", fontsize=11)
    ax1.set_ylabel(r"Required Membrane Area $(\mathrm{m}^2)$", fontsize=11)
    ax1.set_title("A. Membrane Sizing Scale-Up Curve", fontsize=12, fontweight="bold", loc="left")
    ax1.grid(True, linestyle="--", alpha=0.5)
    ax1.legend(loc="upper left", frameon=True)

    # Panel 2: Total Operating Energy vs Batch Size
    ax2 = axes[1]
    ax2.plot(batch_sizes, energy_curve, color="#ff7f0e", lw=2.5, label="Pumping Energy Consumption")
    ax2.scatter([batch_volume_L], [energy_kwh], color="crimson", s=100, zorder=5, label=f"Baseline Energy ($E = {energy_kwh:.2f}\\,\\mathrm{{kWh}}$)")
    ax2.set_xlabel("Batch Volume (L)", fontsize=11)
    ax2.set_ylabel("Pumping Energy (kWh)", fontsize=11)
    ax2.set_title("B. Energy Consumption Scaling", fontsize=12, fontweight="bold", loc="left")
    ax2.grid(True, linestyle="--", alpha=0.5)
    ax2.legend(loc="upper left", frameon=True)

    plt.tight_layout()
    out_p = Path(output_png)
    out_p.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(out_p, bbox_inches="tight")
    plt.close()
    print(f"Successfully generated techno-economics summary: {out_p}")

    return {
        "batch_volume_L": batch_volume_L,
        "required_membrane_area_m2": float(req_membrane_area_m2),
        "installed_cassettes": num_cassettes,
        "installed_area_m2": float(installed_area_m2),
        "energy_consumed_kwh": float(energy_kwh),
        "capex_membrane_usd": float(capex_membrane_usd),
        "opex_energy_usd": float(opex_energy_usd),
    }


if __name__ == "__main__":
    calculate_techno_economics()
