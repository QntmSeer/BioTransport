"""Techno-Economic Analysis and Membrane Sizing (Minimalist Styling)."""

from __future__ import annotations
from pathlib import Path
from typing import Any
import numpy as np
import matplotlib.pyplot as plt
from biotransport.theme import apply_minimalist_theme, PALETTE


def calculate_techno_economics(
    batch_volume_L: float = 100.0,
    concentration_factor: float = 10.0,
    target_time_hours: float = 3.0,
    avg_flux_lmh: float = 45.0,
    membrane_cost_per_m2: float = 850.0,
    electricity_cost_per_kwh: float = 0.15,
    cassette_area_m2: float = 0.5,
    output_png: Path | str = "data/techno_economics_summary.png",
) -> dict[str, Any]:
    apply_minimalist_theme()
    v_permeate_L = batch_volume_L * (1.0 - 1.0 / concentration_factor)
    req_membrane_area_m2 = v_permeate_L / (avg_flux_lmh * target_time_hours)
    num_cassettes = int(np.ceil(req_membrane_area_m2 / cassette_area_m2))
    installed_area_m2 = num_cassettes * cassette_area_m2

    q_feed_m3_s = (5.0 * avg_flux_lmh * installed_area_m2 / 1000.0) / 3600.0
    tmp_pa = 200_000.0
    pump_eff = 0.65
    pump_power_kw = (q_feed_m3_s * tmp_pa / pump_eff) / 1000.0
    energy_kwh = pump_power_kw * target_time_hours

    capex_membrane_usd = installed_area_m2 * membrane_cost_per_m2
    opex_energy_usd = energy_kwh * electricity_cost_per_kwh

    batch_sizes = np.linspace(10.0, 1000.0, 50)
    area_curve = (batch_sizes * (1.0 - 1.0 / concentration_factor)) / (avg_flux_lmh * target_time_hours)
    energy_curve = ((5.0 * avg_flux_lmh * area_curve / 1000.0) / 3600.0 * tmp_pa / pump_eff / 1000.0) * target_time_hours

    fig, axes = plt.subplots(1, 2, figsize=(13.5, 5.2), dpi=300)
    fig.suptitle("Bioprocess Techno-Economics & Membrane Scale-Up Analysis", fontsize=14, color=PALETTE["slate_dark"])

    # Panel 1: Membrane Area Scaling
    ax1 = axes[0]
    ax1.plot(batch_sizes, area_curve, color=PALETTE["blue"], lw=2.2, label="Required Area")
    ax1.scatter([batch_volume_L], [req_membrane_area_m2], color="#ffffff", edgecolors=PALETTE["copper"], s=100, lw=2.0, zorder=5, label=f"Baseline ({req_membrane_area_m2:.2f}\\,\\mathrm{{m}}^2)")
    ax1.set_xlabel("Batch Volume (L)")
    ax1.set_ylabel(r"Required Membrane Area $(\mathrm{m}^2)$")
    ax1.set_title("A. Membrane Sizing Scale-Up Curve", loc="left", color=PALETTE["slate_dark"])
    ax1.legend(loc="upper left", frameon=True)

    # Panel 2: Total Operating Energy
    ax2 = axes[1]
    ax2.plot(batch_sizes, energy_curve, color=PALETTE["teal"], lw=2.2, label="Pumping Energy")
    ax2.scatter([batch_volume_L], [energy_kwh], color="#ffffff", edgecolors=PALETTE["teal"], s=100, lw=2.0, zorder=5, label=f"Baseline ({energy_kwh:.2f}\\,\\mathrm{{kWh}})")
    ax2.set_xlabel("Batch Volume (L)")
    ax2.set_ylabel("Pumping Energy (kWh)")
    ax2.set_title("B. Energy Consumption Scaling", loc="left", color=PALETTE["slate_dark"])
    ax2.legend(loc="upper left", frameon=True)

    plt.tight_layout()
    out_p = Path(output_png)
    out_p.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(out_p, bbox_inches="tight")
    plt.close()
    print(f"Successfully generated minimalist techno-economics summary: {out_p}")

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
