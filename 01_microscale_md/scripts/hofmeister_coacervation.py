"""Hofmeister Salt Series & Ionic Strength Dependent Coacervation Engine."""

from __future__ import annotations
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt


def simulate_hofmeister_phase_behavior(
    output_png: Path | str = "data/hofmeister_salt_screening.png",
) -> dict[str, float]:
    """Models salt-induced transition temperature depression and electrostatic screening."""
    salts = {
        r"$\mathrm{(NH_4)_2SO_4}$ (Kosmotrope)": {"k_salt": 28.5, "color": "#1f77b4", "linestyle": "-"},
        r"$\mathrm{NaH_2PO_4}$": {"k_salt": 19.2, "color": "#ff7f0e", "linestyle": "-"},
        r"$\mathrm{NaCl}$ (Standard)": {"k_salt": 13.5, "color": "#2ca02c", "linestyle": "-"},
        r"$\mathrm{NaSCN}$ (Chaotrope)": {"k_salt": -4.2, "color": "#d62728", "linestyle": "--"},
    }

    salt_conc_M = np.linspace(0.0, 1.5, 100)
    t_t0_k = 308.15  # 35 °C baseline transition temp

    fig, axes = plt.subplots(1, 2, figsize=(14, 5.5), dpi=300)
    fig.suptitle("Hofmeister Series & Ionic Strength Modulated Coacervation", fontsize=14, fontweight="bold")

    # Panel 1: Transition Temperature Shift T_t([Salt])
    ax1 = axes[0]
    for name, props in salts.items():
        t_t = t_t0_k - props["k_salt"] * salt_conc_M
        ax1.plot(
            salt_conc_M,
            t_t - 273.15,
            label=f"{name} ($k = {props['k_salt']:+.1f}\\,\\mathrm{{K/M}}$)",
            color=props["color"],
            linestyle=props["linestyle"],
            lw=2.5,
        )

    # Operating line for ambient / room temperature
    ax1.axhline(y=22.0, color="gray", linestyle=":", lw=1.8, label="Ambient Temperature ($22^\\circ\\mathrm{C}$)")
    ax1.set_xlabel(r"Salt Concentration $[\mathrm{Salt}]\ (\mathrm{mol} \cdot \mathrm{L}^{-1})$", fontsize=11)
    ax1.set_ylabel(r"Transition Temperature $T_t\ (^\circ\mathrm{C})$", fontsize=11)
    ax1.set_title(r"A. Hofmeister $T_t$ Phase Shift ($T_t = T_{t,0} - k_{\mathrm{salt}}[\mathrm{Salt}]$)", fontsize=12, fontweight="bold", loc="left")
    ax1.grid(True, linestyle="--", alpha=0.5)
    ax1.legend(loc="upper right", frameon=True, fontsize=9)

    # Panel 2: Debye Screening Length & Condensate Particle Sizing
    ax2 = axes[1]
    # Debye length lambda_D = 0.304 / sqrt(I) nm for 1:1 monovalent electrolyte in water at 298K
    i_strength = np.linspace(0.01, 1.5, 100)
    debye_length_nm = 0.304 / np.sqrt(i_strength)

    # Hydrophobic condensation efficiency: particle size increases with Debye screening & dehydration
    r_h_effective_nm = 12.65 * (1.0 + 0.45 * np.sqrt(i_strength))

    ax2.plot(i_strength, debye_length_nm, "m-", lw=2.5, label=r"Debye Length $\kappa^{-1}\ (\mathrm{nm})$")
    ax2.plot(i_strength, r_h_effective_nm, "b-", lw=2.5, label=r"Effective Condensate Radius $R_h\ (\mathrm{nm})$")
    ax2.set_xlabel(r"Ionic Strength $I\ (\mathrm{mol} \cdot \mathrm{L}^{-1})$", fontsize=11)
    ax2.set_ylabel(r"Characteristic Length Scale $(\mathrm{nm})$", fontsize=11)
    ax2.set_title("B. Electrostatic Screening & Condensate Swelling", fontsize=12, fontweight="bold", loc="left")
    ax2.grid(True, linestyle="--", alpha=0.5)
    ax2.legend(loc="upper right", frameon=True, fontsize=9)

    plt.tight_layout()
    out_p = Path(output_png)
    out_p.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(out_p, bbox_inches="tight")
    plt.close()
    print(f"Successfully generated Hofmeister screening summary: {out_p}")

    return {
        "baseline_t_t_k": t_t0_k,
        "nacl_k_salt_k_per_m": 13.5,
        "t_t_at_1M_nacl_c": (t_t0_k - 13.5 * 1.0) - 273.15,
    }


if __name__ == "__main__":
    simulate_hofmeister_phase_behavior()
