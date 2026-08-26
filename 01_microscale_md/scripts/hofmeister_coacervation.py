"""Hofmeister Salt Series & Ionic Strength Dependent Coacervation Engine (Minimalist Styling)."""

from __future__ import annotations
from pathlib import Path
import sys
import numpy as np
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent / "02_continuum_transport" / "python"))
from biotransport.theme import apply_minimalist_theme, PALETTE


def simulate_hofmeister_phase_behavior(
    output_png: Path | str = "data/hofmeister_salt_screening.png",
) -> dict[str, float]:
    apply_minimalist_theme()
    salts = {
        r"$\mathrm{(NH_4)_2SO_4}$ (Kosmotrope)": {"k_salt": 28.5, "color": PALETTE["teal"], "linestyle": "-"},
        r"$\mathrm{NaH_2PO_4}$": {"k_salt": 19.2, "color": PALETTE["blue"], "linestyle": "-"},
        r"$\mathrm{NaCl}$ (Standard)": {"k_salt": 13.5, "color": PALETTE["slate_dark"], "linestyle": "-"},
        r"$\mathrm{NaSCN}$ (Chaotrope)": {"k_salt": -4.2, "color": PALETTE["copper"], "linestyle": "--"},
    }

    salt_conc_M = np.linspace(0.0, 1.5, 100)
    t_t0_k = 308.15

    fig, axes = plt.subplots(1, 2, figsize=(14, 5.2), dpi=300)
    fig.suptitle("Hofmeister Series & Ionic Strength Modulated Coacervation", fontsize=14, color=PALETTE["slate_dark"])

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
            lw=2.2,
        )

    ax1.axhline(y=22.0, color=PALETTE["slate_light"], linestyle=":", lw=1.5, label="Ambient ($22^\\circ\\mathrm{C}$)")
    ax1.set_xlabel(r"Salt Concentration $[\mathrm{Salt}]\ (\mathrm{mol} \cdot \mathrm{L}^{-1})$")
    ax1.set_ylabel(r"Transition Temperature $T_t\ (^\circ\mathrm{C})$")
    ax1.set_title(r"A. Hofmeister $T_t$ Phase Shift ($T_t = T_{t,0} - k_{\mathrm{salt}}[\mathrm{Salt}]$)", loc="left", color=PALETTE["slate_dark"])
    ax1.legend(loc="upper right", frameon=True)

    # Panel 2: Debye Screening Length & Condensate Sizing
    ax2 = axes[1]
    i_strength = np.linspace(0.01, 1.5, 100)
    debye_length_nm = 0.304 / np.sqrt(i_strength)
    r_h_effective_nm = 12.65 * (1.0 + 0.45 * np.sqrt(i_strength))

    ax2.plot(i_strength, debye_length_nm, color=PALETTE["slate_med"], lw=2.2, label=r"Debye Length $\kappa^{-1}\ (\mathrm{nm})$")
    ax2.plot(i_strength, r_h_effective_nm, color=PALETTE["teal"], lw=2.2, label=r"Effective Radius $R_h\ (\mathrm{nm})$")
    ax2.set_xlabel(r"Ionic Strength $I\ (\mathrm{mol} \cdot \mathrm{L}^{-1})$")
    ax2.set_ylabel(r"Characteristic Length Scale $(\mathrm{nm})$")
    ax2.set_title("B. Electrostatic Screening & Condensate Sizing", loc="left", color=PALETTE["slate_dark"])
    ax2.legend(loc="upper right", frameon=True)

    plt.tight_layout()
    out_p = Path(output_png)
    out_p.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(out_p, bbox_inches="tight")
    plt.close()
    print(f"Successfully generated minimalist Hofmeister summary: {out_p}")

    return {
        "baseline_t_t_k": t_t0_k,
        "nacl_k_salt_k_per_m": 13.5,
        "t_t_at_1M_nacl_c": (t_t0_k - 13.5 * 1.0) - 273.15,
    }


if __name__ == "__main__":
    simulate_hofmeister_phase_behavior()
