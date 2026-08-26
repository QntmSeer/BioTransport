"""Critical Micelle Concentration (CMC) & Critical Salt Concentration (CSC) Analysis."""

from __future__ import annotations
from pathlib import Path
import sys
import numpy as np
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent / "02_continuum_transport" / "python"))
from biotransport.theme import apply_minimalist_theme, PALETTE


def analyze_cmc_csc_phase_boundaries(
    output_png: Path | str = "data/cmc_csc_phase_boundaries.png",
) -> dict[str, float]:
    """Calculates Critical Micelle Concentration (CMC) and Critical Salt Concentration (CSC)
    for amphiphilic and responsive polypeptide coacervation.
    """
    apply_minimalist_theme()

    # Temperature range (15 °C to 45 °C)
    temps_c = np.linspace(15.0, 45.0, 100)
    temps_k = temps_c + 273.15

    # 1. Critical Salt Concentration: CSC(T) = (T_t0 - T) / k_salt (mol/L)
    t_t0_c = 35.0  # Transition temp at 0M salt
    k_salt = 13.5  # K / M (for NaCl)
    csc_nacl_m = np.clip((t_t0_c - temps_c) / (k_salt / 1.0), 0.0, 2.5)

    # 2. Critical Micelle Concentration: ln(CMC) = DeltaG_mic / (R * T)
    # Standard free energy of micellization DeltaG_mic = DeltaH - T*DeltaS
    # Hydrophobic effect: DeltaH > 0 (endothermic unimer desolvation) and DeltaS >> 0 (release of clathrate water)
    delta_h_kj = 25.0  # kJ/mol (endothermic hydrophobic desolvation)
    delta_s_j_k = 150.0  # J/(mol*K) (favorable conformational entropy of freed water)
    r_gas = 8.314e-3  # kJ/(mol*K)

    delta_g_kj = delta_h_kj - temps_k * (delta_s_j_k * 1e-3)
    # CMC in g/L: CMC = exp(DeltaG / RT) * C_standard (with C_standard = 45 g/L)
    cmc_g_l = np.exp(delta_g_kj / (r_gas * temps_k)) * 2500.0

    # 3. 2D Phase Diagram: Salt Concentration [NaCl] vs Polypeptide Concentration C
    # Fixed operating temperature T = 25 °C
    t_op_c = 25.0
    t_op_k = t_op_c + 273.15
    csc_at_top = max(0.0, (t_t0_c - t_op_c) / (k_salt / 1.0))  # ~0.74 M NaCl
    delta_g_top = delta_h_kj - t_op_k * (delta_s_j_k * 1e-3)
    cmc_at_top = np.exp(delta_g_top / (r_gas * t_op_k)) * 2500.0  # ~0.85 g/L

    salt_grid = np.linspace(0.0, 1.5, 100)
    conc_grid = np.linspace(0.05, 15.0, 100)
    S_mesh, C_mesh = np.meshgrid(salt_grid, conc_grid)

    fig, axes = plt.subplots(1, 2, figsize=(14, 5.4), dpi=300)
    fig.suptitle("Thermodynamic Phase Boundaries: CMC & CSC Regimes", fontsize=14, color=PALETTE["slate_dark"])

    # Panel 1: CMC & CSC Temperature Dependences
    ax1 = axes[0]
    ax1.plot(temps_c, csc_nacl_m, color=PALETTE["copper"], lw=2.2, label=r"Critical Salt Concentration $\mathrm{CSC}(T)$")
    ax1.set_xlabel(r"Solution Temperature $T\ (^\circ\mathrm{C})$")
    ax1.set_ylabel(r"Critical Salt Concentration $[\mathrm{NaCl}]_{\mathrm{crit}}\ (\mathrm{M})$", color=PALETTE["copper"])
    ax1.tick_params(axis="y", labelcolor=PALETTE["copper"])

    ax1_twin = ax1.twinx()
    ax1_twin.plot(temps_c, cmc_g_l, color=PALETTE["teal"], linestyle="--", lw=2.2, label=r"Critical Micelle Conc. $\mathrm{CMC}(T)$")
    ax1_twin.set_ylabel(r"$\mathrm{CMC}\ (\mathrm{g} \cdot \mathrm{L}^{-1})$", color=PALETTE["teal"])
    ax1_twin.tick_params(axis="y", labelcolor=PALETTE["teal"])
    ax1_twin.grid(False)

    ax1.set_title(r"A. Temperature Modulated $\mathrm{CSC}$ & $\mathrm{CMC}$", loc="left", color=PALETTE["slate_dark"])

    # Panel 2: 2D State Diagram at 25 °C (Monomer vs Micelle vs Coacervate)
    ax2 = axes[1]

    # Boundary 1: CMC line (horizontal at C = CMC)
    # Boundary 2: CSC line (vertical at [Salt] = CSC)
    # Phase regions:
    # Region I: Unimer Monomer (C < CMC, [Salt] < CSC)
    # Region II: Micelles / Multimers (C > CMC, [Salt] < CSC)
    # Region III: Macroscopic Liquid Coacervate ([Salt] > CSC and C > 0.1 g/L)

    ax2.axhline(y=cmc_at_top, color=PALETTE["teal"], linestyle="--", lw=2.0, label=f"CMC at $25^\\circ\\mathrm{{C}}$ (${cmc_at_top:.2f}\\,\\mathrm{{g/L}}$)")
    ax2.axvline(x=csc_at_top, color=PALETTE["copper"], linestyle="--", lw=2.0, label=f"CSC at $25^\\circ\\mathrm{{C}}$ (${csc_at_top:.2f}\\,\\mathrm{{M}}$)")

    # Shaded regions
    ax2.fill_between([0.0, csc_at_top], 0.0, cmc_at_top, color=PALETTE["slate_light"], alpha=0.15, label="Zone I: Soluble Unimer Monomer")
    ax2.fill_between([0.0, csc_at_top], cmc_at_top, 15.0, color=PALETTE["teal"], alpha=0.12, label="Zone II: Core-Shell Micelles / Oligomers")
    ax2.fill_between([csc_at_top, 1.5], 0.0, 15.0, color=PALETTE["copper"], alpha=0.12, label="Zone III: Macroscopic Liquid Coacervate")

    # Operating setpoint
    ax2.scatter([0.15], [5.0], color="#ffffff", edgecolors=PALETTE["slate_dark"], s=130, marker="o", lw=1.8, zorder=6, label=r"Baseline TFF Feed ($0.15\,\mathrm{M},\ 5\,\mathrm{g/L}$)")

    ax2.set_xlabel(r"Salt Concentration $[\mathrm{NaCl}]\ (\mathrm{mol} \cdot \mathrm{L}^{-1})$")
    ax2.set_ylabel(r"Polypeptide Concentration $C\ (\mathrm{g} \cdot \mathrm{L}^{-1})$")
    ax2.set_title(r"B. State Diagram at $25^\circ\mathrm{C}$ ($[\mathrm{Salt}]$ vs $C_{\mathrm{bulk}}$)", loc="left", color=PALETTE["slate_dark"])
    ax2.set_xlim([0.0, 1.5])
    ax2.set_ylim([0.0, 15.0])
    ax2.legend(loc="upper left", frameon=True, fontsize=8.5)

    plt.tight_layout()
    out_p = Path(output_png)
    out_p.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(out_p, bbox_inches="tight")
    plt.close()
    print(f"Successfully generated CMC & CSC phase boundary summary: {out_p}")

    return {
        "csc_at_25C_M": float(csc_at_top),
        "cmc_at_25C_g_L": float(cmc_at_top),
        "t_t0_C": t_t0_c,
    }


if __name__ == "__main__":
    analyze_cmc_csc_phase_boundaries()
