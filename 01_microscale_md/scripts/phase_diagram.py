"""Flory-Huggins / Voorn-Overbeek Phase Diagram (Minimalist Academic Styling)."""

from __future__ import annotations
from pathlib import Path
import sys
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as opt

sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent / "02_continuum_transport" / "python"))
from biotransport.theme import apply_minimalist_theme, PALETTE


def flory_huggins_phase_diagram(
    n_monomers: int = 200,
    t_transition_k: float = 308.15,
    output_png: Path | str = "data/coacervation_phase_diagram.png",
) -> dict[str, float]:
    apply_minimalist_theme()
    phi_c = 1.0 / (1.0 + np.sqrt(n_monomers))
    chi_c = 0.5 * (1.0 + 1.0 / np.sqrt(n_monomers)) ** 2
    chi_slope = 0.045
    temps_k = np.linspace(290.0, 340.0, 100)

    phi_spinodal_dilute = []
    phi_spinodal_dense = []
    phi_binodal_dilute = []
    phi_binodal_dense = []
    valid_temps = []

    for t in temps_k:
        chi = chi_c + chi_slope * (t - t_transition_k)
        if chi <= chi_c:
            continue

        valid_temps.append(t)
        a = 2.0 * chi * n_monomers
        b = -(2.0 * chi * n_monomers + 1.0 - n_monomers)
        c = 1.0
        discriminant = b**2 - 4.0 * a * c
        if discriminant >= 0:
            s1 = (-b - np.sqrt(discriminant)) / (2.0 * a)
            s2 = (-b + np.sqrt(discriminant)) / (2.0 * a)
            phi_spinodal_dilute.append(min(s1, s2))
            phi_spinodal_dense.append(max(s1, s2))
        else:
            phi_spinodal_dilute.append(phi_c)
            phi_spinodal_dense.append(phi_c)

        def common_tangent(p):
            p1, p2 = p
            if p1 <= 1e-6 or p1 >= 0.99 or p2 <= 1e-6 or p2 >= 0.99 or p1 >= p2:
                return [1e3, 1e3]
            mu_diff = (
                (1.0 / n_monomers) * np.log(p1)
                - np.log(1.0 - p1)
                + chi * (1.0 - 2.0 * p1)
                - ((1.0 / n_monomers) * np.log(p2) - np.log(1.0 - p2) + chi * (1.0 - 2.0 * p2))
            )
            pi_diff = (
                -(1.0 / n_monomers) * p1
                - np.log(1.0 - p1)
                - p1
                - chi * (p1**2)
                - (-(1.0 / n_monomers) * p2 - np.log(1.0 - p2) - p2 - chi * (p2**2))
            )
            return [mu_diff, pi_diff]

        try:
            sol = opt.fsolve(common_tangent, [phi_spinodal_dilute[-1] * 0.4, min(0.85, phi_spinodal_dense[-1] * 1.2)])
            if 0 < sol[0] < sol[1] < 1.0:
                phi_binodal_dilute.append(sol[0])
                phi_binodal_dense.append(sol[1])
            else:
                phi_binodal_dilute.append(phi_spinodal_dilute[-1])
                phi_binodal_dense.append(phi_spinodal_dense[-1])
        except Exception:
            phi_binodal_dilute.append(phi_spinodal_dilute[-1])
            phi_binodal_dense.append(phi_spinodal_dense[-1])

    valid_temps = np.array(valid_temps)
    rho_polymer_g_l = 1350.0

    fig, ax = plt.subplots(figsize=(8.5, 5.8), dpi=300)

    ax.plot(
        np.array(phi_binodal_dilute) * rho_polymer_g_l,
        valid_temps,
        color=PALETTE["teal"],
        lw=2.4,
        label=r"Binodal (Coexistence Curve)",
    )
    ax.plot(
        np.array(phi_binodal_dense) * rho_polymer_g_l,
        valid_temps,
        color=PALETTE["teal"],
        lw=2.4,
    )
    ax.plot(
        np.array(phi_spinodal_dilute) * rho_polymer_g_l,
        valid_temps,
        color=PALETTE["copper"],
        linestyle="--",
        lw=1.8,
        label=r"Spinodal Boundary",
    )
    ax.plot(
        np.array(phi_spinodal_dense) * rho_polymer_g_l,
        valid_temps,
        color=PALETTE["copper"],
        linestyle="--",
        lw=1.8,
    )

    ax.scatter(
        [phi_c * rho_polymer_g_l],
        [t_transition_k],
        color="#ffffff",
        edgecolors=PALETTE["slate_dark"],
        s=110,
        lw=1.8,
        zorder=6,
        label=f"Critical Point ($T_c = {t_transition_k:.1f}\\,\\mathrm{{K}}$, $C_c = {phi_c * rho_polymer_g_l:.1f}\\,\\mathrm{{g/L}}$)",
    )

    ax.fill_betweenx(
        valid_temps,
        np.array(phi_binodal_dilute) * rho_polymer_g_l,
        np.array(phi_binodal_dense) * rho_polymer_g_l,
        color=PALETTE["teal"],
        alpha=0.08,
        label="2-Phase Coacervate Regime",
    )

    ax.set_xlabel(r"Polypeptide Concentration $C\ (\mathrm{g} \cdot \mathrm{L}^{-1})$")
    ax.set_ylabel(r"Temperature $T\ (\mathrm{K})$")
    ax.set_title("ELP Coacervation Phase Diagram (Flory-Huggins Thermodynamics)", loc="left", color=PALETTE["slate_dark"])
    ax.legend(loc="upper left", frameon=True)

    out_p = Path(output_png)
    out_p.parent.mkdir(parents=True, exist_ok=True)
    plt.tight_layout()
    plt.savefig(out_p, bbox_inches="tight")
    plt.close()
    print(f"Successfully generated minimalist phase diagram: {out_p}")

    return {
        "critical_temperature_K": t_transition_k,
        "critical_concentration_g_L": float(phi_c * rho_polymer_g_l),
        "critical_chi": float(chi_c),
    }


if __name__ == "__main__":
    flory_huggins_phase_diagram()
