"""Liquid-Liquid Condensate Droplet Coalescence and Ripening Kinetics (Minimalist Styling)."""

from __future__ import annotations
from pathlib import Path
import sys
import numpy as np
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent / "02_continuum_transport" / "python"))
from biotransport.theme import apply_minimalist_theme, PALETTE


def simulate_droplet_coalescence(
    gamma_interfacial_un_m: float = 45.0,
    eta_condensate_pa_s: float = 1.2,
    output_png: Path | str = "data/droplet_coalescence_kinetics.png",
) -> dict[str, float]:
    apply_minimalist_theme()
    times_s = np.linspace(0.0, 1800.0, 100)

    gamma_n_m = gamma_interfacial_un_m * 1e-6
    v_cap_m_s = gamma_n_m / eta_condensate_pa_s

    r0_nm = 12.65
    k_lsw_nm3_s = 450.0
    r_ostwald_nm = (r0_nm**3 + k_lsw_nm3_s * times_s) ** (1.0 / 3.0)

    shear_rate = 4000.0
    alpha_coalescence = 2.5e-7
    r_shear_nm = r0_nm * (1.0 + alpha_coalescence * shear_rate * times_s) ** 0.85

    r_grid = np.linspace(5.0, 300.0, 200)
    time_snapshots = [0, 300, 900, 1800]
    minimal_colors = [PALETTE["slate_light"], PALETTE["teal_light"], PALETTE["teal"], PALETTE["slate_dark"]]

    fig, axes = plt.subplots(1, 2, figsize=(14, 5.2), dpi=300)
    fig.suptitle("Bioparticle Condensate Coalescence & Droplet Growth Kinetics", fontsize=14, color=PALETTE["slate_dark"])

    # Panel 1: Droplet Radius Evolution
    ax1 = axes[0]
    t_min = times_s / 60.0
    ax1.plot(t_min, r_ostwald_nm, color=PALETTE["slate_med"], linestyle="--", lw=2.0, label=r"Ostwald Ripening ($R \propto t^{1/3}$)")
    ax1.plot(t_min, r_shear_nm, color=PALETTE["teal"], lw=2.4, label=r"Shear Coalescence ($\dot{\gamma} = 4,000\,\mathrm{s}^{-1}$)")
    ax1.set_xlabel("Incubation / Shearing Time (min)")
    ax1.set_ylabel(r"Mean Droplet Radius $R\ (\mathrm{nm})$")
    ax1.set_title(r"A. Droplet Growth Trajectory ($v_{\mathrm{cap}} = \gamma / \eta$)", loc="left", color=PALETTE["slate_dark"])
    ax1.legend(loc="upper left", frameon=True)

    # Panel 2: Evolving Droplet Size Distribution
    ax2 = axes[1]
    for t_snap, col in zip(time_snapshots, minimal_colors):
        mean_r = float(r0_nm * (1.0 + alpha_coalescence * shear_rate * t_snap) ** 0.85)
        sigma = mean_r * 0.28
        p_dist = (1.0 / (sigma * np.sqrt(2 * np.pi))) * np.exp(-0.5 * ((r_grid - mean_r) / sigma) ** 2)
        ax2.plot(r_grid, p_dist, color=col, lw=2.0, label=f"$t = {t_snap // 60}\\,\\mathrm{{min}}$ ($\\bar{{R}} = {mean_r:.1f}\\,\\mathrm{{nm}}$)")
        ax2.fill_between(r_grid, p_dist, color=col, alpha=0.12)

    ax2.set_xlabel(r"Droplet Radius $R\ (\mathrm{nm})$")
    ax2.set_ylabel(r"Probability Density $P(R)$")
    ax2.set_title("B. Droplet Size Distribution Evolution", loc="left", color=PALETTE["slate_dark"])
    ax2.legend(loc="upper right", frameon=True)

    plt.tight_layout()
    out_p = Path(output_png)
    out_p.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(out_p, bbox_inches="tight")
    plt.close()
    print(f"Successfully generated minimalist droplet coalescence summary: {out_p}")

    return {
        "capillary_velocity_um_s": float(v_cap_m_s * 1e6),
        "initial_r0_nm": r0_nm,
        "final_r_shear_nm": float(r_shear_nm[-1]),
    }


if __name__ == "__main__":
    simulate_droplet_coalescence()
