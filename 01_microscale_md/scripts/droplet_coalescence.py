"""Liquid-Liquid Condensate Droplet Coalescence and Ripening Kinetics Module."""

from __future__ import annotations
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt


def simulate_droplet_coalescence(
    gamma_interfacial_un_m: float = 45.0,  # Interfacial tension in uN/m
    eta_condensate_pa_s: float = 1.2,  # Condensate viscosity in Pa.s
    output_png: Path | str = "data/droplet_coalescence_kinetics.png",
) -> dict[str, float]:
    """Models capillary droplet fusion, Ostwald ripening, and shear-induced coalescence."""
    times_s = np.linspace(0.0, 1800.0, 100)  # 0 to 30 min

    # Capillary velocity v_cap = gamma / eta (m/s)
    gamma_n_m = gamma_interfacial_un_m * 1e-6
    v_cap_m_s = gamma_n_m / eta_condensate_pa_s

    # 1. Ostwald Ripening (Lifshitz-Slyozov-Wagner theory): R(t)^3 = R0^3 + K_LSW * t
    r0_nm = 12.65
    k_lsw_nm3_s = 450.0  # nm^3 / s
    r_ostwald_nm = (r0_nm**3 + k_lsw_nm3_s * times_s) ** (1.0 / 3.0)

    # 2. Shear-Accelerated Coalescence (at gamma_dot = 4000 s^-1)
    # Collision frequency scales with shear rate: dR/dt = alpha * gamma_dot * R
    shear_rate = 4000.0
    alpha_coalescence = 2.5e-7
    r_shear_nm = r0_nm * (1.0 + alpha_coalescence * shear_rate * times_s) ** 0.85

    # 3. Dynamic Droplet Size Distribution P(R, t) at t = 0, 5 min, 15 min, 30 min
    r_grid = np.linspace(5.0, 300.0, 200)
    time_snapshots = [0, 300, 900, 1800]  # seconds
    colors = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728"]

    fig, axes = plt.subplots(1, 2, figsize=(14, 5.5), dpi=300)
    fig.suptitle("Bioparticle Condensate Coalescence & Droplet Growth Kinetics", fontsize=14, fontweight="bold")

    # Panel 1: Droplet Radius Evolution R(t)
    ax1 = axes[0]
    t_min = times_s / 60.0
    ax1.plot(t_min, r_ostwald_nm, "b--", lw=2.2, label=r"Ostwald Ripening (Diffusion, $R \propto t^{1/3}$)")
    ax1.plot(t_min, r_shear_nm, "r-", lw=2.5, label=r"Shear-Induced Coalescence ($\dot{\gamma} = 4,000\,\mathrm{s}^{-1}$)")
    ax1.set_xlabel("Incubation / Shearing Time (min)", fontsize=11)
    ax1.set_ylabel(r"Mean Droplet Radius $R\ (\mathrm{nm})$", fontsize=11)
    ax1.set_title(r"A. Droplet Growth Trajectory ($v_{\mathrm{cap}} = \gamma / \eta$)", fontsize=12, fontweight="bold", loc="left")
    ax1.grid(True, linestyle="--", alpha=0.5)
    ax1.legend(loc="upper left", frameon=True)

    # Panel 2: Evolving Droplet Size Distribution P(R)
    ax2 = axes[1]
    for t_snap, col in zip(time_snapshots, colors):
        mean_r = float(r0_nm * (1.0 + alpha_coalescence * shear_rate * t_snap) ** 0.85)
        sigma = mean_r * 0.28
        p_dist = (1.0 / (sigma * np.sqrt(2 * np.pi))) * np.exp(-0.5 * ((r_grid - mean_r) / sigma) ** 2)
        ax2.plot(r_grid, p_dist, color=col, lw=2.2, label=f"$t = {t_snap // 60}\\,\\mathrm{{min}}$ ($\\bar{{R}} = {mean_r:.1f}\\,\\mathrm{{nm}}$)")
        ax2.fill_between(r_grid, p_dist, color=col, alpha=0.15)

    ax2.set_xlabel(r"Droplet Radius $R\ (\mathrm{nm})$", fontsize=11)
    ax2.set_ylabel(r"Probability Density $P(R)$", fontsize=11)
    ax2.set_title("B. Dynamic Droplet Size Distribution Evolution", fontsize=12, fontweight="bold", loc="left")
    ax2.grid(True, linestyle="--", alpha=0.5)
    ax2.legend(loc="upper right", frameon=True, fontsize=9)

    plt.tight_layout()
    out_p = Path(output_png)
    out_p.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(out_p, bbox_inches="tight")
    plt.close()
    print(f"Successfully generated droplet coalescence summary: {out_p}")

    return {
        "capillary_velocity_um_s": float(v_cap_m_s * 1e6),
        "initial_r0_nm": r0_nm,
        "final_r_shear_nm": float(r_shear_nm[-1]),
    }


if __name__ == "__main__":
    simulate_droplet_coalescence()
