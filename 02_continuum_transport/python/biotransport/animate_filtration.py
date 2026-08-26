"""Render animated GIF of continuum filtration time-series dynamics."""

from pathlib import Path
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter
import numpy as np
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from biotransport.bridge import MdBridgeModel, ProcessSimulator


def main():
    model = MdBridgeModel.load_json("data/sample_md_params.json")
    sim = ProcessSimulator(model)
    res = sim.simulate_filtration(tmp_pa=200000.0, total_time_s=3600.0, n_steps=35)

    time_min = np.array(res["time_s"]) / 60.0
    flux = np.array(res["flux_lmh"])
    cw = np.array(res["wall_conc_g_l"])
    rc = np.array(res["cake_resistance_m_inv"]) / 1e12

    y_um = np.linspace(0.0, 20.0, 80)
    decay_length = 4.5

    fig, axes = plt.subplots(1, 2, figsize=(13, 5.2), dpi=140)
    fig.suptitle("Multiscale Bioparticle Transport and Fouling Dynamics (Time-Series)", fontsize=13, fontweight="bold")

    def update(frame_idx):
        ax_flux, ax_profile = axes[0], axes[1]
        ax_flux.clear()
        ax_profile.clear()

        t_sub = time_min[:frame_idx + 1]
        f_sub = flux[:frame_idx + 1]

        ax_flux.plot(time_min, flux, color="gray", alpha=0.3, linestyle="--", label="60-min Trend")
        ax_flux.plot(t_sub, f_sub, color="#1f77b4", lw=2.5, label="Permeate Flux J(t)")
        ax_flux.scatter([t_sub[-1]], [f_sub[-1]], color="#d62728", s=50, zorder=5)
        ax_flux.fill_between(t_sub, f_sub, color="#1f77b4", alpha=0.15)
        ax_flux.set_xlim(0, 60)
        ax_flux.set_ylim(0, max(flux) * 1.1)
        ax_flux.set_xlabel("Filtration Time (min)", fontsize=11)
        ax_flux.set_ylabel("Permeate Flux (L / m^2 / h)", fontsize=11)
        ax_flux.set_title(f"A. Permeate Flux Decline (t = {t_sub[-1]:.1f} min)", fontsize=11, fontweight="bold", loc="left")
        ax_flux.grid(True, linestyle="--", alpha=0.5)
        ax_flux.legend(loc="upper right")

        cw_curr = cw[frame_idx]
        c_profile = 5.0 + (cw_curr - 5.0) * np.exp(-y_um / decay_length)
        ax_profile.plot(y_um, c_profile, color="#ff7f0e", lw=2.5, label="Concentration C(y, t)")
        ax_profile.fill_between(y_um, c_profile, color="#ff7f0e", alpha=0.2)
        ax_profile.axhline(y=420.0, color="darkred", linestyle=":", lw=1.8, label="Gel Limit C_gel")
        ax_profile.scatter([0.0], [cw_curr], color="#d62728", s=50, zorder=5, label=f"Cw = {cw_curr:.1f} g/L")
        ax_profile.set_xlim(0, 20)
        ax_profile.set_ylim(0, 450)
        ax_profile.set_xlabel("Distance from Membrane y (um)", fontsize=11)
        ax_profile.set_ylabel("Local Concentration (g/L)", fontsize=11)
        ax_profile.set_title(f"B. Boundary Layer Profile (Cake Rc = {rc[frame_idx]:.2f} x 10^12 m^-1)", fontsize=11, fontweight="bold", loc="left")
        ax_profile.grid(True, linestyle="--", alpha=0.5)
        ax_profile.legend(loc="upper right")
        plt.tight_layout()

    anim = FuncAnimation(fig, update, frames=len(time_min), interval=100)
    out_p = Path("data/multiscale_filtration_timeseries.gif")
    out_p.parent.mkdir(parents=True, exist_ok=True)
    anim.save(out_p, writer=PillowWriter(fps=8))
    plt.close()
    print(f"Successfully generated: {out_p}")


if __name__ == "__main__":
    main()
