"""Render animated GIF of continuum filtration time-series dynamics (Minimalist Styling)."""

from pathlib import Path
import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from biotransport.bridge import MdBridgeModel, ProcessSimulator
from biotransport.theme import apply_minimalist_theme, PALETTE


def animate_filtration(
    params_json: Path | str = "data/sample_md_params.json",
    output_gif: Path | str = "data/multiscale_filtration_timeseries.gif",
) -> None:
    apply_minimalist_theme()
    model = MdBridgeModel.load_json(params_json)
    sim = ProcessSimulator(model)
    res = sim.simulate_filtration(tmp_pa=200000.0, total_time_s=3600.0, n_steps=35)

    time_min = np.array(res["time_s"]) / 60.0
    flux = np.array(res["flux_lmh"])
    cw = np.array(res["wall_conc_g_l"])
    rc = np.array(res["cake_resistance_m_inv"]) / 1e12

    y_um = np.linspace(0.0, 20.0, 80)
    decay_length = 4.5

    fig, axes = plt.subplots(1, 2, figsize=(13, 5.2), dpi=140)
    fig.suptitle("Multiscale Bioparticle Transport & Boundary Layer Dynamics", fontsize=13, color=PALETTE["slate_dark"])

    def update(frame_idx):
        ax_flux, ax_profile = axes[0], axes[1]
        ax_flux.clear()
        ax_profile.clear()

        t_sub = time_min[:frame_idx + 1]
        f_sub = flux[:frame_idx + 1]

        ax_flux.plot(time_min, flux, color=PALETTE["slate_light"], alpha=0.4, linestyle="--", label="Full Run Trend")
        ax_flux.plot(t_sub, f_sub, color=PALETTE["slate_dark"], lw=2.2, label=r"Permeate Flux $J(t)$")
        ax_flux.scatter([t_sub[-1]], [f_sub[-1]], color=PALETTE["copper"], s=45, zorder=5)
        ax_flux.fill_between(t_sub, f_sub, color=PALETTE["slate_dark"], alpha=0.10)
        ax_flux.set_xlim(0, 60)
        ax_flux.set_ylim(0, max(flux) * 1.1)
        ax_flux.set_xlabel("Filtration Time (min)")
        ax_flux.set_ylabel(r"Permeate Flux $(\mathrm{L} \cdot \mathrm{m}^{-2} \cdot \mathrm{h}^{-1})$")
        ax_flux.set_title(f"A. Permeate Flux Decline ($t = {t_sub[-1]:.1f}\\,\\mathrm{{min}}$)", loc="left", color=PALETTE["slate_dark"])
        ax_flux.legend(loc="upper right", frameon=True)

        cw_curr = cw[frame_idx]
        c_profile = 5.0 + (cw_curr - 5.0) * np.exp(-y_um / decay_length)
        ax_profile.plot(y_um, c_profile, color=PALETTE["teal"], lw=2.2, label=r"Concentration $C(y, t)$")
        ax_profile.fill_between(y_um, c_profile, color=PALETTE["teal"], alpha=0.12)
        ax_profile.axhline(y=420.0, color=PALETTE["copper"], linestyle=":", lw=1.5, label=r"$C_{\mathrm{gel}} = 420\,\mathrm{g/L}$")
        ax_profile.scatter([0.0], [cw_curr], color=PALETTE["teal"], s=45, zorder=5, label=f"$C_w = {cw_curr:.1f}\\,\\mathrm{{g/L}}$")
        ax_profile.set_xlim(0, 20)
        ax_profile.set_ylim(0, 450)
        ax_profile.set_xlabel(r"Distance from Membrane $y\ (\mu\mathrm{m})$")
        ax_profile.set_ylabel(r"Local Concentration $(\mathrm{g} \cdot \mathrm{L}^{-1})$")
        ax_profile.set_title(f"B. Boundary Layer Profile ($R_c = {rc[frame_idx]:.2f} \\times 10^{{12}}\\,\\mathrm{{m}}^{{-1}}$)", loc="left", color=PALETTE["slate_dark"])
        ax_profile.legend(loc="upper right", frameon=True)
        plt.tight_layout()

    anim = FuncAnimation(fig, update, frames=len(time_min), interval=100)
    out_p = Path(output_gif)
    out_p.parent.mkdir(parents=True, exist_ok=True)
    anim.save(out_p, writer=PillowWriter(fps=8))
    plt.close()
    print(f"Successfully generated: {out_p}")


main = animate_filtration

if __name__ == "__main__":
    main()
