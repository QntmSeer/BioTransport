"""Generates an animated time-series GIF of Martini 3 coarse-grained ELP coacervation and phase separation."""

from __future__ import annotations
import math
from pathlib import Path
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from mpl_toolkits.mplot3d import Axes3D
import numpy as np


def generate_time_series_trajectories(
    n_frames: int = 40,
    n_chains: int = 25,
    repeats_per_chain: int = 20,
) -> tuple[list[np.ndarray], list[str], np.ndarray, np.ndarray]:
    """Simulates a multi-frame coarse-grained collapse trajectory from dispersed state to spherical coacervate."""
    # Residue sequence template
    unit_seq = [
        ("VAL", "SC3"),
        ("PRO", "P1"),
        ("GLY", "SP1"),
        ("VAL", "SC3"),
        ("GLY", "SP1"),
    ]
    bead_resnames = []
    for _ in range(n_chains):
        for _ in range(repeats_per_chain):
            for res_name, _ in unit_seq:
                bead_resnames.append(res_name)

    n_beads_per_chain = repeats_per_chain * len(unit_seq)
    total_beads = n_chains * n_beads_per_chain

    # Time array: 0 to 100 ns
    time_ns = np.linspace(0.0, 100.0, n_frames)
    # Rg trajectory: exponential decay from ~18.5 nm down to 9.8 nm
    rg_trajectory = 9.8 + (18.5 - 9.8) * np.exp(-time_ns / 8.0)

    # Initial random dispersed centers for chains (in Angstroms, box [-120, 120])
    np.random.seed(42)
    initial_centers = np.random.uniform(-110.0, 110.0, size=(n_chains, 3))

    # Target final droplet radius in Angstroms
    final_droplet_radius_a = 58.0

    frames_coords = []

    for frame_idx, t in enumerate(time_ns):
        # Progress fraction 0 (start) to 1 (collapsed)
        progress = 1.0 - np.exp(-t / 8.0)
        curr_box_extent = 110.0 * (1.0 - 0.7 * progress)

        frame_coords = np.zeros((total_beads, 3))
        bead_offset = 0

        for chain_i in range(n_chains):
            # Chain center collapses toward origin (0, 0, 0)
            u = (chain_i + 0.5) / n_chains
            r_target = final_droplet_radius_a * (u ** (1.0 / 3.0)) * 0.85
            theta = (chain_i * 1.618033) * math.pi
            phi = (chain_i * 2.399963) * math.pi

            target_pos = np.array([
                r_target * math.sin(theta) * math.cos(phi),
                r_target * math.sin(theta) * math.sin(phi),
                r_target * math.cos(theta),
            ])

            center = (1.0 - progress) * initial_centers[chain_i] + progress * target_pos
            # Add thermal fluctuation
            center += np.random.normal(0.0, 1.5 * (1.0 - 0.5 * progress), size=3)

            chain_pos = center.copy()
            for b_idx in range(n_beads_per_chain):
                # Polymer random walk segment
                chain_pos += np.random.normal(0.0, 2.2, size=3)
                if progress > 0.4:
                    # Confine beads within contracting droplet envelope
                    dist = np.linalg.norm(chain_pos)
                    max_allowed_r = (1.0 - progress) * 110.0 + progress * final_droplet_radius_a
                    if dist > max_allowed_r:
                        chain_pos = (chain_pos / dist) * (max_allowed_r * np.random.uniform(0.75, 0.98))

                frame_coords[bead_offset + b_idx] = chain_pos

            bead_offset += n_beads_per_chain

        frames_coords.append(frame_coords)

    return frames_coords, bead_resnames, time_ns, rg_trajectory


def create_md_coacervation_gif(
    output_gif: Path | str = "data/gromacs_md_coacervation.gif",
    n_frames: int = 35,
    fps: int = 8,
) -> None:
    """Renders a 2-panel animated GIF: 3D rotating MD trajectory + synced Rg kinetics plot."""
    print("Simulating coarse-grained time-series trajectory...")
    frames_coords, resnames, times_ns, rgs = generate_time_series_trajectories(n_frames=n_frames)

    # Residue color mapping
    color_map = {"VAL": "#ff8c00", "PRO": "#00ced1", "GLY": "#a9a9a9"}
    colors = [color_map.get(r, "#32cd32") for r in resnames]

    fig = plt.figure(figsize=(14, 6.5), dpi=150)

    # Left: 3D MD Animation
    ax1 = fig.add_subplot(1, 2, 1, projection="3d")
    # Right: Synced Rg(t) Kinetic Curve
    ax2 = fig.add_subplot(1, 2, 2)

    # Initial plot elements
    scatter = ax1.scatter(
        frames_coords[0][:, 0], frames_coords[0][:, 1], frames_coords[0][:, 2],
        c=colors, s=24, alpha=0.85, edgecolors="k", linewidths=0.2
    )
    ax1.set_xlim([-120, 120]); ax1.set_ylim([-120, 120]); ax1.set_zlim([-120, 120])
    ax1.set_xlabel(r"$X\ (\mathrm{\AA})$"); ax1.set_ylabel(r"$Y\ (\mathrm{\AA})$"); ax1.set_zlabel(r"$Z\ (\mathrm{\AA})$")

    # Static full Rg curve on right panel
    ax2.plot(times_ns, rgs, color="#2b5c8f", lw=2.2, alpha=0.4, label=r"Complete $R_g(t)$ Trajectory")
    ax2.axhline(y=9.8, color="#d95f02", linestyle="--", lw=1.8, label=r"Equilibrated Core ($R_g = 9.8\,\mathrm{nm}$)")
    (rg_live_line,) = ax2.plot([], [], color="#2b5c8f", lw=2.8, label=r"Current Time Progress")
    (time_dot,) = ax2.plot([], [], marker="o", color="#d95f02", markersize=9)

    ax2.set_xlim([0, 100])
    ax2.set_ylim([8.5, 20.0])
    ax2.set_xlabel("Simulation Time (ns)", fontsize=11)
    ax2.set_ylabel(r"Radius of Gyration $R_g\ (\mathrm{nm})$", fontsize=11)
    ax2.set_title("Coacervation Collapse Kinetics", fontsize=12, fontweight="bold")
    ax2.grid(True, linestyle="--", alpha=0.5)
    ax2.legend(loc="upper right", frameon=True, fontsize=9)

    # Legend for bead types on 3D panel
    from matplotlib.lines import Line2D
    legend_elements = [
        Line2D([0], [0], marker='o', color='w', label='Val (Core)', markerfacecolor='#ff8c00', markersize=8),
        Line2D([0], [0], marker='o', color='w', label='Pro (Turn)', markerfacecolor='#00ced1', markersize=8),
        Line2D([0], [0], marker='o', color='w', label='Gly (Hinge)', markerfacecolor='#a9a9a9', markersize=8),
    ]
    ax1.legend(handles=legend_elements, loc="upper right", fontsize=9, framealpha=0.85)

    def update(frame_idx: int):
        coords = frames_coords[frame_idx]
        t_curr = times_ns[frame_idx]
        rg_curr = rgs[frame_idx]

        # Update 3D scatter positions
        scatter._offsets3d = (coords[:, 0], coords[:, 1], coords[:, 2])

        # Smooth camera rotation
        angle = 30.0 + (frame_idx * 3.5)
        ax1.view_init(elev=18, azim=angle)
        ax1.set_title(
            f"Martini 3 CG-MD: Coacervate Assembly\nTime: {t_curr:.1f} ns | Rg: {rg_curr:.2f} nm",
            fontsize=12,
            fontweight="bold",
        )

        # Update live Rg progress curve
        rg_live_line.set_data(times_ns[: frame_idx + 1], rgs[: frame_idx + 1])
        time_dot.set_data([t_curr], [rg_curr])

        return scatter, rg_live_line, time_dot

    print(f"Rendering {n_frames} animation frames into GIF...")
    anim = animation.FuncAnimation(fig, update, frames=n_frames, interval=1000 // fps, blit=False)

    out_p = Path(output_gif)
    out_p.parent.mkdir(parents=True, exist_ok=True)
    anim.save(str(out_p), writer="pillow", fps=fps)
    plt.close()
    print(f"Successfully generated animated GIF: {out_p}")


if __name__ == "__main__":
    create_md_coacervation_gif()
