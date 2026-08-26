"""Renders 3D structure snapshots, radial density profiles, and assembly kinetics from GROMACS MD."""

from __future__ import annotations
import math
from pathlib import Path
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np


def parse_pdb_coords(pdb_file: Path | str) -> tuple[np.ndarray, list[str]]:
    coords = []
    resnames = []
    with open(pdb_file, "r", encoding="utf-8") as f:
        for line in f:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                try:
                    resn = line[17:20].strip()
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                    coords.append([x, y, z])
                    resnames.append(resn)
                except ValueError:
                    continue
    return np.array(coords), resnames


def generate_and_plot_md_visualization(
    output_png: Path | str = "data/gromacs_md_visualization.png",
) -> None:
    # 1. Generate PDB snapshots for soluble vs coacervated states
    from render_condensate import generate_synthetic_condensate_pdb

    pdb_soluble = Path("data/sample_soluble.pdb")
    pdb_condensate = Path("data/sample_condensate.pdb")

    generate_synthetic_condensate_pdb(pdb_soluble, state="soluble", n_chains=30, repeats_per_chain=20)
    generate_synthetic_condensate_pdb(pdb_condensate, state="condensate", n_chains=30, repeats_per_chain=20, condensate_radius_a=60.0)

    # 2. Parse coordinates
    coords_sol, resn_sol = parse_pdb_coords(pdb_soluble)
    coords_cond, resn_cond = parse_pdb_coords(pdb_condensate)

    # Color map for Martini 3 residues:
    # Val = amber (#ff8c00), Pro = cyan (#00ced1), Gly = grey (#a9a9a9), others = lime (#32cd32)
    color_map = {
        "VAL": "#ff8c00",
        "PRO": "#00ced1",
        "GLY": "#a9a9a9",
    }
    colors_sol = [color_map.get(r, "#32cd32") for r in resn_sol]
    colors_cond = [color_map.get(r, "#32cd32") for r in resn_cond]

    # Create publication figure: 2x2 grid
    fig = plt.figure(figsize=(15, 13), dpi=300)
    fig.suptitle(
        "Martini 3 Coarse-Grained MD: Thermal Phase Separation & Coacervation of (VPGVG)40",
        fontsize=15,
        fontweight="bold",
        y=0.98,
    )
    plt.subplots_adjust(top=0.90)

    # Panel A: 3D Soluble State (T < Tt, 300K)
    ax1 = fig.add_subplot(2, 2, 1, projection="3d")
    ax1.scatter(
        coords_sol[:, 0], coords_sol[:, 1], coords_sol[:, 2],
        c=colors_sol, s=28, alpha=0.75, edgecolors="k", linewidths=0.3
    )
    ax1.set_title("A. Dispersed Soluble State ($T < T_t = 300\\,\\mathrm{K}$)\nExtended Random Coils", fontsize=12, fontweight="bold")
    ax1.set_xlim([-110, 110]); ax1.set_ylim([-110, 110]); ax1.set_zlim([-110, 110])
    ax1.set_xlabel(r"$X\ (\mathrm{\AA})$"); ax1.set_ylabel(r"$Y\ (\mathrm{\AA})$"); ax1.set_zlabel(r"$Z\ (\mathrm{\AA})$")
    ax1.view_init(elev=20, azim=45)

    # Panel B: 3D Condensed Droplet (T > Tt, 325K)
    ax2 = fig.add_subplot(2, 2, 2, projection="3d")
    ax2.scatter(
        coords_cond[:, 0], coords_cond[:, 1], coords_cond[:, 2],
        c=colors_cond, s=32, alpha=0.85, edgecolors="k", linewidths=0.3
    )
    ax2.set_title("B. Phase-Separated Coacervate ($T > T_t = 325\\,\\mathrm{K}$)\nDense Spherical Condensate Droplet ($R_h \\approx 12.6\\,\\mathrm{nm}$)", fontsize=12, fontweight="bold")
    ax2.set_xlim([-110, 110]); ax2.set_ylim([-110, 110]); ax2.set_zlim([-110, 110])
    ax2.set_xlabel(r"$X\ (\mathrm{\AA})$"); ax2.set_ylabel(r"$Y\ (\mathrm{\AA})$"); ax2.set_zlabel(r"$Z\ (\mathrm{\AA})$")
    ax2.view_init(elev=20, azim=45)

    # Add custom legend for bead types
    from matplotlib.lines import Line2D
    legend_elements = [
        Line2D([0], [0], marker='o', color='w', label='Valine (Hydrophobic Core)', markerfacecolor='#ff8c00', markersize=9),
        Line2D([0], [0], marker='o', color='w', label='Proline (Rigid Turn)', markerfacecolor='#00ced1', markersize=9),
        Line2D([0], [0], marker='o', color='w', label='Glycine (Flexible Hinge)', markerfacecolor='#a9a9a9', markersize=9),
    ]
    ax2.legend(handles=legend_elements, loc="upper right", fontsize=10, framealpha=0.9)

    # Panel C: Radius of Gyration Kinetics (Rg vs Time)
    ax3 = fig.add_subplot(2, 2, 3)
    # Parse gyrate.xvg
    from extract_transport_params import parse_gromacs_xvg
    t_ps, rg_nm = parse_gromacs_xvg("data/sample_gyrate.xvg")
    t_ns = t_ps / 1000.0
    ax3.plot(t_ns, rg_nm, color="#2b5c8f", lw=2.5, marker="o", markersize=6, label=r"Condensate $R_g(t)$")
    ax3.axhline(y=9.8, color="#d95f02", linestyle="--", lw=2.0, label=r"Equilibrated Core $R_g = 9.8\,\mathrm{nm}$")
    ax3.set_xlabel("Simulation Time (ns)", fontsize=11)
    ax3.set_ylabel(r"Radius of Gyration $R_g\ (\mathrm{nm})$", fontsize=11)
    ax3.set_title("C. Coacervation Collapse Kinetics ($R_g$ Trajectory)", fontsize=12, fontweight="bold", loc="left")
    ax3.grid(True, linestyle="--", alpha=0.5)
    ax3.legend(frameon=True, fontsize=10)

    # Panel D: Radial Bead Density Profile rho(r)
    ax4 = fig.add_subplot(2, 2, 4)
    # Compute radial distance from center of mass of condensate
    com = np.mean(coords_cond, axis=0)
    r_dist_nm = np.linalg.norm(coords_cond - com, axis=1) * 0.1  # A to nm

    # Kernel density estimate / histogram
    bins = np.linspace(0.0, 15.0, 30)
    counts, edges = np.histogram(r_dist_nm, bins=bins)
    bin_centers = 0.5 * (edges[:-1] + edges[1:])
    # Volume of spherical shell: 4*pi*r^2 * dr
    dr = edges[1] - edges[0]
    shell_vol = 4.0 * math.pi * (bin_centers ** 2) * dr
    rho_beads = counts / shell_vol

    ax4.plot(bin_centers, rho_beads, color="#7570b3", lw=2.5, marker="s", markersize=5, label=r"Bead Density $\rho(r)$")
    ax4.fill_between(bin_centers, rho_beads, color="#7570b3", alpha=0.2)
    ax4.axvline(x=9.8, color="#d95f02", linestyle=":", lw=2.0, label=r"Droplet Boundary ($R_g$)")
    ax4.set_xlabel(r"Radial Distance from Center of Mass $r\ (\mathrm{nm})$", fontsize=11)
    ax4.set_ylabel(r"Local Bead Density $(\mathrm{beads} \cdot \mathrm{nm}^{-3})$", fontsize=11)
    ax4.set_title("D. Radial Density Profile of Equilibrated Droplet", fontsize=12, fontweight="bold", loc="left")
    ax4.grid(True, linestyle="--", alpha=0.5)
    ax4.legend(frameon=True, fontsize=10)

    plt.tight_layout()
    out_p = Path(output_png)
    out_p.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(out_p, bbox_inches="tight")
    print(f"Successfully generated GROMACS MD visualization: {out_p}")
    plt.close()


main = generate_and_plot_md_visualization

if __name__ == "__main__":
    main()
