"""Orchestrates the GROMACS simulation workflow for Martini 3 ELP assembly."""

from __future__ import annotations
import argparse
import subprocess
import sys
from pathlib import Path


def run_cmd(cmd: list[str], dry_run: bool = False) -> None:
    print(f"[CMD] {' '.join(cmd)}")
    if not dry_run:
        res = subprocess.run(cmd, capture_output=True, text=True)
        if res.returncode != 0:
            print(f"Error executing command: {res.stderr}", file=sys.stderr)
            sys.exit(res.returncode)


def run_gromacs_pipeline(
    work_dir: Path,
    tpr_prefix: str = "elp_quench",
    dry_run: bool = True,
) -> None:
    work_dir = Path(work_dir)
    mdp_dir = work_dir / "mdp"
    top_dir = work_dir / "topologies"

    print("=== Stage 1: Energy Minimization ===")
    run_cmd(["gmx", "grompp", "-f", str(mdp_dir / "em.mdp"), "-c", "initial.gro", "-p", "topol.top", "-o", "em.tpr"], dry_run)
    run_cmd(["gmx", "mdrun", "-v", "-deffnm", "em"], dry_run)

    print("=== Stage 2: NPT Equilibration (300K) ===")
    run_cmd(["gmx", "grompp", "-f", str(mdp_dir / "npt_equilibration.mdp"), "-c", "em.gro", "-p", "topol.top", "-o", "npt.tpr"], dry_run)
    run_cmd(["gmx", "mdrun", "-v", "-deffnm", "npt"], dry_run)

    print("=== Stage 3: Thermal Quench & Coacervation (325K) ===")
    run_cmd(["gmx", "grompp", "-f", str(mdp_dir / "quench_coacervation.mdp"), "-c", "npt.gro", "-p", "topol.top", "-o", f"{tpr_prefix}.tpr"], dry_run)
    run_cmd(["gmx", "mdrun", "-v", "-deffnm", tpr_prefix], dry_run)

    print("=== Stage 4: Analysis & Parameter Extraction ===")
    run_cmd(["gmx", "gyrate", "-f", f"{tpr_prefix}.xtc", "-s", f"{tpr_prefix}.tpr", "-o", "gyrate.xvg"], dry_run)
    run_cmd(["gmx", "msd", "-f", f"{tpr_prefix}.xtc", "-s", f"{tpr_prefix}.tpr", "-o", "msd.xvg"], dry_run)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="GROMACS Martini 3 pipeline runner")
    parser.add_argument("--dry-run", action="store_true", default=True, help="Print commands without executing")
    args = parser.parse_args()

    run_gromacs_pipeline(Path(__file__).parent.parent, dry_run=args.dry_run)
