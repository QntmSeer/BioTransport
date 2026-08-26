"""Automated Martini 3 Coarse-Grained Topology Generator for Elastin-Like Polypeptides (ELPs).

Constructs GROMACS .itp files for repetitive (VPGXG)_n sequences using standard
Martini 3 bead mappings (BB, SC1, SC2) and elastic network models (ENM) or bonded parameters.
"""

from __future__ import annotations
import argparse
from pathlib import Path

# Standard Martini 3 amino acid mapping dictionary
MARTINI3_BEAD_MAPPING: dict[str, list[tuple[str, str, float, float]]] = {
    # ResName: [(BeadName, BeadType, Mass_amu, Charge_e), ...]
    "VAL": [("BB", "P2", 72.0, 0.0), ("SC1", "SC3", 43.0, 0.0)],
    "PRO": [("BB", "P1", 97.0, 0.0)],
    "GLY": [("BB", "SP1", 57.0, 0.0)],
    "ALA": [("BB", "P2", 71.0, 0.0)],
    "LYS": [("BB", "P2", 72.0, 0.0), ("SC1", "C3", 42.0, 0.0), ("SC2", "TQ1p", 30.0, 1.0)],
    "GLU": [("BB", "P2", 72.0, 0.0), ("SC1", "Qa", 58.0, -1.0)],
}


def generate_elp_itp(
    repeat_unit: str = "VPGVG",
    repeats: int = 40,
    output_path: Path | str = "elp_vpgvg.itp",
    molecule_name: str = "ELP_VPGVG40",
) -> str:
    """Generates a complete Martini 3 .itp topology string for an ELP chain."""
    full_sequence = []
    # Expand 1-letter to 3-letter codes
    code_map = {"V": "VAL", "P": "PRO", "G": "GLY", "A": "ALA", "K": "LYS", "E": "GLU"}
    for _ in range(repeats):
        for char in repeat_unit.upper():
            if char in code_map:
                full_sequence.append(code_map[char])

    lines = [
        f"; Martini 3 Topology for {molecule_name}",
        f"; Sequence: ({repeat_unit})_{repeats} ({len(full_sequence)} residues)",
        "; Generated automatically for Multiscale Bioparticle Transport Framework",
        "",
        "[ moleculetype ]",
        f"; Name nrexcl",
        f"{molecule_name} 1",
        "",
        "[ atoms ]",
        "; nr type resnr residue atom cgnr charge mass",
    ]

    atom_idx = 1
    res_idx = 1
    bb_indices = []

    for resname in full_sequence:
        beads = MARTINI3_BEAD_MAPPING.get(resname, [("BB", "P2", 72.0, 0.0)])
        for b_name, b_type, b_mass, b_charge in beads:
            lines.append(
                f"{atom_idx:>5} {b_type:>5} {res_idx:>5} {resname:>5} {b_name:>5} {atom_idx:>5} {b_charge:>7.3f} {b_mass:>7.3f}"
            )
            if b_name == "BB":
                bb_indices.append(atom_idx)
            atom_idx += 1
        res_idx += 1

    lines.extend([
        "",
        "[ bonds ]",
        "; Backbone sequential bonds (standard Martini 3 BB-BB equilibrium ~0.35 nm, k=7000 kJ/mol/nm^2)",
        "; ai aj funct length_nm force_const",
    ])

    for i in range(len(bb_indices) - 1):
        lines.append(f"{bb_indices[i]:>5} {bb_indices[i+1]:>5} 1 0.350 7000.0")

    lines.extend([
        "",
        "[ angles ]",
        "; Backbone angles (127.0 deg, k=30.0 kJ/mol/rad^2)",
        "; ai aj ak funct angle_deg force_const",
    ])

    for i in range(len(bb_indices) - 2):
        lines.append(f"{bb_indices[i]:>5} {bb_indices[i+1]:>5} {bb_indices[i+2]:>5} 2 127.0 30.0")

    content = "\n".join(lines) + "\n"

    out_p = Path(output_path)
    out_p.parent.mkdir(parents=True, exist_ok=True)
    out_p.write_text(content, encoding="utf-8")
    return content


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate Martini 3 ELP topology")
    parser.add_argument("--repeat", type=str, default="VPGVG", help="ELP repeat sequence (e.g. VPGVG)")
    parser.add_argument("--n", type=int, default=40, help="Number of repeats (default: 40)")
    parser.add_argument("--out", type=str, default="elp_vpgvg.itp", help="Output .itp file path")
    args = parser.parse_args()

    generate_elp_itp(repeat_unit=args.repeat, repeats=args.n, output_path=args.out)
    print(f"Successfully generated ELP topology: {args.out}")
