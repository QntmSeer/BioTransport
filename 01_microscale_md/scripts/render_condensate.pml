# PyMOL Rendering Script for Martini 3 Coarse-Grained ELP Condensates
# Usage: pymol render_condensate.pml (or File -> Run Script inside PyMOL)

# 1. Environment & Ray-tracing quality
reinitialize
set ray_shadows, 1
set ray_trace_mode, 1
set ray_trace_fog, 0
set antialias, 2
set ambient, 0.25
set direct, 0.70
set reflect, 0.40
bg_color white

# 2. Load system (replace with your .gro or .pdb structure)
# load ../sample_condensate.pdb, elp_system

# 3. Representation: Spheres for CG beads
hide everything
show spheres

# 4. Bead size scaling for Martini 3 (standard van der Waals diameter ~0.47 nm = 2.35 A radius)
alter all, vdw=2.35
rebuild

# 5. Color palette by amino acid hydrophobicity:
# - Valine (hydrophobic core): Amber / Warm Orange
color 0xffa500, resn VAL
# - Proline (rigid turn): Cyan / Teal
color 0x00ced1, resn PRO
# - Glycine (flexible hinge): Light Grey
color 0xd3d3d3, resn GLY
# - Alanine / Lysine / Glutamate (if present)
color 0x32cd32, resn ALA
color 0x4169e1, resn LYS
color 0xff4500, resn GLU

# - Backbone beads (BB): Slight transparency or highlight
set sphere_transparency, 0.15, name BB

# 6. Center and orient camera
orient
zoom all, 5

# 7. Render high-resolution publication image
# ray 2400, 2400
# png elp_martini3_condensate_render.png, dpi=300
