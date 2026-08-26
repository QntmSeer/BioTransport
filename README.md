<p align="center">
  <img src="docs/assets/logo.png" width="180" alt="BioTransport Logo" />
</p>

<h1 align="center">BioTransport</h1>
<h3 align="center">High-Performance Multiscale Bioparticle Transport & Phase-Separation Engine</h3>

<p align="center">
  <a href="LICENSE"><img src="https://img.shields.io/badge/License-MIT-blue.svg" alt="License: MIT"></a>
  <a href="01_microscale_md/"><img src="https://img.shields.io/badge/MD-GROMACS%20Martini%203-0f766e" alt="MD"></a>
  <a href="02_continuum_transport/biotransport-rs/"><img src="https://img.shields.io/badge/Engine-Rust%201.75+-1d4ed8" alt="Rust Engine"></a>
  <a href="02_continuum_transport/python/"><img src="https://img.shields.io/badge/Python-3.10%2B-1e293b" alt="Python"></a>
  <a href="scripts/master_run_pipeline.py"><img src="https://img.shields.io/badge/Modules-15%20Verified-047857" alt="Modules"></a>
</p>

---

## 1. Overview & Multiscale Paradigm

**BioTransport** is a first-principles multiscale computational physics and bioprocess engineering suite. It directly couples **Microscale Coarse-Grained Molecular Dynamics (Martini 3 / GROMACS)** of phase-separating biopolypeptides with **Macroscale Tangential Flow Filtration (TFF)** and industrial purification operations.

```
       MICROSCALE MD                    CONTINUUM TRANSPORT                 INDUSTRIAL SCALE-UP
 ┌──────────────────────┐             ┌──────────────────────┐             ┌─────────────────────┐
 │ • Martini 3 (VPGVG)₄₀│             │ • biotransport-rs    │             │ • Fed-Batch (10x)   │
 │ • Coacervate Core Rg │  ────────>  │ • Darcy-Starling     │  ────────>  │ • Diafiltration (8DV│
 │ • Hofmeister Series  │             │ • Virial Osmotic Π   │             │ • CIP Reusability   │
 │ • CMC & CSC Regimes  │             │ • Compressible Cake  │             │ • SEC Pareto Basin  │
 └──────────────────────┘             └──────────────────────┘             └─────────────────────┘
```

1. **Microscale Coacervation MD (`01_microscale_md`):** Simulates sequence-dependent peptide collapse, droplet assembly ($R_g \approx 9.8\,\text{nm}$, $R_h \approx 12.65\,\text{nm}$), Flory-Huggins binodal phase coexistence, Hofmeister salt screening ($T_t([\mathrm{NaCl}])$), CMC/CSC thermodynamic regimes, and capillary droplet ripening ($v_{\text{cap}} = \gamma/\eta$).
2. **Standardized Closure Contract (`data/sample_md_params.json`):** Bridges MD closures ($R_h$, $\rho_p$, $D_0$, $B_2$, $C_{\text{gel}}$, $\mu(C)$) to continuum solvers.
3. **Macroscale Continuum Engine (`02_continuum_transport/biotransport-rs`):** Pure Rust 1D PDE boundary layer solver coupled with Darcy-Starling compressible cake kinetics ($r_c \propto \Delta P^{0.45}$) and Lévêque back-diffusion.
4. **Bioprocess Purification & Scale-Up:** Sizing Fed-Batch up-concentration ($10\times\mathrm{VCF}$), Constant-Volume Diafiltration ($7\times\mathrm{DV}$, $>99.9\%$ impurity clearance), multi-cycle Cleaning-in-Place (CIP, $97.3\%$ flux recovery), and 2D Specific Energy Consumption ($\text{SEC}$) Pareto optimization.

---

## 2. The 15 Unified Simulation Modules

| Domain | # | Module | Generated Asset |
| :--- | :--- | :--- | :--- |
| **Microscale Thermodynamics** | **1** | Martini 3 CG-MD Droplet Assembly & $\rho(r)$ | `data/gromacs_md_visualization.png` |
| | **2** | 3D Rotating MD Trajectory GIF (20 ns) | `data/gromacs_md_coacervation.gif` |
| | **3** | Hofmeister Series Salt Screening ($T_t([\mathrm{Salt}])$) | `data/hofmeister_salt_screening.png` |
| | **4** | Critical Micelle (CMC) & Critical Salt (CSC) Phase Boundaries | `data/cmc_csc_phase_boundaries.png` |
| | **5** | Droplet Coalescence & Ripening Kinetics ($v_{\text{cap}} = \gamma/\eta$) | `data/droplet_coalescence_kinetics.png` |
| | **6** | Flory-Huggins Binodal & Spinodal Phase Diagram | `data/coacervation_phase_diagram.png` |
| **Continuum Transport** | **7** | TFF Boundary Layer & Fouling Dynamics ($J(t), C_w, R_c$) | `data/filtration_summary_figure.png` |
| | **8** | Continuum Filtration Dynamics Animated GIF | `data/multiscale_filtration_timeseries.gif` |
| | **9** | Parametric Limiting Flux Profiles across Shear Rates $\dot{\gamma}$ | `data/limiting_flux_curves.png` |
| | **10** | 2D Optimization Map with Iso-Flux & SEC Contours | `data/intense_optimization_landscape.png` |
| **Bioprocess Operations** | **11** | Hermia's 4-Mechanism Fouling Diagnostic Breakdown | `data/hermia_fouling_analysis.png` |
| | **12** | Fed-Batch Concentration ($10\times\mathrm{VCF}$) + Diafiltration | `data/fed_batch_concentration_diafiltration.png` |
| | **13** | Constant-Volume Diafiltration (CVD) Buffer Exchange | `data/diafiltration_summary.png` |
| | **14** | Cleaning-in-Place (CIP) & Fouling Reversibility Breakdown | `data/cip_fouling_reversibility.png` |
| **Industrial Scale-Up** | **15** | Techno-Economics (100L Batch) & Engineering Summary | `data/techno_economics_summary.png`<br>`data/MULTISCALE_SIMULATION_REPORT.md` |

---

## 3. Directory Layout

```text
BioTransport/
├── README.md
├── LICENSE
├── pyproject.toml
├── docs/assets/logo.png
├── data/
│   ├── sample_md_params.json
│   ├── MULTISCALE_SIMULATION_REPORT.md
│   └── *.png / *.gif
├── 01_microscale_md/
│   ├── topologies/
│   ├── mdp/
│   └── scripts/
│       ├── plot_md_snapshots.py
│       ├── animate_md_trajectory.py
│       ├── hofmeister_coacervation.py
│       ├── cmc_csc_analysis.py
│       ├── droplet_coalescence.py
│       ├── phase_diagram.py
│       └── extract_transport_params.py
├── 02_continuum_transport/
│   ├── biotransport-rs/               # High-Performance Rust Crate
│   │   ├── Cargo.toml
│   │   └── src/
│   └── python/biotransport/           # Python API & Physics Modules
│       ├── bridge.py
│       ├── theme.py                   # Minimalist Academic Palette
│       ├── visualize.py
│       ├── parameter_sweep.py
│       ├── hermia_fouling.py
│       ├── fed_batch_purification.py
│       ├── diafiltration.py
│       ├── cip_reversibility.py
│       └── techno_economics.py
└── scripts/
    ├── master_run_pipeline.py         # Unified 15-Module Master Runner
    ├── generate_multiscale_report.py
    └── sync_and_run_on_agni.ps1       # HPC Synchronization & Execution
```

---

## 4. Quickstart

### Execute Entire 15-Module Suite

```bash
# Python Environment
pip install -e .

# Run all 15 stages end-to-end:
python scripts/master_run_pipeline.py
```

### HPC Workstation Deployment (Agni)

```powershell
# Sync and execute remotely on Agni HPC:
.\scripts\sync_and_run_on_agni.ps1
```

---

## 5. References

* **Macromolecules (2023).** *"Coacervate Formation of Elastin-like Polypeptides in Explicit Aqueous Solution Using Coarse-Grained Molecular Dynamics Simulations."* *Macromolecules*, 56(18), 7543–7555.
* **Bacchin, P., et al. (2006).** *"Critical and sustainable fluxes in membrane filtration: Definitions, literature review, and new analysis."* *J. Membr. Sci.*, 281(1-2), 42–69.
* **Kjølbye, L. R., et al. (2024).** *"Martini 3 building blocks for lipid nanoparticle design."* *J. Chem. Theory Comput.*
