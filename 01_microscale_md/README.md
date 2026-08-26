# 01_microscale_md: Coarse-Grained MD Topologies & Simulation Protocol

This directory contains the **Martini 3 coarse-grained molecular dynamics setup** for simulating thermal phase separation and coacervation of Elastin-Like Polypeptides (ELPs).

---

## 1. Scientific Role in the Multiscale Architecture

In the multiscale pipeline, this module serves as the **microscale input stage** that parameterizes macroscopic closures for the continuum transport engine (`02_continuum_transport`):

$$\text{Martini 3 CG-MD} \xrightarrow{\quad\text{Extract Physical Closures}\quad} \left[ R_h,\ D_0,\ \rho_p,\ B_2,\ C_{\mathrm{gel}},\ \mu(C) \right] \xrightarrow{\quad\text{Bridge JSON}\quad} \text{Continuum TFF Engine}$$

---

## 2. Directory Contents

* **`topologies/`:**
  * `martini_v3.0.0.itp`: Core Martini 3 coarse-grained force-field definition.
  * `elp_vpgvg.itp`: Coarse-grained topology for $(VPGVG)_n$ repeats with calibrated bead types (Val = P2, Pro = P1, Gly = SP2).
  * `build_elp_topology.py`: Procedural generator for custom ELP sequence block topologies.
* **`mdp/`:**
  * `em.mdp`: Steepest descent energy minimization.
  * `npt_equilibration.mdp`: 300 K isothermal-isobaric equilibration (Parrinello-Rahman, $1\,\text{bar}$).
  * `quench_coacervation.mdp`: Thermal quench (325 K) driving hydrophobic unimer collapse and droplet nucleation.
  * `shear_nemd.mdp`: Non-equilibrium Lees-Edwards shearing for shear-induced droplet alignment.
* **`scripts/`:**
  * `run_md_pipeline.py`: Automated orchestration script to run `gmx grompp` and `gmx mdrun` end-to-end.
  * `extract_transport_params.py`: Parses GROMACS XVG outputs (`gmx gyrate`, `gmx msd`) to extract $R_g, R_h, D_0, B_2$ into `transport_params.json`.
  * `phase_diagram.py`: Numerical root-finder for Flory-Huggins binodal/spinodal coexistence curves.
  * `hofmeister_coacervation.py`: Calculates salt-induced $T_t$ phase shifts and Debye screening lengths.
  * `cmc_csc_analysis.py`: Calculates Critical Micelle (CMC) and Critical Salt Concentration (CSC) phase boundaries.
  * `droplet_coalescence.py`: Models capillary droplet fusion ($v_{\text{cap}} = \gamma/\eta$) and Ostwald ripening kinetics ($R \propto t^{1/3}$).

---

## 3. Execution Modes: Production vs. Synthetic Benchmark

* **Production GROMACS Execution (requires local/HPC GROMACS installation):**
  ```bash
  python scripts/run_md_pipeline.py --run
  ```
* **Synthetic Verification Benchmark:**
  For rapid testing of the continuum solvers without running multi-hour GPU simulations, `data/sample_md_params.json` and `sample_gyrate.xvg` provide synthetic, literature-calibrated baseline parameters ($R_g = 9.8\,\text{nm}$, $R_h = 12.65\,\text{nm}$, $D_0 = 2.05 \times 10^{-11}\,\text{m}^2/\text{s}$).
