# 01_microscale_md: Coarse-Grained MD of Bioparticle Self-Assembly

This module implements the **Martini 3 Coarse-Grained Molecular Dynamics (CG-MD)** simulation workflow for modeling bioparticle coacervation, liquid-liquid phase separation (LLPS), and shear-induced alignment in GROMACS.

---

## 1. Simulation Protocols

| Stage | MDP File | Purpose | Key Parameters |
| :--- | :--- | :--- | :--- |
| **Energy Minimization** | `mdp/em.mdp` | Remove steric clashes | `integrator = steep`, `emtol = 100.0` |
| **NPT Equilibration** | `mdp/npt_equilibration.mdp` | Density equilibration at $T < T_t$ | `v-rescale` (300K), `Parrinello-Rahman` (1 bar), $\tau_p = 12.0\,\text{ps}$ |
| **Thermal Quench / Assembly** | `mdp/quench_coacervation.mdp` | Triggers phase separation ($T > T_t$) | `v-rescale` (325K), $dt = 20\,\text{fs}$, 500 ns production |
| **Couette Shear Flow** | `mdp/shear_nemd.mdp` | Non-equilibrium MD shear | `deform = vx 0 0 0 0 0`, shear rate $\dot{\gamma} = 10^7\,\text{s}^{-1}$ |

---

## 2. Automated Parameter Extraction Pipeline

The script `scripts/extract_transport_params.py` parses trajectory and energy output files to compute the macroscopic constitutive closures needed for continuum modeling:

1. **Radius of Gyration ($R_g$) & Hydrodynamic Radius ($R_h$):**
   $$R_h \approx 1.28 \cdot R_g \quad \text{(for globular coacervates / Kirkwood-Riseman approximation)}$$
2. **Mean Squared Displacement & Diffusivity ($D_0$):**
   $$\text{MSD}(t) = \langle |\mathbf{r}(t) - \mathbf{r}(0)|^2 \rangle \implies D_0 = \frac{1}{6} \frac{d(\text{MSD})}{dt}$$
3. **Condensate Density ($\rho_p$):**
   $$\rho_p = \frac{\sum m_i}{\frac{4}{3}\pi R_g^3}$$

The output is exported to `data/sample_md_params.json`.
