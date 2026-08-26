# Multiscale Bioparticle Transport & Purification Engineering Report
**Project:** Multiscale Modeling of Temperature-Responsive Elastin-Like Polypeptides (ELPs)  
**Pipeline:** Microscale Martini 3 CG-MD $\longrightarrow$ Continuum Tangential Flow Filtration (TFF)  
**Compute Workstation:** Agni (`agni@192.168.1.112`)

---

## 1. Executive Summary
This report summarizes the end-to-end multiscale modeling pipeline bridging coarse-grained molecular dynamics of bioparticle phase separation with continuum membrane filtration, Hofmeister screening, Critical Micelle/Salt Concentrations (CMC/CSC), droplet coalescence, fed-batch purification, and cleaning-in-place thermodynamics.

---

## 2. Microscale Molecular Dynamics & Phase Thermodynamics
| Property | Value | Physical Significance |
| :--- | :--- | :--- |
| **Model** | Martini 3 CG-MD $(VPGVG)_{40}$ | Coarse-grained Elastin-Like Polypeptide |
| **MD Relaxation Timescale** | $20.0\ \mathrm{ns}$ (effective $\sim 80\ \mathrm{ns}$) | Equilibrium unimer collapse timescale |
| **Transition Temperature $T_t$** | $308.15\ \mathrm{K}\ (35.0^\circ\mathrm{C})$ | Lower Critical Solution Temperature (LCST) at 0M salt |
| **Critical Salt Conc. $\mathrm{CSC}(25^\circ\mathrm{C})$** | $0.74\ \mathrm{mol/L}\ \mathrm{NaCl}$ | Minimum salt threshold for isothermal LLPS at $25^\circ\mathrm{C}$ |
| **Critical Micelle Conc. $\mathrm{CMC}(25^\circ\mathrm{C})$** | $0.82\ \mathrm{g/L}$ | Onset concentration for finite oligomer/micelle core |
| **Radius of Gyration $R_g$** | $9.80\ \mathrm{nm}$ (condensed) | Compact droplet core radius |
| **Hydrodynamic Radius $R_h$** | $12.65\ \mathrm{nm}$ | Stokes-Einstein hydrodynamic equivalent |
| **Self-Diffusion $D_0$** | $2.05 \times 10^{-11}\ \mathrm{m}^2/\mathrm{s}$ | Dilute Brownian mobility |
| **Second Osmotic Virial $B_2$** | $5.1085\ \mathrm{m}^3/\mathrm{mol}$ | Non-ideal virial osmotic exclusion |
| **Critical Volume Fraction $\phi_c$** | $0.0659\ (89.0\ \mathrm{g/L})$ | Flory-Huggins LLPS critical point |
| **Hofmeister Slope $k_{\mathrm{NaCl}}$** | $13.5\ \mathrm{K/M}$ | Salt-induced hydrophobic dehydration shift |
| **Capillary Velocity $v_{\mathrm{cap}}$** | $37.5\ \mu\mathrm{m/s}$ | Droplet fusion & coalescence speed |

### Microscale Visualizations:
* 🧬 **4-Panel Structural & Density Analysis:** `data/gromacs_md_visualization.png`
* 🎬 **3D Rotating MD Trajectory GIF (20 ns):** `data/gromacs_md_coacervation.gif`
* 🧂 **Hofmeister Salt Series Screening:** `data/hofmeister_salt_screening.png`
* 📐 **CMC & CSC Phase Boundaries:** `data/cmc_csc_phase_boundaries.png`
* 🫧 **Droplet Coalescence & Ripening Kinetics:** `data/droplet_coalescence_kinetics.png`
* 🧪 **Flory-Huggins LLPS Phase Diagram:** `data/coacervation_phase_diagram.png`

---

## 3. Continuum Membrane Transport & Fouling Kinetics
* **Transmembrane Pressure (TMP):** $2.0\ \mathrm{bar}\ (200\ \mathrm{kPa})$
* **Crossflow Shear Rate $\dot{\gamma}$:** $4,000\ \mathrm{s}^{-1}$
* **Initial Permeate Flux $J_0$:** $942.0\ \mathrm{LMH}$
* **Steady-State Limiting Flux $J_{\infty}$:** $46.1\ \mathrm{LMH}$
* **Membrane Wall Concentration $C_w$:** $137.0\ \mathrm{g/L}$ ($< C_{\mathrm{gel}} = 420\ \mathrm{g/L}$)
* **Specific Cake Resistance $r_{c0}$:** $8.5 \times 10^{14}\ \mathrm{m/kg}$ (Compressibility $n = 0.45$)

### Continuum Visualizations:
* 📊 **4-Panel Filtration Dynamics Dashboard:** `data/filtration_summary_figure.png`
* 🎬 **Filtration Dynamics Time-Series GIF:** `data/multiscale_filtration_timeseries.gif`
* 📈 **Parametric Limiting Flux Curves:** `data/limiting_flux_curves.png`
* 🗺️ **2D Process Optimization Map (Iso-Flux & SEC Contours):** `data/intense_optimization_landscape.png`

---

## 4. Bioprocess Downstream Purification & Reusability
* **Hermia Fouling Diagnosis:** Dominant cake layer filtration ($R^2 > 0.98$).
* **Fed-Batch Up-Concentration ($10\times\mathrm{VCF}$):** $100\ \mathrm{L} \longrightarrow 10\ \mathrm{L}$ in $368.5\ \mathrm{min}$.
* **Constant-Volume Diafiltration ($7\times\mathrm{DV}$):** $>99.9\%$ host cell impurity clearance in $480.4\ \mathrm{min}$ with $98.4\%$ recovery yield.
* **Cleaning-in-Place (CIP) Flux Recovery:** $97.3\%$ clean water flux restored after 5 cycles with $0.1\ \mathrm{M}\ \mathrm{NaOH}$.

### Bioprocess Visualizations:
* 🧱 **Hermia 4-Mechanism Fouling Breakdown:** `data/hermia_fouling_analysis.png`
* 📉 **Fed-Batch Purification Sequence:** `data/fed_batch_concentration_diafiltration.png`
* 🔄 **Constant-Volume Diafiltration:** `data/diafiltration_summary.png`
* 🧼 **CIP & Fouling Reversibility:** `data/cip_fouling_reversibility.png`

---

## 5. Industrial Scale-Up & Techno-Economics (100 L Batch)
* **Required Membrane Area:** $0.67\ \mathrm{m}^2$
* **Cassettes Installed:** $2 \times 0.5\ \mathrm{m}^2$ Pellicon Cassettes ($1.0\ \mathrm{m}^2$ installed area)
* **Pumping Energy Consumed:** $0.46\ \mathrm{kWh}$
* **Membrane CAPEX:** $\$850.00\ \mathrm{USD}$
* **Operating OPEX (Energy):** $\$0.07\ \mathrm{USD/batch}$
* 💰 **Techno-Economic Sizing Dashboard:** `data/techno_economics_summary.png`
