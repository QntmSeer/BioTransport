"""Interactive Streamlit Web Application for Multiscale Bioparticle Transport & Membrane Fouling."""

import sys
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
import streamlit as st

# Add python package directory to path
sys.path.insert(0, str(Path(__file__).parent.parent / "02_continuum_transport" / "python"))

try:
    from biotransport.bridge import MdBridgeModel, ProcessSimulator
    from biotransport.visualize import plot_multiscale_results
except ImportError:
    st.error("Could not import biotransport module. Run from repository root.")

st.set_page_config(
    page_title="Multiscale Bioparticle Transport Simulator",
    page_icon="🧪",
    layout="wide",
)

st.title("🧪 Multiscale Bioparticle Transport & Membrane Fouling Simulator")
st.markdown(
    """
    **Coupling Coarse-Grained Molecular Dynamics (Martini 3 / GROMACS) to Macroscale Continuum Filtration Models.**
    Adjust microscale physical properties and macroscale process conditions to observe real-time membrane fouling dynamics.
    """
)

# Sidebar: Controls
st.sidebar.header("1. Microscale Parameters (from CG-MD)")
rh_nm = st.sidebar.slider("Hydrodynamic Radius $R_h$ (nm)", min_value=3.0, max_value=30.0, value=12.6, step=0.5)
rho_p = st.sidebar.slider("Condensate Density $\\rho_p$ ($\\text{kg/m}^3$)", min_value=1000.0, max_value=1300.0, value=1150.0, step=10.0)
d0_exp = st.sidebar.slider("Diffusivity $\\log_{10}(D_0\\ [\\text{m}^2/\\text{s}])$", min_value=-12.0, max_value=-10.0, value=-10.76, step=0.05)
d0 = 10.0 ** d0_exp
b2_exp = st.sidebar.slider("Osmotic Virial $\\log_{10}(B_2\\ [\\text{m}^3/\\text{mol}])$", min_value=-5.0, max_value=-3.0, value=-4.0, step=0.1)
b2 = 10.0 ** b2_exp
comp_n = st.sidebar.slider("Cake Compressibility Exponent $n$", min_value=0.0, max_value=0.9, value=0.45, step=0.05)

st.sidebar.header("2. Process Operating Conditions")
tmp_bar = st.sidebar.slider("Transmembrane Pressure TMP (bar)", min_value=0.2, max_value=3.5, value=1.5, step=0.1)
tmp_pa = tmp_bar * 100_000.0
shear_rate = st.sidebar.slider("Crossflow Shear Rate $\\dot{\\gamma}$ ($\\text{s}^{-1}$)", min_value=500, max_value=15000, value=4000, step=500)
bulk_conc = st.sidebar.slider("Feed Bulk Concentration $C_b$ (g/L)", min_value=1.0, max_value=50.0, value=10.0, step=1.0)
total_time_min = st.sidebar.slider("Filtration Duration (minutes)", min_value=10, max_value=180, value=60, step=10)

# Build dynamic model instance
md_dict = {
    "metadata": {"source": "Streamlit Interactive UI", "model": "ELP_Dynamic"},
    "thermodynamics": {
        "temperature_K": 315.15,
        "molecular_weight_g_mol": 45000.0,
        "transition_temperature_Tt_K": 308.15,
    },
    "microscale_properties": {
        "radius_of_gyration_Rg_nm": rh_nm / 1.28,
        "hydrodynamic_radius_Rh_nm": rh_nm,
        "particle_density_kg_m3": rho_p,
        "diffusion_coefficient_D0_m2_s": d0,
        "osmotic_virial_B2_m3_mol": b2,
        "compressibility_exponent_n": comp_n,
        "gel_concentration_g_L": 420.0,
    },
}

md_model = MdBridgeModel(**md_dict)
sim = ProcessSimulator(md_model)

# Run Simulation
results = sim.simulate_filtration(
    tmp_pa=tmp_pa,
    bulk_conc_g_l=bulk_conc,
    shear_rate_s_inv=float(shear_rate),
    total_time_s=float(total_time_min * 60),
    n_steps=100,
)

# Key Metrics Display
col1, col2, col3, col4 = st.columns(4)
col1.metric("Initial Flux", f"{results['initial_flux_lmh']:.1f} LMH")
col2.metric("Final Flux", f"{results['final_flux_lmh']:.1f} LMH", delta=f"-{results['flux_decline_percent']:.1f}%")
col3.metric("Permeate Yield", f"{results['permeate_volume_l_m2'][-1]:.2f} L/m²")
col4.metric("Max Wall Concentration", f"{max(results['wall_conc_g_l']):.1f} g/L")

# Plots
st.subheader("Dynamic Fouling & Transport Profiles")

t_min = np.array(results["time_s"]) / 60.0
fig, axes = plt.subplots(1, 3, figsize=(16, 4.5), dpi=150)

# Panel 1: Flux Decline
axes[0].plot(t_min, results["flux_lmh"], color="#1f77b4", lw=2.5, label="Flux $J(t)$")
axes[0].fill_between(t_min, results["flux_lmh"], color="#1f77b4", alpha=0.15)
axes[0].set_xlabel("Time (min)")
axes[0].set_ylabel("Permeate Flux (LMH)")
axes[0].set_title("Permeate Flux Decline")
axes[0].grid(True, linestyle="--", alpha=0.5)

# Panel 2: Concentration Polarization
axes[1].plot(t_min, results["wall_conc_g_l"], color="#d62728", lw=2.5, label="Wall Conc $C_w(t)$")
axes[1].axhline(y=420.0, color="darkred", linestyle=":", lw=1.8, label="Gel Limit")
axes[1].set_xlabel("Time (min)")
axes[1].set_ylabel("Wall Concentration (g/L)")
axes[1].set_title("Concentration Polarization")
axes[1].grid(True, linestyle="--", alpha=0.5)
axes[1].legend()

# Panel 3: Cake Resistance
axes[2].plot(t_min, np.array(results["cake_resistance_m_inv"]) / 1e12, color="#2ca02c", lw=2.5, label="Cake $R_c(t)$")
axes[2].set_xlabel("Time (min)")
axes[2].set_ylabel("Cake Resistance ($10^{12}\\ \\mathrm{m}^{-1}$)")
axes[2].set_title("Cake Layer Resistance Growth")
axes[2].grid(True, linestyle="--", alpha=0.5)

st.pyplot(fig)
