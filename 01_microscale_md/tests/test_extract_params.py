import math
import sys
from pathlib import Path
import pytest

# Add scripts directory to path
sys.path.insert(0, str(Path(__file__).parent.parent / "scripts"))
from extract_transport_params import (
    compute_hydrodynamic_radius_from_rg,
    compute_stokes_einstein_diffusion,
    extract_parameters_from_synthetic_or_traj,
)


def test_hydrodynamic_radius_scaling():
    rg = 10.0
    rh_globular = compute_hydrodynamic_radius_from_rg(rg, "globular_condensate")
    assert rh_globular == pytest.approx(math.sqrt(5.0 / 3.0) * 10.0, rel=1e-3)
    assert rh_globular > rg

    rh_coil = compute_hydrodynamic_radius_from_rg(rg, "random_coil")
    assert rh_coil < rg


def test_stokes_einstein_diffusion():
    # For a 10 nm radius particle in water (0.001 Pa.s) at 300 K:
    # D0 = (1.38e-23 * 300) / (6 * pi * 0.001 * 10e-9) ≈ 2.196e-11 m^2/s
    d0 = compute_stokes_einstein_diffusion(rh_nm=10.0, temperature_k=300.0, viscosity_pa_s=1e-3)
    assert d0 == pytest.approx(2.196e-11, rel=1e-2)


def test_extract_parameters_schema():
    params = extract_parameters_from_synthetic_or_traj(
        temperature_k=315.15,
        measured_rg_nm=12.0,
    )
    assert "thermodynamics" in params
    assert "microscale_properties" in params
    props = params["microscale_properties"]
    assert props["hydrodynamic_radius_Rh_nm"] > props["radius_of_gyration_Rg_nm"]
    assert 1000.0 <= props["particle_density_kg_m3"] <= 1300.0
    assert props["diffusion_coefficient_D0_m2_s"] > 0.0
    assert props["gel_concentration_g_L"] > 0.0
