import json
import sys
from pathlib import Path
import pytest

sys.path.insert(0, str(Path(__file__).parent.parent))
from biotransport.bridge import MdBridgeModel, ProcessSimulator, run_continuum_simulation


@pytest.fixture
def sample_json(tmp_path: Path) -> Path:
    data = {
        "metadata": {"source": "Test", "model": "ELP_Test", "solvent": "Water"},
        "thermodynamics": {
            "temperature_K": 315.15,
            "molecular_weight_g_mol": 45000.0,
            "transition_temperature_Tt_K": 308.15,
        },
        "microscale_properties": {
            "radius_of_gyration_Rg_nm": 9.8,
            "hydrodynamic_radius_Rh_nm": 12.6,
            "particle_density_kg_m3": 1150.0,
            "diffusion_coefficient_D0_m2_s": 1.74e-11,
            "osmotic_virial_B2_m3_mol": 1.2e-4,
            "compressibility_exponent_n": 0.45,
            "gel_concentration_g_L": 400.0,
        },
    }
    file_p = tmp_path / "test_params.json"
    with open(file_p, "w", encoding="utf-8") as f:
        json.dump(data, f)
    return file_p


def test_load_md_bridge_schema(sample_json: Path):
    model = MdBridgeModel.load_json(sample_json)
    assert model.thermodynamics.temperature_K == 315.15
    assert model.microscale_properties.hydrodynamic_radius_Rh_nm == 12.6


def test_process_simulation_runs(sample_json: Path):
    results = run_continuum_simulation(sample_json, tmp_pa=150_000.0)
    assert "flux_lmh" in results
    assert "wall_conc_g_l" in results
    assert len(results["flux_lmh"]) > 0
    assert results["initial_flux_lmh"] >= results["final_flux_lmh"]
    assert results["flux_decline_percent"] >= 0.0
