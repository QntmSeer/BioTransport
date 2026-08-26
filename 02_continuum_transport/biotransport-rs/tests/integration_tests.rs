//! Integration tests for biotransport-rs.

use biotransport::{
    compute_compressible_cake_resistance, compute_uncompressed_specific_cake_resistance,
    compute_virial_osmotic_pressure, simulate_tff_filtration, solve_boundary_layer_profile,
    MdBridgeParameters, MicroscaleProperties, ProcessOperatingConditions, Thermodynamics,
};

fn get_test_md_params() -> MdBridgeParameters {
    MdBridgeParameters {
        metadata: None,
        thermodynamics: Thermodynamics {
            temperature_k: 315.15,
            molecular_weight_g_mol: 45000.0,
            transition_temperature_tt_k: Some(308.15),
        },
        microscale_properties: MicroscaleProperties {
            radius_of_gyration_rg_nm: 9.8,
            hydrodynamic_radius_rh_nm: 12.6,
            particle_density_kg_m3: 1150.0,
            diffusion_coefficient_d0_m2_s: 1.74e-11,
            osmotic_virial_b2_m3_mol: 1.2e-4,
            compressibility_exponent_n: 0.45,
            gel_concentration_g_l: 400.0,
        },
    }
}

#[test]
fn test_end_to_end_multiscale_bridge() {
    let md = get_test_md_params();
    let ops = ProcessOperatingConditions {
        transmembrane_pressure_pa: 200_000.0, // 2 bar
        membrane_resistance_m_inv: 1.0e12,
        crossflow_shear_rate_s_inv: 5000.0,
        bulk_concentration_g_l: 15.0,
        total_time_s: 3600.0,
        time_steps: 50,
        cake_porosity: 0.42,
    };

    let summary = simulate_tff_filtration(&md, &ops);

    assert!(summary.initial_flux_lmh > 0.0);
    assert!(summary.final_flux_lmh > 0.0);
    assert!(summary.initial_flux_lmh >= summary.final_flux_lmh);
    assert!(summary.total_permeate_collected_l_m2 > 0.0);
    assert!(summary.maximum_wall_concentration_g_l >= ops.bulk_concentration_g_l);
    assert!(summary.specific_cake_resistance_m_kg > 1.0e12);
}

#[test]
fn test_boundary_layer_monotonicity() {
    let profile = solve_boundary_layer_profile(10.0, 1.5e-5, 2.5e-5, 2.0e-11, 20);
    for window in profile.concentration_g_l.windows(2) {
        assert!(window[0] >= window[1]); // Decreasing from wall to bulk
    }
}
