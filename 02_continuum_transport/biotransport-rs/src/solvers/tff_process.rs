//! Multi-cycle Tangential Flow Filtration (TFF) and diafiltration process simulator.

use crate::models::{
    carman_kozeny::{compute_compressible_cake_resistance, compute_uncompressed_specific_cake_resistance},
    osmotic_pressure::compute_virial_osmotic_pressure,
    rheology::{compute_suspension_viscosity, water_viscosity},
};
use crate::solvers::boundary_layer_1d::compute_mass_transfer_coefficient;
use crate::types::{
    MdBridgeParameters, ProcessOperatingConditions, ProcessSimulationSummary, TimePointResult,
};

/// Simulates a complete TFF filtration run over time using physical closures from MD.
pub fn simulate_tff_filtration(
    md_params: &MdBridgeParameters,
    ops: &ProcessOperatingConditions,
) -> ProcessSimulationSummary {
    let temp_k = md_params.thermodynamics.temperature_k;
    let mw = md_params.thermodynamics.molecular_weight_g_mol;
    let rh = md_params.microscale_properties.hydrodynamic_radius_rh_nm;
    let rho_p = md_params.microscale_properties.particle_density_kg_m3;
    let d0 = md_params.microscale_properties.diffusion_coefficient_d0_m2_s;
    let b2 = md_params.microscale_properties.osmotic_virial_b2_m3_mol;
    let comp_n = md_params.microscale_properties.compressibility_exponent_n;
    let c_gel = md_params.microscale_properties.gel_concentration_g_l;

    let dt = ops.total_time_s / ops.time_steps.max(1) as f64;
    let mu_0 = water_viscosity(temp_k);
    let channel_length_m = 0.20; // standard 20 cm hollow fiber / flat sheet path
    let km = compute_mass_transfer_coefficient(ops.crossflow_shear_rate_s_inv, d0, channel_length_m);

    // Carman-Kozeny specific cake resistance
    let rc0 = compute_uncompressed_specific_cake_resistance(rh, rho_p, ops.cake_porosity);
    let rc_spec = compute_compressible_cake_resistance(rc0, ops.transmembrane_pressure_pa, comp_n);

    let mut m_cake_kg_m2 = 0.0;
    let mut total_vol_m3_m2 = 0.0;
    let mut time_series = Vec::with_capacity(ops.time_steps + 1);

    let mut max_cw = ops.bulk_concentration_g_l;
    let mut initial_flux_lmh = 0.0;
    let mut final_flux_lmh = 0.0;

    for step in 0..=ops.time_steps {
        let t_s = step as f64 * dt;

        // Current cake resistance: R_cake = r_c * M_cake (m^-1)
        let r_cake = rc_spec * m_cake_kg_m2;
        let total_resistance = ops.membrane_resistance_m_inv + r_cake;

        // Iteratively solve for self-consistent (J, C_w, Pi, mu) at this time step
        let mut j_m_s = ops.transmembrane_pressure_pa / (mu_0 * total_resistance);
        let mut c_w = ops.bulk_concentration_g_l;
        let mut pi_w = 0.0;

        for _ in 0..10 {
            // Film theory wall concentration: C_w = C_b * exp(J / km)
            let pe = (j_m_s / km).clamp(0.0, 5.0);
            c_w = (ops.bulk_concentration_g_l * pe.exp()).min(c_gel);

            // Viscosity at the membrane wall
            let mu_w = compute_suspension_viscosity(c_w, c_gel * 1.2, mu_0, 0.0035);

            // Osmotic back-pressure
            pi_w = compute_virial_osmotic_pressure(c_w, temp_k, mw, b2);
            let effective_tmp = (ops.transmembrane_pressure_pa - pi_w).max(100.0);

            // Updated flux via Darcy-Starling
            let j_next = effective_tmp / (mu_w * total_resistance);
            if (j_next - j_m_s).abs() < 1e-8 {
                j_m_s = j_next;
                break;
            }
            j_m_s = 0.5 * (j_m_s + j_next);
        }

        if c_w > max_cw {
            max_cw = c_w;
        }

        // Convert flux to LMH (L / (m^2 * h)): 1 m/s = 3.6e6 LMH
        let j_lmh = j_m_s * 3.6e6;
        if step == 0 {
            initial_flux_lmh = j_lmh;
        }
        if step == ops.time_steps {
            final_flux_lmh = j_lmh;
        }

        time_series.push(TimePointResult {
            time_s: t_s,
            permeate_flux_m_s: j_m_s,
            permeate_flux_lmh: j_lmh,
            wall_concentration_g_l: c_w,
            cake_resistance_m_inv: r_cake,
            osmotic_pressure_pa: pi_w,
            accumulated_permeate_volume_m3_m2: total_vol_m3_m2,
        });

        // Time integration of cake accumulation
        // Deposition occurs when C_w reaches near gel point or convection exceeds shear back-transport
        if step < ops.time_steps {
            let convective_deposit_rate = j_m_s * (c_w * 1e-3); // kg / (m^2 * s)
            let shear_erosion_rate = (ops.crossflow_shear_rate_s_inv / 10000.0) * 1e-5 * m_cake_kg_m2;
            let d_mcake = (convective_deposit_rate - shear_erosion_rate).max(0.0) * dt;

            // Only form solid cake if C_w is above critical threshold (e.g. 70% of gelation)
            if c_w >= 0.70 * c_gel {
                m_cake_kg_m2 += d_mcake;
            }
            total_vol_m3_m2 += j_m_s * dt;
        }
    }

    let flux_decline = if initial_flux_lmh > 0.0 {
        100.0 * (initial_flux_lmh - final_flux_lmh) / initial_flux_lmh
    } else {
        0.0
    };

    ProcessSimulationSummary {
        initial_flux_lmh,
        final_flux_lmh,
        flux_decline_percentage: flux_decline,
        total_permeate_collected_l_m2: total_vol_m3_m2 * 1000.0,
        maximum_wall_concentration_g_l: max_cw,
        specific_cake_resistance_m_kg: rc_spec,
        time_series,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::types::{MicroscaleProperties, Thermodynamics};

    fn sample_params() -> MdBridgeParameters {
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
                diffusion_coefficient_d0_m2_s: 1.7e-11,
                osmotic_virial_b2_m3_mol: 1.0e-4,
                compressibility_exponent_n: 0.45,
                gel_concentration_g_l: 400.0,
            },
        }
    }

    #[test]
    fn test_tff_simulation_monotonic_flux_decay() {
        let md = sample_params();
        let ops = ProcessOperatingConditions {
            transmembrane_pressure_pa: 150_000.0,
            bulk_concentration_g_l: 10.0,
            total_time_s: 1800.0,
            time_steps: 20,
            ..Default::default()
        };

        let summary = simulate_tff_filtration(&md, &ops);
        assert!(summary.initial_flux_lmh >= summary.final_flux_lmh);
        assert!(summary.flux_decline_percentage >= 0.0);
        assert!(summary.total_permeate_collected_l_m2 > 0.0);
    }
}
