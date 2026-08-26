//! Type definitions and JSON data contract schemas.

use serde::{Deserialize, Serialize};

/// Input parameters bridged directly from coarse-grained molecular dynamics.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct MdBridgeParameters {
    pub metadata: Option<Metadata>,
    pub thermodynamics: Thermodynamics,
    pub microscale_properties: MicroscaleProperties,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Metadata {
    pub source: String,
    pub model: String,
    pub solvent: String,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Thermodynamics {
    #[serde(rename = "temperature_K")]
    pub temperature_k: f64,
    #[serde(rename = "molecular_weight_g_mol")]
    pub molecular_weight_g_mol: f64,
    #[serde(rename = "transition_temperature_Tt_K")]
    pub transition_temperature_tt_k: Option<f64>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct MicroscaleProperties {
    #[serde(rename = "radius_of_gyration_Rg_nm")]
    pub radius_of_gyration_rg_nm: f64,
    #[serde(rename = "hydrodynamic_radius_Rh_nm")]
    pub hydrodynamic_radius_rh_nm: f64,
    #[serde(rename = "particle_density_kg_m3")]
    pub particle_density_kg_m3: f64,
    #[serde(rename = "diffusion_coefficient_D0_m2_s")]
    pub diffusion_coefficient_d0_m2_s: f64,
    #[serde(rename = "osmotic_virial_B2_m3_mol")]
    pub osmotic_virial_b2_m3_mol: f64,
    #[serde(rename = "compressibility_exponent_n", default = "default_compressibility")]
    pub compressibility_exponent_n: f64,
    #[serde(rename = "gel_concentration_g_L", default = "default_gel_concentration")]
    pub gel_concentration_g_l: f64,
}

fn default_compressibility() -> f64 {
    0.45
}

fn default_gel_concentration() -> f64 {
    400.0
}

/// Operational parameters for the membrane filtration / TFF process.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ProcessOperatingConditions {
    /// Transmembrane pressure in Pascals (e.g. 100,000 to 300,000 Pa).
    pub transmembrane_pressure_pa: f64,
    /// Clean membrane hydraulic resistance in m^-1 (e.g. 1.0e12 m^-1).
    pub membrane_resistance_m_inv: f64,
    /// Crossflow shear rate in s^-1 (e.g. 2000 to 10000 s^-1).
    pub crossflow_shear_rate_s_inv: f64,
    /// Feed bulk concentration in g/L (kg/m^3).
    pub bulk_concentration_g_l: f64,
    /// Total filtration time in seconds.
    pub total_time_s: f64,
    /// Number of output time steps.
    pub time_steps: usize,
    /// Porosity of the deposited cake layer (dimensionless, 0.35 to 0.55).
    pub cake_porosity: f64,
}

impl Default for ProcessOperatingConditions {
    fn default() -> Self {
        Self {
            transmembrane_pressure_pa: 150_000.0, // 1.5 bar
            membrane_resistance_m_inv: 1.0e12,
            crossflow_shear_rate_s_inv: 4000.0,
            bulk_concentration_g_l: 10.0,
            total_time_s: 3600.0, // 1 hour
            time_steps: 100,
            cake_porosity: 0.40,
        }
    }
}

/// Simulation output results at a given time step.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TimePointResult {
    pub time_s: f64,
    pub permeate_flux_m_s: f64,
    pub permeate_flux_lmh: f64, // L / (m^2 * h)
    pub wall_concentration_g_l: f64,
    pub cake_resistance_m_inv: f64,
    pub osmotic_pressure_pa: f64,
    pub accumulated_permeate_volume_m3_m2: f64,
}

/// Overall summary of the TFF / diafiltration process run.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ProcessSimulationSummary {
    pub initial_flux_lmh: f64,
    pub final_flux_lmh: f64,
    pub flux_decline_percentage: f64,
    pub total_permeate_collected_l_m2: f64,
    pub maximum_wall_concentration_g_l: f64,
    pub specific_cake_resistance_m_kg: f64,
    pub time_series: Vec<TimePointResult>,
}
