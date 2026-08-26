//! Osmotic pressure equations of state (Van 't Hoff, Second Virial, and Carnahan-Starling).

/// Computes osmotic pressure (Pa) using the second virial expansion:
///
/// Pi(C) = (R * T / M_w) * C + (R * T * B_2 / M_w^2) * C^2
///
/// where:
/// - C is mass concentration in g/L = kg/m^3
/// - M_w is molecular weight in g/mol = kg/kmol (converted to kg/mol: M_w * 1e-3)
/// - B_2 is osmotic second virial coefficient in m^3/mol
/// - R is ideal gas constant (8.314462 J / (mol * K))
pub fn compute_virial_osmotic_pressure(
    c_g_l: f64,
    temperature_k: f64,
    mw_g_mol: f64,
    b2_m3_mol: f64,
) -> f64 {
    let r_gas = 8.314462618;
    let mw_kg_mol = mw_g_mol * 1e-3;
    let c_kg_m3 = c_g_l.max(0.0);

    // Van 't Hoff first term: (R * T / Mw) * C
    let term1 = (r_gas * temperature_k / mw_kg_mol) * c_kg_m3;
    // Second virial term: (R * T * B2 / Mw^2) * C^2
    let term2 = (r_gas * temperature_k * b2_m3_mol / (mw_kg_mol * mw_kg_mol)) * c_kg_m3.powi(2);

    (term1 + term2).max(0.0)
}

/// Carnahan-Starling hard-sphere osmotic pressure equation:
///
/// Z = (1 + phi + phi^2 - phi^3) / (1 - phi)^3
/// Pi = n_v * k_B * T * Z
pub fn compute_carnahan_starling_osmotic_pressure(
    volume_fraction_phi: f64,
    temperature_k: f64,
    rh_nm: f64,
) -> f64 {
    let k_b = 1.380649e-23;
    let phi = volume_fraction_phi.clamp(0.0, 0.64); // random close packing limit
    let z = (1.0 + phi + phi.powi(2) - phi.powi(3)) / (1.0 - phi).powi(3);

    let rh_m = rh_nm * 1e-9;
    let v_particle = (4.0 / 3.0) * std::f64::consts::PI * rh_m.powi(3);
    let number_density = phi / v_particle;

    number_density * k_b * temperature_k * z
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_virial_osmotic_pressure_zero() {
        let pi = compute_virial_osmotic_pressure(0.0, 300.0, 45000.0, 1e-4);
        assert_eq!(pi, 0.0);
    }

    #[test]
    fn test_virial_osmotic_pressure_scaling() {
        let pi_low = compute_virial_osmotic_pressure(10.0, 300.0, 45000.0, 1e-4);
        let pi_high = compute_virial_osmotic_pressure(100.0, 300.0, 45000.0, 1e-4);
        assert!(pi_high > 10.0 * pi_low); // Nonlinear due to second virial term
    }
}
