//! Carman-Kozeny hydrodynamic cake resistance and compressibility closures.

/// Computes the uncompressed specific cake resistance r_c0 (m / kg)
/// using the classic Carman-Kozeny equation:
///
/// r_c0 = 180 * (1 - epsilon)^2 / (rho_p * d_p^2 * epsilon^3)
///
/// where:
/// - epsilon is the cake void fraction / porosity (0 < epsilon < 1)
/// - rho_p is particle density in kg/m^3
/// - d_p is particle diameter in meters (2 * Rh)
pub fn compute_uncompressed_specific_cake_resistance(
    rh_nm: f64,
    rho_p_kg_m3: f64,
    porosity: f64,
) -> f64 {
    let dp_m = 2.0 * rh_nm * 1e-9;
    let eps = porosity.clamp(0.05, 0.95);
    let numerator = 180.0 * (1.0 - eps).powi(2);
    let denominator = rho_p_kg_m3 * dp_m.powi(2) * eps.powi(3);
    numerator / denominator
}

/// Computes the pressure-dependent compressible specific cake resistance:
///
/// r_c(Delta P) = r_c0 * (Delta P / P_ref)^n
///
/// where:
/// - Delta P is transmembrane pressure (Pa)
/// - P_ref is standard reference pressure (100 kPa = 1 bar)
/// - n is the compressibility index (n = 0 for rigid beads, n > 0 for soft condensates/vesicles)
pub fn compute_compressible_cake_resistance(
    rc0: f64,
    transmembrane_pressure_pa: f64,
    compressibility_exponent: f64,
) -> f64 {
    let p_ref = 100_000.0; // 1 bar reference pressure
    let p_ratio = (transmembrane_pressure_pa / p_ref).max(0.1);
    rc0 * p_ratio.powf(compressibility_exponent)
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn test_carman_kozeny_scaling() {
        let rc0 = compute_uncompressed_specific_cake_resistance(12.5, 1150.0, 0.40);
        assert!(rc0 > 1.0e13);
        assert!(rc0 < 1.0e16);

        // Larger particles should yield lower specific resistance (r_c ~ 1/dp^2)
        let rc0_larger = compute_uncompressed_specific_cake_resistance(25.0, 1150.0, 0.40);
        assert_relative_eq!(rc0 / rc0_larger, 4.0, epsilon = 1e-3);
    }

    #[test]
    fn test_compressibility() {
        let rc0 = 1.0e14;
        let rc_compressed = compute_compressible_cake_resistance(rc0, 200_000.0, 0.5);
        // (200k / 100k)^0.5 = sqrt(2) ≈ 1.4142
        assert_relative_eq!(rc_compressed / rc0, std::f64::consts::SQRT_2, epsilon = 1e-4);
    }
}
