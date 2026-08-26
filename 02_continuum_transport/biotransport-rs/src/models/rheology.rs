//! Suspension rheology and concentration-dependent viscosity models.

/// Computes suspension dynamic viscosity (Pa.s) using the Krieger-Dougherty model:
///
/// mu(C) = mu_0 * (1 - C / C_max)^(- [eta] * C_max)
///
/// where:
/// - mu_0 is solvent viscosity at temperature T (e.g. ~0.00089 Pa.s for water at 25°C)
/// - C is bioparticle concentration (g/L)
/// - C_max is maximum packing / gelation concentration (g/L)
/// - [eta] is intrinsic viscosity (approx 0.003 - 0.005 L/g for globular biopolymers)
pub fn compute_suspension_viscosity(
    c_g_l: f64,
    c_max_g_l: f64,
    solvent_viscosity_pa_s: f64,
    intrinsic_viscosity_l_g: f64,
) -> f64 {
    let c = c_g_l.max(0.0);
    let c_max = c_max_g_l.max(c + 1.0);
    let phi_rel = (c / c_max).clamp(0.0, 0.98);
    let exponent = intrinsic_viscosity_l_g * c_max;

    solvent_viscosity_pa_s * (1.0 - phi_rel).powf(-exponent)
}

/// Temperature-dependent water viscosity via Vogel-Fulcher-Tammann / standard empirical fit.
pub fn water_viscosity(temperature_k: f64) -> f64 {
    let t_c = temperature_k - 273.15;
    // Standard empirical Andrade fit for water:
    let a = 2.414e-5; // Pa.s
    let b = 247.8;
    let c = 140.0;
    a * 10.0_f64.powf(b / (t_c + c))
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn test_water_viscosity_at_20c() {
        let mu_20 = water_viscosity(293.15);
        // Approx 1.002 mPa.s at 20 C
        assert_relative_eq!(mu_20, 1.002e-3, epsilon = 5e-5);
    }

    #[test]
    fn test_suspension_viscosity_increases_with_concentration() {
        let mu_0 = 1.0e-3;
        let mu_10 = compute_suspension_viscosity(10.0, 400.0, mu_0, 0.004);
        let mu_100 = compute_suspension_viscosity(100.0, 400.0, mu_0, 0.004);
        assert!(mu_100 > mu_10);
        assert!(mu_10 > mu_0);
    }
}
