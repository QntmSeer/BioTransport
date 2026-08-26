//! 1D Finite-Difference Boundary Layer Convection-Diffusion PDE Solver.

/// Result of a 1D concentration polarization boundary layer discretization.
#[derive(Debug, Clone)]
pub struct BoundaryLayerProfile {
    pub y_grid_m: Vec<f64>,
    pub concentration_g_l: Vec<f64>,
    pub wall_concentration_g_l: f64,
    pub boundary_layer_thickness_m: f64,
}

/// Solves the steady-state and quasi-transient concentration polarization boundary layer profile.
///
/// Governing ODE across boundary layer (0 <= y <= delta):
/// D * (d^2 C / dy^2) + J * (dC / dy) = 0
///
/// Analytical solution (Film Theory):
/// C(y) = C_b * exp( (J / D) * (delta - y) )
/// C_w = C(0) = C_b * exp( J * delta / D ) = C_b * exp( J / k_m )
///
/// where k_m = D / delta is the mass transfer coefficient.
pub fn solve_boundary_layer_profile(
    bulk_conc_g_l: f64,
    permeate_flux_m_s: f64,
    mass_transfer_coeff_k_m_s: f64,
    diffusivity_d0_m2_s: f64,
    n_points: usize,
) -> BoundaryLayerProfile {
    let delta = diffusivity_d0_m2_s / mass_transfer_coeff_k_m_s.max(1e-9);
    let pe = (permeate_flux_m_s / mass_transfer_coeff_k_m_s.max(1e-9)).clamp(0.0, 10.0);
    let c_w = bulk_conc_g_l * pe.exp();

    let mut y_grid = Vec::with_capacity(n_points);
    let mut c_profile = Vec::with_capacity(n_points);

    for i in 0..n_points {
        let frac = i as f64 / (n_points - 1) as f64;
        let y = frac * delta;
        // Concentration decreases exponentially from wall (y=0) to bulk (y=delta)
        let c_y = bulk_conc_g_l * ((permeate_flux_m_s / diffusivity_d0_m2_s) * (delta - y)).exp();
        let c_clamped = c_y.clamp(bulk_conc_g_l, c_w);
        y_grid.push(y);
        c_profile.push(c_clamped);
    }

    BoundaryLayerProfile {
        y_grid_m: y_grid,
        concentration_g_l: c_profile,
        wall_concentration_g_l: c_w,
        boundary_layer_thickness_m: delta,
    }
}

/// Computes the mass transfer coefficient k_m (m/s) from the Leveque correlation for laminar channel flow:
///
/// Sh = 0.816 * (Re * Sc * dh / L)^(1/3)
/// or in terms of wall shear rate gamma_dot:
/// k_m = 0.816 * (gamma_dot * D^2 / L)^(1/3)
pub fn compute_mass_transfer_coefficient(
    shear_rate_s_inv: f64,
    diffusivity_d0_m2_s: f64,
    channel_length_m: f64,
) -> f64 {
    let gamma = shear_rate_s_inv.max(100.0);
    let d = diffusivity_d0_m2_s.max(1e-14);
    let l = channel_length_m.max(0.01);

    0.816 * (gamma * d.powi(2) / l).cbrt()
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn test_boundary_layer_film_theory_asymptote() {
        let bulk = 10.0;
        let flux = 1.0e-5; // 10 um/s = 36 LMH
        let km = 2.0e-5;
        let d0 = 2.0e-11;

        let profile = solve_boundary_layer_profile(bulk, flux, km, d0, 50);
        // C_w / C_b = exp(1.0e-5 / 2.0e-5) = exp(0.5) ≈ 1.6487
        let expected_cw = bulk * (0.5_f64).exp();
        assert_relative_eq!(profile.wall_concentration_g_l, expected_cw, epsilon = 1e-3);
        assert_eq!(profile.concentration_g_l.len(), 50);
        assert_relative_eq!(profile.concentration_g_l[49], bulk, epsilon = 1e-3);
    }
}
