"""Data contract schema and reference numerical solver for bioparticle transport."""

from __future__ import annotations
import json
import math
from pathlib import Path
from typing import Any, Optional
import numpy as np
from pydantic import BaseModel, Field


class Thermodynamics(BaseModel):
    temperature_K: float = Field(default=315.15, description="Temperature in Kelvin")
    molecular_weight_g_mol: float = Field(default=45000.0, description="Molar mass in g/mol")
    transition_temperature_Tt_K: Optional[float] = None


class MicroscaleProperties(BaseModel):
    radius_of_gyration_Rg_nm: float
    hydrodynamic_radius_Rh_nm: float
    particle_density_kg_m3: float = Field(default=1150.0)
    diffusion_coefficient_D0_m2_s: float
    osmotic_virial_B2_m3_mol: float = Field(default=1.0e-4)
    compressibility_exponent_n: float = Field(default=0.45)
    gel_concentration_g_L: float = Field(default=400.0)


class MdBridgeModel(BaseModel):
    metadata: Optional[dict[str, Any]] = None
    thermodynamics: Thermodynamics
    microscale_properties: MicroscaleProperties

    @classmethod
    def load_json(cls, path: Path | str) -> MdBridgeModel:
        with open(path, "r", encoding="utf-8") as f:
            data = json.load(f)
        return cls(**data)


class ProcessSimulator:
    """Simulates TFF filtration and membrane fouling based on MD closures."""

    def __init__(self, md_data: MdBridgeModel):
        self.md = md_data
        self.props = md_data.microscale_properties
        self.thermo = md_data.thermodynamics

    def water_viscosity(self, temp_k: float) -> float:
        t_c = temp_k - 273.15
        return float(2.414e-5 * (10.0 ** (247.8 / (t_c + 140.0))))

    def specific_cake_resistance(self, tmp_pa: float, porosity: float = 0.40) -> float:
        dp_m = 2.0 * self.props.hydrodynamic_radius_Rh_nm * 1e-9
        rho = self.props.particle_density_kg_m3
        eps = np.clip(porosity, 0.1, 0.9)
        rc0 = 180.0 * ((1.0 - eps) ** 2) / (rho * (dp_m ** 2) * (eps ** 3))
        p_ratio = max(0.1, tmp_pa / 100_000.0)
        return float(rc0 * (p_ratio ** self.props.compressibility_exponent_n))

    def virial_osmotic_pressure(self, c_g_l: float) -> float:
        r_gas = 8.314462618
        mw_kg_mol = self.thermo.molecular_weight_g_mol * 1e-3
        t = self.thermo.temperature_K
        b2 = self.props.osmotic_virial_B2_m3_mol
        term1 = (r_gas * t / mw_kg_mol) * c_g_l
        term2 = (r_gas * t * b2 / (mw_kg_mol ** 2)) * (c_g_l ** 2)
        return float(max(0.0, term1 + term2))

    def mass_transfer_coefficient(self, shear_rate_s_inv: float, channel_length_m: float = 0.20) -> float:
        d0 = max(1e-14, self.props.diffusion_coefficient_D0_m2_s)
        return float(0.816 * ((shear_rate_s_inv * (d0 ** 2) / channel_length_m) ** (1.0 / 3.0)))

    def simulate_filtration(
        self,
        tmp_pa: float = 150_000.0,
        bulk_conc_g_l: float = 10.0,
        membrane_rm: float = 1.0e12,
        shear_rate_s_inv: float = 4000.0,
        total_time_s: float = 3600.0,
        n_steps: int = 100,
        porosity: float = 0.40,
    ) -> dict[str, Any]:
        dt = total_time_s / max(1, n_steps)
        mu_0 = self.water_viscosity(self.thermo.temperature_K)
        km = self.mass_transfer_coefficient(shear_rate_s_inv)
        rc_spec = self.specific_cake_resistance(tmp_pa, porosity)
        c_gel = self.props.gel_concentration_g_L

        time_vec = np.linspace(0, total_time_s, n_steps + 1)
        flux_vec = []
        cw_vec = []
        rc_vec = []
        vol_vec = []

        m_cake = 0.0
        total_vol = 0.0
        # Physical osmotic ceiling where Pi(C) = 0.88 * TMP (prevents unphysical shutoff)
        r_gas = 8.314462
        temp_k = self.thermo.temperature_K
        mw = self.thermo.molecular_weight_g_mol / 1000.0
        b2 = self.props.osmotic_virial_B2_m3_mol
        a_coef = (r_gas * temp_k * b2 / (mw**2))
        b_coef = (r_gas * temp_k / mw)
        c_coef = -0.88 * tmp_pa
        c_osm_max = (-b_coef + np.sqrt(max(0.0, b_coef**2 - 4 * a_coef * c_coef))) / (2 * max(1e-12, a_coef))
        c_max_allowed = min(c_gel, float(c_osm_max))

        c_w = bulk_conc_g_l

        for t in time_vec:
            r_cake = rc_spec * m_cake
            total_r = membrane_rm + r_cake

            phi_rel = np.clip(c_w / (c_gel * 1.2), 0.0, 0.95)
            mu_w = mu_0 * ((1.0 - phi_rel) ** (-1.2))
            pi_w = self.virial_osmotic_pressure(c_w)
            eff_tmp = max(1000.0, tmp_pa - pi_w)
            j_m_s = eff_tmp / (mu_w * total_r)

            j_lmh = j_m_s * 3.6e6
            flux_vec.append(j_lmh)
            cw_vec.append(c_w)
            rc_vec.append(r_cake)
            vol_vec.append(total_vol * 1000.0)

            # Unsteady polarization evolution: relaxation toward film-theory target C_b*exp(J/km)
            pe = np.clip(j_m_s / km, 0.0, 3.8)
            c_w_target = min(c_max_allowed, float(bulk_conc_g_l * np.exp(pe)))
            tau_cp = 180.0  # s (boundary layer polarization relaxation time constant)
            c_w += (c_w_target - c_w) * (1.0 - np.exp(-dt / tau_cp))

            # Dynamic cake layer mass accumulation in crossflow
            j_crit = 0.30 * km
            excess_flux = max(0.0, j_m_s - j_crit)
            sticking_prob = 0.015 * ((c_w / c_gel) ** 1.5)
            deposit_rate = sticking_prob * excess_flux * bulk_conc_g_l
            shear_erosion = (shear_rate_s_inv / 4000.0) * 1.5e-3 * m_cake
            m_cake += max(0.0, deposit_rate - shear_erosion) * dt
            total_vol += j_m_s * dt

        return {
            "time_s": time_vec.tolist(),
            "flux_lmh": flux_vec,
            "wall_conc_g_l": cw_vec,
            "cake_resistance_m_inv": rc_vec,
            "permeate_volume_l_m2": vol_vec,
            "initial_flux_lmh": flux_vec[0],
            "final_flux_lmh": flux_vec[-1],
            "flux_decline_percent": float(100.0 * (flux_vec[0] - flux_vec[-1]) / max(1e-5, flux_vec[0])),
            "flux_decline_pct": float(100.0 * (flux_vec[0] - flux_vec[-1]) / max(1e-5, flux_vec[0])),
            "total_permeate_l_m2": vol_vec[-1],
            "max_wall_conc_g_l": max(cw_vec),
            "specific_cake_resistance": rc_spec,
        }


def run_continuum_simulation(
    params_json: Path | str,
    tmp_pa: float = 150_000.0,
    total_time_s: float = 3600.0,
    shear_rate_s_inv: float = 4000.0,
    bulk_conc_g_l: float = 5.0,
    membrane_rm: float = 1.0e12,
    n_steps: int = 100,
    porosity: float = 0.40,
) -> dict[str, Any]:
    md = MdBridgeModel.load_json(params_json)
    sim = ProcessSimulator(md)
    return sim.simulate_filtration(
        tmp_pa=tmp_pa,
        total_time_s=total_time_s,
        shear_rate_s_inv=shear_rate_s_inv,
        bulk_conc_g_l=bulk_conc_g_l,
        membrane_rm=membrane_rm,
        n_steps=n_steps,
        porosity=porosity,
    )
