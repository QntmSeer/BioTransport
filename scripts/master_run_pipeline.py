"""Master end-to-end pipeline runner for the multiscale bioparticle transport project.

Executes all 15 unified modules:
1. Microscale CG-MD Coacervation simulation & 4-panel structure analysis.
2. Microscale 3D animated GIF generation (20 ns trajectory).
3. Hofmeister Salt Series & Ionic Screening Phase Shift.
4. Critical Micelle Concentration (CMC) & Critical Salt Concentration (CSC) Thermodynamics.
5. Droplet Coalescence & Ostwald Ripening Kinetics.
6. Flory-Huggins LLPS Phase Coexistence & Binodal Diagram.
7. Continuum transport & boundary layer dynamics simulation.
8. Continuum filtration animated time-series GIF generation.
9. Parametric Limiting Flux vs. TMP curve generation across shear rates.
10. Rigorous 2D Continuum Process Optimization & Design Space Map.
11. Hermia's 4-Mechanism Fouling Diagnostic.
12. Full Fed-Batch Concentration (10x VCF) + Diafiltration Sequence.
13. Constant-Volume Diafiltration (CVD) & Buffer Exchange kinetics.
14. Cleaning-in-Place (CIP) & Fouling Reversibility Breakdown.
15. Techno-Economics, Membrane Sizing & Engineering Summary Report.
"""

from __future__ import annotations
import sys
import time
from pathlib import Path

# Add python source trees
project_root = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(project_root / "01_microscale_md" / "scripts"))
sys.path.insert(0, str(project_root / "02_continuum_transport" / "python"))
sys.path.insert(0, str(project_root / "scripts"))

from plot_md_snapshots import generate_and_plot_md_visualization as run_md_snapshots
from animate_md_trajectory import create_md_coacervation_gif as run_md_animation
from hofmeister_coacervation import simulate_hofmeister_phase_behavior
from cmc_csc_analysis import analyze_cmc_csc_phase_boundaries
from droplet_coalescence import simulate_droplet_coalescence
from phase_diagram import flory_huggins_phase_diagram
from biotransport.bridge import MdBridgeModel, ProcessSimulator, run_continuum_simulation
from biotransport.visualize import plot_multiscale_results, plot_limiting_flux_curves
from biotransport.animate_filtration import main as run_filtration_animation
from biotransport.parameter_sweep import run_intense_parameter_sweep
from biotransport.hermia_fouling import analyze_hermia_fouling
from biotransport.fed_batch_purification import simulate_fed_batch_purification
from biotransport.diafiltration import simulate_diafiltration
from biotransport.cip_reversibility import simulate_cip_reversibility
from biotransport.techno_economics import calculate_techno_economics
from generate_multiscale_report import generate_report


def main() -> None:
    print("=" * 78)
    print("         MULTISCALE BIOPARTICLE TRANSPORT MASTER PIPELINE")
    print("   Microscale CG-MD Coacervation + Continuum TFF Purification Platform")
    print("=" * 78)
    start_total = time.time()

    data_dir = project_root / "data"
    data_dir.mkdir(parents=True, exist_ok=True)
    params_json = data_dir / "sample_md_params.json"

    # Stage 1: Microscale CG-MD Structural & Density Analysis
    print("\n[1/15] Microscale CG-MD Structural & Radial Density Profile...")
    t0 = time.time()
    run_md_snapshots()
    print(f"       -> Output: data/gromacs_md_visualization.png ({time.time() - t0:.2f}s)")

    # Stage 2: Microscale 3D Rotating Time-Series GIF (20 ns timescale)
    print("\n[2/15] Microscale 3D Coacervation Animated Trajectory GIF (20 ns)...")
    t0 = time.time()
    run_md_animation()
    print(f"       -> Output: data/gromacs_md_coacervation.gif ({time.time() - t0:.2f}s)")

    # Stage 3: Hofmeister Salt Screening & T_t Phase Shift
    print("\n[3/15] Hofmeister Salt Series & Ionic Screening Thermodynamics...")
    t0 = time.time()
    simulate_hofmeister_phase_behavior(output_png=data_dir / "hofmeister_salt_screening.png")
    print(f"       -> Output: data/hofmeister_salt_screening.png ({time.time() - t0:.2f}s)")

    # Stage 4: Critical Micelle (CMC) & Critical Salt (CSC) Phase Boundaries
    print("\n[4/15] Critical Micelle (CMC) & Critical Salt (CSC) Phase Thermodynamics...")
    t0 = time.time()
    cmc_csc_res = analyze_cmc_csc_phase_boundaries(output_png=data_dir / "cmc_csc_phase_boundaries.png")
    print(f"       -> Output: data/cmc_csc_phase_boundaries.png ({time.time() - t0:.2f}s)")
    print(f"          CSC at 25°C: {cmc_csc_res['csc_at_25C_M']:.2f} M NaCl | CMC at 25°C: {cmc_csc_res['cmc_at_25C_g_L']:.2f} g/L")

    # Stage 5: Condensate Droplet Coalescence & Ripening Kinetics
    print("\n[5/15] Condensate Droplet Coalescence & Ostwald Ripening Kinetics...")
    t0 = time.time()
    coalesce_res = simulate_droplet_coalescence(output_png=data_dir / "droplet_coalescence_kinetics.png")
    print(f"       -> Output: data/droplet_coalescence_kinetics.png ({time.time() - t0:.2f}s)")
    print(f"          v_cap: {coalesce_res['capillary_velocity_um_s']:.1f} um/s | Final R: {coalesce_res['final_r_shear_nm']:.1f} nm")

    # Stage 6: Flory-Huggins Phase Diagram
    print("\n[6/15] Flory-Huggins LLPS Binodal & Spinodal Phase Coexistence...")
    t0 = time.time()
    flory_huggins_phase_diagram(output_png=data_dir / "coacervation_phase_diagram.png")
    print(f"       -> Output: data/coacervation_phase_diagram.png ({time.time() - t0:.2f}s)")

    # Stage 7: Continuum TFF Transport Dynamics Simulation
    print("\n[7/15] Continuum Tangential Flow Filtration (TFF) Dynamics...")
    t0 = time.time()
    results = run_continuum_simulation(
        params_json=params_json,
        tmp_pa=150_000.0,
        total_time_s=3600.0,
    )
    with open(data_dir / "continuum_filtration_results.json", "w", encoding="utf-8") as f:
        import json
        json.dump(results, f, indent=2)
    plot_multiscale_results(results, data_dir / "filtration_summary_figure.png")
    print(f"       -> Output: data/filtration_summary_figure.png ({time.time() - t0:.2f}s)")
    print(f"          Initial Flux: {results['initial_flux_lmh']:.1f} LMH | Final Flux: {results['final_flux_lmh']:.1f} LMH")

    # Stage 8: Continuum Filtration Dynamics Animated GIF
    print("\n[8/15] Continuum Filtration Dynamics Time-Series GIF...")
    t0 = time.time()
    run_filtration_animation(str(params_json), str(data_dir / "multiscale_filtration_timeseries.gif"))
    print(f"       -> Output: data/multiscale_filtration_timeseries.gif ({time.time() - t0:.2f}s)")

    # Stage 9: Parametric Limiting Flux Curves
    print("\n[9/15] Parametric Limiting Flux Curves across Shear Rates...")
    t0 = time.time()
    model = MdBridgeModel.load_json(params_json)
    sim = ProcessSimulator(model)
    plot_limiting_flux_curves(sim, data_dir / "limiting_flux_curves.png")
    print(f"       -> Output: data/limiting_flux_curves.png ({time.time() - t0:.2f}s)")

    # Stage 10: Rigorous 2D Design Space Map
    print("\n[10/15] 2D Process Optimization Map with SEC Pareto Contours...")
    t0 = time.time()
    run_intense_parameter_sweep(
        params_json=str(params_json),
        output_png=str(data_dir / "intense_optimization_landscape.png"),
        sample_density=50,
    )
    print(f"       -> Output: data/intense_optimization_landscape.png ({time.time() - t0:.2f}s)")

    # Stage 11: Hermia 4-Mechanism Fouling Diagnostic
    print("\n[11/15] Hermia's 4-Mechanism Fouling Diagnostics...")
    t0 = time.time()
    hermia_scores = analyze_hermia_fouling(results, output_png=data_dir / "hermia_fouling_analysis.png")
    print(f"        -> Output: data/hermia_fouling_analysis.png ({time.time() - t0:.2f}s)")
    for mech, r2 in hermia_scores.items():
        print(f"           - {mech}: R² = {r2:.4f}")

    # Stage 12: Fed-Batch Purification Sequence (10x VCF + 7x DV)
    print("\n[12/15] Fed-Batch Concentration (10x VCF) + Diafiltration Sequence...")
    t0 = time.time()
    fed_batch_res = simulate_fed_batch_purification(
        v0_L=100.0,
        target_vcf=10.0,
        target_dv=7.0,
        output_png=data_dir / "fed_batch_concentration_diafiltration.png",
    )
    print(f"        -> Output: data/fed_batch_concentration_diafiltration.png ({time.time() - t0:.2f}s)")
    print(f"           Concentration: {fed_batch_res['concentration_duration_min']:.1f} min | Diafiltration: {fed_batch_res['diafiltration_duration_min']:.1f} min")

    # Stage 13: Constant-Volume Diafiltration (CVD)
    print("\n[13/15] Constant-Volume Diafiltration (CVD) Buffer Exchange Kinetics...")
    t0 = time.time()
    dia_res = simulate_diafiltration(v_retentate_L=10.0, target_dv=8.0, output_png=data_dir / "diafiltration_summary.png")
    print(f"        -> Output: data/diafiltration_summary.png ({time.time() - t0:.2f}s)")
    print(f"           Impurity Remaining: {dia_res['final_impurity_remaining_pct']:.2e}% | Yield: {dia_res['final_product_yield_pct']:.1f}%")

    # Stage 14: Cleaning-in-Place (CIP) & Fouling Reversibility Breakdown
    print("\n[14/15] Cleaning-in-Place (CIP) & Fouling Reversibility Breakdown...")
    t0 = time.time()
    cip_res = simulate_cip_reversibility(output_png=data_dir / "cip_fouling_reversibility.png")
    print(f"        -> Output: data/cip_fouling_reversibility.png ({time.time() - t0:.2f}s)")
    print(f"           Reversible Fouling: {cip_res['reversible_fouling_pct']:.1f}% | Irreversible: {cip_res['irreversible_fouling_pct']:.1f}% | 5-Cycle FR: {cip_res['cycle_5_flux_recovery_pct']:.1f}%")

    # Stage 15: Techno-Economics, Membrane Sizing & Engineering Summary
    print("\n[15/15] Industrial Techno-Economics, Membrane Sizing & Engineering Report...")
    t0 = time.time()
    econ_res = calculate_techno_economics(batch_volume_L=100.0, output_png=data_dir / "techno_economics_summary.png")
    generate_report(data_dir / "MULTISCALE_SIMULATION_REPORT.md")
    print(f"        -> Output: data/techno_economics_summary.png ({time.time() - t0:.2f}s)")
    print(f"        -> Output: data/MULTISCALE_SIMULATION_REPORT.md")
    print(f"           100L Batch Req Area: {econ_res['required_membrane_area_m2']:.2f} m^2 ({econ_res['installed_cassettes']} Cassettes)")

    total_elapsed = time.time() - start_total
    print("\n" + "=" * 78)
    print(f"  COMPLETE MULTISCALE SUITE (15 MODULES) EXECUTED SUCCESSFULLY in {total_elapsed:.2f}s")
    print("=" * 78)


if __name__ == "__main__":
    main()
