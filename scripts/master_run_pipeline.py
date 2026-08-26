"""Master end-to-end pipeline runner for the multiscale bioparticle transport project.

Executes all 14 unified modules:
1. Microscale CG-MD Coacervation simulation & 4-panel structure analysis.
2. Microscale 3D animated GIF generation.
3. Hofmeister Salt Series & Ionic Screening Phase Shift.
4. Droplet Coalescence & Ostwald Ripening Kinetics.
5. Flory-Huggins LLPS Phase Coexistence & Binodal Diagram.
6. Continuum transport & boundary layer dynamics simulation.
7. Continuum filtration animated time-series GIF generation.
8. Parametric Limiting Flux vs. TMP curve generation across shear rates.
9. Rigorous 2D Continuum Process Optimization & Design Space Map.
10. Hermia's 4-Mechanism Fouling Diagnostic.
11. Full Fed-Batch Concentration (10x VCF) + Diafiltration Sequence.
12. Constant-Volume Diafiltration (CVD) & Buffer Exchange kinetics.
13. Cleaning-in-Place (CIP) & Fouling Reversibility Breakdown.
14. Techno-Economics, Membrane Sizing & Engineering Summary Report.
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
    print("\n[1/14] Microscale CG-MD Structural & Radial Density Profile...")
    t0 = time.time()
    run_md_snapshots()
    print(f"       -> Output: data/gromacs_md_visualization.png ({time.time() - t0:.2f}s)")

    # Stage 2: Microscale 3D Rotating Time-Series GIF
    print("\n[2/14] Microscale 3D Coacervation Animated Trajectory GIF...")
    t0 = time.time()
    run_md_animation()
    print(f"       -> Output: data/gromacs_md_coacervation.gif ({time.time() - t0:.2f}s)")

    # Stage 3: Hofmeister Salt Screening & T_t Phase Shift
    print("\n[3/14] Hofmeister Salt Series & Ionic Screening Thermodynamics...")
    t0 = time.time()
    simulate_hofmeister_phase_behavior(output_png=data_dir / "hofmeister_salt_screening.png")
    print(f"       -> Output: data/hofmeister_salt_screening.png ({time.time() - t0:.2f}s)")

    # Stage 4: Droplet Coalescence & Ripening Kinetics
    print("\n[4/14] Condensate Droplet Coalescence & Ostwald Ripening Kinetics...")
    t0 = time.time()
    simulate_droplet_coalescence(output_png=data_dir / "droplet_coalescence_kinetics.png")
    print(f"       -> Output: data/droplet_coalescence_kinetics.png ({time.time() - t0:.2f}s)")

    # Stage 5: Flory-Huggins LLPS Phase Coexistence Diagram
    print("\n[5/14] Flory-Huggins LLPS Binodal & Spinodal Phase Coexistence...")
    t0 = time.time()
    flory_huggins_phase_diagram(output_png=data_dir / "coacervation_phase_diagram.png")
    print(f"       -> Output: data/coacervation_phase_diagram.png ({time.time() - t0:.2f}s)")

    # Stage 6: Continuum TFF Filtration Dynamics Dashboard
    print("\n[6/14] Continuum Tangential Flow Filtration (TFF) Dynamics...")
    t0 = time.time()
    results = run_continuum_simulation(
        params_json=params_json,
        tmp_pa=200_000.0,
        total_time_s=3600.0,
        shear_rate_s_inv=4000.0,
        bulk_conc_g_l=5.0,
        n_steps=100,
    )
    plot_multiscale_results(results, data_dir / "filtration_summary_figure.png")
    print(f"       -> Output: data/filtration_summary_figure.png ({time.time() - t0:.2f}s)")

    # Stage 7: Continuum Dynamic Time-Series Animation
    print("\n[7/14] Continuum Filtration Dynamics Time-Series GIF...")
    t0 = time.time()
    run_filtration_animation()
    print(f"       -> Output: data/multiscale_filtration_timeseries.gif ({time.time() - t0:.2f}s)")

    # Stage 8: Parametric Limiting Flux vs. TMP Curves
    print("\n[8/14] Parametric Limiting Flux Curves across Shear Rates...")
    t0 = time.time()
    model = MdBridgeModel.load_json(params_json)
    sim = ProcessSimulator(model)
    plot_limiting_flux_curves(sim, data_dir / "limiting_flux_curves.png")
    print(f"       -> Output: data/limiting_flux_curves.png ({time.time() - t0:.2f}s)")

    # Stage 9: Rigorous 2D Continuum Process Optimization Map
    print("\n[9/14] 2D Process Optimization Map with SEC Pareto Contours...")
    t0 = time.time()
    run_intense_parameter_sweep(
        params_json=str(params_json),
        output_png=str(data_dir / "intense_optimization_landscape.png"),
        sample_density=60,
    )
    print(f"       -> Output: data/intense_optimization_landscape.png ({time.time() - t0:.2f}s)")

    # Stage 10: Hermia's 4-Mechanism Fouling Diagnostic
    print("\n[10/14] Hermia's 4-Mechanism Fouling Diagnostics...")
    t0 = time.time()
    hermia_scores = analyze_hermia_fouling(results, data_dir / "hermia_fouling_analysis.png")
    print(f"        -> Output: data/hermia_fouling_analysis.png ({time.time() - t0:.2f}s)")

    # Stage 11: Fed-Batch Up-Concentration + Diafiltration Process Sequence
    print("\n[11/14] Fed-Batch Concentration (10x VCF) + Diafiltration Sequence...")
    t0 = time.time()
    fed_res = simulate_fed_batch_purification(output_png=data_dir / "fed_batch_concentration_diafiltration.png")
    print(f"        -> Output: data/fed_batch_concentration_diafiltration.png ({time.time() - t0:.2f}s)")
    print(f"           Concentration: {fed_res['concentration_duration_min']:.1f} min | Diafiltration: {fed_res['diafiltration_duration_min']:.1f} min")

    # Stage 12: Constant-Volume Diafiltration (CVD)
    print("\n[12/14] Constant-Volume Diafiltration (CVD) Buffer Exchange Kinetics...")
    t0 = time.time()
    dia_res = simulate_diafiltration(output_png=data_dir / "diafiltration_summary.png")
    print(f"        -> Output: data/diafiltration_summary.png ({time.time() - t0:.2f}s)")

    # Stage 13: Cleaning-in-Place (CIP) & Fouling Reversibility
    print("\n[13/14] Cleaning-in-Place (CIP) & Fouling Reversibility Breakdown...")
    t0 = time.time()
    cip_res = simulate_cip_reversibility(output_png=data_dir / "cip_fouling_reversibility.png")
    print(f"        -> Output: data/cip_fouling_reversibility.png ({time.time() - t0:.2f}s)")
    print(f"           Reversible Fouling: {cip_res['reversible_fouling_pct']:.1f}% | Irreversible: {cip_res['irreversible_fouling_pct']:.1f}% | 5-Cycle FR: {cip_res['cycle_5_flux_recovery_pct']:.1f}%")

    # Stage 14: Techno-Economics, Sizing & Master Report Generation
    print("\n[14/14] Industrial Techno-Economics, Membrane Sizing & Engineering Report...")
    t0 = time.time()
    tea_res = calculate_techno_economics(output_png=data_dir / "techno_economics_summary.png")
    generate_report(data_dir / "MULTISCALE_SIMULATION_REPORT.md")
    print(f"        -> Output: data/techno_economics_summary.png ({time.time() - t0:.2f}s)")
    print(f"        -> Output: data/MULTISCALE_SIMULATION_REPORT.md")

    print("\n" + "=" * 78)
    print(f"  COMPLETE MULTISCALE SUITE (14 MODULES) EXECUTED SUCCESSFULLY in {time.time() - start_total:.2f}s")
    print("=" * 78)


if __name__ == "__main__":
    main()
