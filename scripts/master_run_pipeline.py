"""Master end-to-end pipeline runner for the multiscale bioparticle transport project.

Executes all 10 unified modules:
1. Microscale CG-MD Coacervation simulation & 4-panel structure analysis.
2. Microscale 3D animated GIF generation.
3. Flory-Huggins LLPS Phase Coexistence & Binodal Diagram.
4. Continuum transport & boundary layer dynamics simulation.
5. Continuum filtration animated time-series GIF generation.
6. Parametric Limiting Flux vs. TMP curve generation across shear rates.
7. High-intensity 2D Honeycomb optimization landscape.
8. Hermia's 4-Mechanism Fouling Diagnostic.
9. Constant-Volume Diafiltration (CVD) & Buffer Exchange kinetics.
10. Techno-Economics, Membrane Sizing & Engineering Summary Report.
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
from phase_diagram import flory_huggins_phase_diagram
from biotransport.bridge import MdBridgeModel, ProcessSimulator, run_continuum_simulation
from biotransport.visualize import plot_multiscale_results, plot_limiting_flux_curves
from biotransport.animate_filtration import main as run_filtration_animation
from biotransport.parameter_sweep import run_intense_parameter_sweep
from biotransport.hermia_fouling import analyze_hermia_fouling
from biotransport.diafiltration import simulate_diafiltration
from biotransport.techno_economics import calculate_techno_economics
from generate_multiscale_report import generate_report


def main() -> None:
    print("=" * 76)
    print("        MULTISCALE BIOPARTICLE TRANSPORT MASTER PIPELINE")
    print("   Microscale CG-MD Coacervation + Continuum TFF Purification Platform")
    print("=" * 76)
    start_total = time.time()

    data_dir = project_root / "data"
    data_dir.mkdir(parents=True, exist_ok=True)
    params_json = data_dir / "sample_md_params.json"

    # Stage 1: Microscale CG-MD Structural & Density Analysis
    print("\n[1/10] Microscale CG-MD Structural & Radial Density Profile...")
    t0 = time.time()
    run_md_snapshots()
    print(f"       -> Output: data/gromacs_md_visualization.png ({time.time() - t0:.2f}s)")

    # Stage 2: Microscale 3D Rotating Time-Series GIF
    print("\n[2/10] Microscale 3D Coacervation Animated Trajectory GIF...")
    t0 = time.time()
    run_md_animation()
    print(f"       -> Output: data/gromacs_md_coacervation.gif ({time.time() - t0:.2f}s)")

    # Stage 3: Flory-Huggins LLPS Phase Coexistence Diagram
    print("\n[3/10] Flory-Huggins LLPS Binodal & Spinodal Phase Coexistence...")
    t0 = time.time()
    flory_huggins_phase_diagram(output_png=data_dir / "coacervation_phase_diagram.png")
    print(f"       -> Output: data/coacervation_phase_diagram.png ({time.time() - t0:.2f}s)")

    # Stage 4: Continuum TFF Filtration Dynamics Dashboard
    print("\n[4/10] Continuum Tangential Flow Filtration (TFF) Dynamics...")
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
    print(f"          Initial Flux: {results['initial_flux_lmh']:.2f} LMH | Final Flux: {results['final_flux_lmh']:.2f} LMH")
    print(f"          Flux Decline: {results['flux_decline_percent']:.1f}% | Max Wall Conc: {results['max_wall_conc_g_l']:.1f} g/L")

    # Stage 5: Continuum Dynamic Time-Series Animation
    print("\n[5/10] Continuum Filtration Dynamics Time-Series GIF...")
    t0 = time.time()
    run_filtration_animation()
    print(f"       -> Output: data/multiscale_filtration_timeseries.gif ({time.time() - t0:.2f}s)")

    # Stage 6: Parametric Limiting Flux vs. TMP Curves
    print("\n[6/10] Parametric Limiting Flux Curves across Shear Rates...")
    t0 = time.time()
    model = MdBridgeModel.load_json(params_json)
    sim = ProcessSimulator(model)
    plot_limiting_flux_curves(sim, data_dir / "limiting_flux_curves.png")
    print(f"       -> Output: data/limiting_flux_curves.png ({time.time() - t0:.2f}s)")

    # Stage 7: Rigorous 2D Continuum Process Optimization & Design Space Map
    print("\n[7/10] Generating 2D Process Optimization Map with SEC Pareto Contours...")
    t0 = time.time()
    run_intense_parameter_sweep(
        params_json=str(params_json),
        output_png=str(data_dir / "intense_optimization_landscape.png"),
        sample_density=60,
    )
    print(f"       -> Output: data/intense_optimization_landscape.png ({time.time() - t0:.2f}s)")

    # Stage 8: Hermia's 4-Mechanism Fouling Diagnostic
    print("\n[8/10] Hermia's 4-Mechanism Fouling Diagnostics...")
    t0 = time.time()
    hermia_scores = analyze_hermia_fouling(results, data_dir / "hermia_fouling_analysis.png")
    print(f"       -> Output: data/hermia_fouling_analysis.png ({time.time() - t0:.2f}s)")
    for mech, score in hermia_scores.items():
        print(f"          {mech}: R^2 = {score:.4f}")

    # Stage 9: Constant-Volume Diafiltration & Buffer Exchange
    print("\n[9/10] Constant-Volume Diafiltration (CVD) Buffer Exchange Kinetics...")
    t0 = time.time()
    dia_res = simulate_diafiltration(output_png=data_dir / "diafiltration_summary.png")
    print(f"       -> Output: data/diafiltration_summary.png ({time.time() - t0:.2f}s)")
    print(f"          Final Yield: {dia_res['final_product_yield_pct']:.2f}% | Impurity Remaining: {dia_res['final_impurity_remaining_pct']:.4f}%")

    # Stage 10: Techno-Economics, Sizing & Master Report Generation
    print("\n[10/10] Industrial Techno-Economics, Membrane Sizing & Engineering Report...")
    t0 = time.time()
    tea_res = calculate_techno_economics(output_png=data_dir / "techno_economics_summary.png")
    generate_report(data_dir / "MULTISCALE_SIMULATION_REPORT.md")
    print(f"        -> Output: data/techno_economics_summary.png ({time.time() - t0:.2f}s)")
    print(f"        -> Output: data/MULTISCALE_SIMULATION_REPORT.md")
    print(f"           100L Batch Req Area: {tea_res['required_membrane_area_m2']:.2f} m^2 ({tea_res['installed_cassettes']} Cassettes)")

    print("\n" + "=" * 76)
    print(f"  COMPLETE MULTISCALE SUITE (10 MODULES) EXECUTED SUCCESSFULLY in {time.time() - start_total:.2f}s")
    print("=" * 76)


if __name__ == "__main__":
    main()
