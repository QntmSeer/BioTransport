from __future__ import annotations
import argparse
import sys
from pathlib import Path

# Add parent directory to path if run as script
sys.path.insert(0, str(Path(__file__).parent.parent))

try:
    from biotransport.bridge import run_continuum_simulation
    from biotransport.visualize import plot_multiscale_results
except ImportError:
    from bridge import run_continuum_simulation
    from visualize import plot_multiscale_results


def main() -> None:
    parser = argparse.ArgumentParser(description="Multiscale Bioparticle Transport Simulator")
    parser.add_argument("--params", type=str, required=True, help="Path to transport_params.json from CG-MD")
    parser.add_argument("--tmp", type=float, default=150_000.0, help="Transmembrane pressure (Pa)")
    parser.add_argument("--time", type=float, default=3600.0, help="Total filtration time (s)")
    parser.add_argument("--steps", type=int, default=100, help="Number of simulation time steps")
    parser.add_argument("--shear", type=float, default=4000.0, help="Crossflow shear rate (s^-1)")
    parser.add_argument("--bulk-conc", type=float, default=5.0, help="Feed bulk bioparticle concentration (g/L)")
    parser.add_argument("--plot", type=str, default=None, help="Optional output PNG path for visualization")
    args = parser.parse_args()

    results = run_continuum_simulation(
        args.params,
        tmp_pa=args.tmp,
        total_time_s=args.time,
        shear_rate_s_inv=args.shear,
        bulk_conc_g_l=args.bulk_conc,
        n_steps=args.steps,
    )

    print("================ Simulation Results (Python Bridge) ================")
    print(f"Initial Flux:     {results['initial_flux_lmh']:.2f} LMH")
    print(f"Final Flux:       {results['final_flux_lmh']:.2f} LMH")
    print(f"Flux Decline:     {results['flux_decline_percent']:.1f}%")
    print(f"Total Permeate:   {results['permeate_volume_l_m2'][-1]:.2f} L/m^2")
    print(f"Max Wall Conc:    {max(results['wall_conc_g_l']):.1f} g/L")
    print("=====================================================================")

    if args.plot:
        plot_multiscale_results(results, args.plot)


if __name__ == "__main__":
    main()
