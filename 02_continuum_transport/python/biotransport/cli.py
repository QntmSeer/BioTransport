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
    parser.add_argument("--plot", type=str, default=None, help="Optional output PNG path for visualization")
    args = parser.parse_args()

    results = run_continuum_simulation(args.params, tmp_pa=args.tmp)

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
