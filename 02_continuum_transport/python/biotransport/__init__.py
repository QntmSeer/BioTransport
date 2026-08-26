"""Biotransport Python package: Multiscale bioparticle transport and fouling modeling."""

from .bridge import MdBridgeModel, ProcessSimulator, run_continuum_simulation
from .visualize import plot_multiscale_results

__version__ = "0.1.0"
__all__ = [
    "MdBridgeModel",
    "ProcessSimulator",
    "run_continuum_simulation",
    "plot_multiscale_results",
]
