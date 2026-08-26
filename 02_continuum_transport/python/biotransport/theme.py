"""Minimalist design theme and palette constants for publication figures."""

import matplotlib.pyplot as plt

# Minimalist Nordic / Modern Academic Palette
PALETTE = {
    "slate_dark": "#1e293b",
    "slate_med": "#475569",
    "slate_light": "#94a3b8",
    "bg_grid": "#e2e8f0",
    "teal": "#0f766e",
    "teal_light": "#14b8a6",
    "blue": "#1d4ed8",
    "indigo": "#4338ca",
    "amber": "#b45309",
    "copper": "#c2410c",
    "sage": "#047857",
    "crimson_muted": "#be123c",
}

def apply_minimalist_theme():
    """Sets clean minimalist rcParams for all Matplotlib plots."""
    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Segoe UI", "Helvetica Neue", "Arial", "DejaVu Sans"],
        "axes.edgecolor": "#94a3b8",
        "axes.linewidth": 0.9,
        "axes.labelcolor": "#1e293b",
        "axes.titlesize": 12,
        "axes.titleweight": "600",
        "axes.labelsize": 11,
        "xtick.color": "#475569",
        "ytick.color": "#475569",
        "xtick.labelsize": 10,
        "ytick.labelsize": 10,
        "grid.color": "#e2e8f0",
        "grid.linestyle": "--",
        "grid.alpha": 0.6,
        "legend.frameon": True,
        "legend.framealpha": 0.92,
        "legend.edgecolor": "#cbd5e1",
        "legend.fontsize": 9.5,
        "figure.titlesize": 14,
        "figure.titleweight": "700",
    })
