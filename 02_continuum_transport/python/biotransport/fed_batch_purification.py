"""Fed-Batch Ultrafiltration Concentration and Diafiltration (Minimalist Styling)."""

from __future__ import annotations
from pathlib import Path
from typing import Any
import numpy as np
import matplotlib.pyplot as plt
from biotransport.theme import apply_minimalist_theme, PALETTE


def simulate_fed_batch_purification(
    v0_L: float = 100.0,
    target_vcf: float = 10.0,
    c_bulk_0_g_L: float = 5.0,
    target_dv: float = 7.0,
    membrane_area_m2: float = 0.65,
    gamma_dot: float = 4000.0,
    output_png: Path | str = "data/fed_batch_concentration_diafiltration.png",
) -> dict[str, Any]:
    apply_minimalist_theme()
    km = 2.06e-6 * (gamma_dot / 4000.0) ** (1.0 / 3.0)
    c_gel = 420.0

    v_target_ret_L = v0_L / target_vcf
    v_ret_current = v0_L
    time_conc_s = [0.0]
    v_ret_profile = [v0_L]
    c_bulk_profile = [c_bulk_0_g_L]
    flux_profile = []

    dt = 30.0
    t_curr = 0.0

    while v_ret_current > v_target_ret_L:
        c_b = c_bulk_0_g_L * (v0_L / v_ret_current)
        j_lim = km * np.log(max(1.05, c_gel / c_b))
        j_lmh = j_lim * 3.6e6 * 0.85
        flux_profile.append(j_lmh)

        q_perm_L_s = (j_lmh * membrane_area_m2) / 3600.0
        dv_perm = q_perm_L_s * dt
        v_ret_current = max(v_target_ret_L, v_ret_current - dv_perm)

        t_curr += dt
        time_conc_s.append(t_curr)
        v_ret_profile.append(v_ret_current)
        c_bulk_profile.append(c_b)

    flux_profile.append(flux_profile[-1])
    t_conc_end = t_curr

    dvs = np.linspace(0.0, target_dv, 50)
    j_dia_lmh = flux_profile[-1]
    q_dia_L_h = j_dia_lmh * membrane_area_m2
    time_dia_hours = dvs * (v_target_ret_L / q_dia_L_h)
    time_dia_s = t_conc_end + time_dia_hours * 3600.0

    total_times_min = np.concatenate([np.array(time_conc_s) / 60.0, time_dia_s / 60.0])
    total_v_ret = np.concatenate([np.array(v_ret_profile), np.full_like(dvs, v_target_ret_L)])
    total_c_bulk = np.concatenate([np.array(c_bulk_profile), np.full_like(dvs, c_bulk_profile[-1])])
    total_flux = np.concatenate([np.array(flux_profile), np.full_like(dvs, j_dia_lmh)])

    c_imp_conc = 2.0 * (v0_L / np.array(v_ret_profile)) ** (1.0 - 0.95)
    c_imp_dia = c_imp_conc[-1] * np.exp(-(1.0 - 0.05) * dvs)
    total_c_imp = np.concatenate([c_imp_conc, c_imp_dia])

    fig, axes = plt.subplots(2, 2, figsize=(13.5, 9.2), dpi=300)
    fig.suptitle("Industrial Fed-Batch Concentration & Diafiltration Dynamics", fontsize=14, color=PALETTE["slate_dark"])

    t_split_min = t_conc_end / 60.0

    # Panel A: Volume & Concentration Factor
    ax1 = axes[0, 0]
    ax1.plot(total_times_min, total_v_ret, color=PALETTE["blue"], lw=2.2, label=r"Retentate Volume $V_{\mathrm{ret}}(t)$")
    ax1.axvline(x=t_split_min, color=PALETTE["slate_light"], linestyle="--", lw=1.5, label="Stage Split")
    ax1.set_xlabel("Process Time (min)")
    ax1.set_ylabel("Retentate Volume (L)")
    ax1.set_title(r"A. Retentate Volume Reduction ($100\,\mathrm{L} \to 10\,\mathrm{L}$)", loc="left", color=PALETTE["slate_dark"])
    ax1.legend(loc="upper right", frameon=True)

    # Panel B: Bulk Concentrations
    ax2 = axes[0, 1]
    ax2.plot(total_times_min, total_c_bulk, color=PALETTE["teal"], lw=2.2, label=r"Bioparticle Product $C_b(t)$")
    ax2.plot(total_times_min, total_c_imp, color=PALETTE["copper"], lw=2.2, label=r"Host Cell Impurity $C_{\mathrm{imp}}(t)$")
    ax2.axvline(x=t_split_min, color=PALETTE["slate_light"], linestyle="--", lw=1.5)
    ax2.set_xlabel("Process Time (min)")
    ax2.set_ylabel(r"Concentration $(\mathrm{g} \cdot \mathrm{L}^{-1})$")
    ax2.set_title(r"B. Bioparticle Up-Concentration ($10\times$) & Clearance", loc="left", color=PALETTE["slate_dark"])
    ax2.legend(loc="center right", frameon=True)

    # Panel C: Permeate Flux J(t)
    ax3 = axes[1, 0]
    ax3.plot(total_times_min, total_flux, color=PALETTE["slate_dark"], lw=2.2, label=r"Permeate Flux $J(t)$")
    ax3.axvline(x=t_split_min, color=PALETTE["slate_light"], linestyle="--", lw=1.5)
    ax3.set_xlabel("Process Time (min)")
    ax3.set_ylabel(r"Permeate Flux $(\mathrm{L} \cdot \mathrm{m}^{-2} \cdot \mathrm{h}^{-1})$")
    ax3.set_title("C. Flux Decline Profile across Process Sequence", loc="left", color=PALETTE["slate_dark"])
    ax3.legend(loc="upper right", frameon=True)

    # Panel D: Process Stage Duration Breakdown
    ax4 = axes[1, 1]
    stage_durations = [t_split_min, (time_dia_s[-1] - t_conc_end) / 60.0]
    bars = ax4.bar(
        ["Stage I: Concentration\n(VCF = 10x)", "Stage II: Diafiltration\n(DV = 7x)"],
        stage_durations,
        color=[PALETTE["teal"], PALETTE["slate_med"]],
        edgecolor=PALETTE["slate_dark"],
        width=0.48,
        lw=0.8,
    )
    for bar, v in zip(bars, stage_durations):
        ax4.text(bar.get_x() + bar.get_width() / 2, v + 8.0, f"{v:.1f} min", ha="center", fontsize=9.5, fontweight="600", color=PALETTE["slate_dark"])
    ax4.set_ylabel("Duration (min)")
    ax4.set_ylim([0, max(stage_durations) * 1.25])
    ax4.set_title("D. Bioprocess Unit Operation Timeline", loc="left", color=PALETTE["slate_dark"])

    plt.tight_layout()
    out_p = Path(output_png)
    out_p.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(out_p, bbox_inches="tight")
    plt.close()
    print(f"Successfully generated minimalist fed-batch purification summary: {out_p}")

    return {
        "concentration_duration_min": float(t_split_min),
        "diafiltration_duration_min": float((time_dia_s[-1] - t_conc_end) / 60.0),
        "total_batch_duration_min": float(total_times_min[-1]),
        "final_product_conc_g_l": float(total_c_bulk[-1]),
        "final_impurity_conc_g_l": float(total_c_imp[-1]),
    }


if __name__ == "__main__":
    simulate_fed_batch_purification()
