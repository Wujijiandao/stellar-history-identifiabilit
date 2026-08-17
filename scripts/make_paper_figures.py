from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data" / "derived"
OUT = ROOT / "figures"
OUT.mkdir(exist_ok=True)

toi = pd.read_csv(DATA / "toi700_rotation_conditioned_histories.csv")
lhs = pd.read_csv(DATA / "lhs1140_rotation_compatible_histories.csv")
activity = pd.read_csv(DATA / "lhs1140_activity_model_comparison.csv")
sens = pd.read_csv(DATA / "lhs1140_weight_sensitivity.csv")

def save(fig, stem):
    fig.savefig(OUT / f"{stem}.pdf", bbox_inches="tight")
    fig.savefig(OUT / f"{stem}.svg", bbox_inches="tight")
    plt.close(fig)

fig, ax = plt.subplots(figsize=(7.2, 5.2))
for mass, g in toi.groupby("stellar_mass_msun"):
    ax.scatter(g.current_euv_flux_erg_s_cm2, g.cumulative_euv_fluence_erg_cm2/1e19,
               s=28, alpha=0.8, label=f"{mass:.2f} $M_\\odot$ grid")
ax.set_xlabel("Model-predicted current EUV flux (erg s$^{-1}$ cm$^{-2}$)")
ax.set_ylabel("Cumulative EUV fluence since 10 Myr ($10^{19}$ erg cm$^{-2}$)")
ax.set_title("TOI-700: present forcing converges while cumulative exposure remains broad")
ax.legend(frameon=False); ax.grid(alpha=0.2); fig.tight_layout()
save(fig, "figure_1_toi700_current_vs_cumulative_euv")

fig, ax = plt.subplots(figsize=(7.2, 5.2))
ax.scatter(lhs.current_euv_flux_erg_s_cm2, lhs.cumulative_euv_fluence_erg_cm2/1e19, s=32, alpha=0.8)
ax.set_xlabel("Model-predicted current EUV flux (erg s$^{-1}$ cm$^{-2}$)")
ax.set_ylabel("Cumulative EUV fluence since 10 Myr ($10^{19}$ erg cm$^{-2}$)")
ax.set_title("LHS 1140: <1% present-EUV spread does not identify cumulative exposure")
ax.grid(alpha=0.2); fig.tight_layout()
save(fig, "figure_2_lhs1140_current_vs_cumulative_euv")

fig, ax = plt.subplots(figsize=(7.2, 5.2))
ax.scatter(lhs.current_x_luminosity_erg_s/1e26, lhs.cumulative_euv_fluence_erg_cm2/1e19, s=32, alpha=0.8)
ax.set_xlabel("Native model current $L_X$ ($10^{26}$ erg s$^{-1}$)")
ax.set_ylabel("Cumulative EUV fluence since 10 Myr ($10^{19}$ erg cm$^{-2}$)")
ax.set_title("LHS 1140: present activity carries little ranking information for history")
ax.grid(alpha=0.2); fig.tight_layout()
save(fig, "figure_3_lhs1140_current_lx_vs_history")

fig, ax = plt.subplots(figsize=(7.2, 5.2))
x = np.arange(len(activity))
y = activity.lx_erg_s.to_numpy()/1e26
ax.scatter(x, y, s=90)
spin_idx = next(i for i, lab in enumerate(activity.mapping.tolist()) if "Spinelli" in lab)
spin = y[spin_idx]
ax.errorbar([spin_idx], [spin], yerr=[[spin-1.13], [1.53-spin]], fmt="none", capsize=5)
ax.set_xticks(x)
ax.set_xticklabels(["Spinelli observed", "MORS native median", "Wright18 central"], rotation=15, ha="right")
ax.set_ylabel("$L_X$ ($10^{26}$ erg s$^{-1}$)")
ax.set_title("LHS 1140 present-activity mapping audit")
ax.grid(axis="y", alpha=0.2); fig.tight_layout()
save(fig, "figure_4_lhs1140_activity_mapping_audit")

fig, ax = plt.subplots(figsize=(7.2, 5.2))
ax.plot(sens.sigma_dex_diagnostic_only, sens.ess_fraction, marker="o", label="Effective-history fraction")
ax.plot(sens.sigma_dex_diagnostic_only, sens.euv_width_ratio_to_unweighted, marker="s",
        label="EUV 5-95% width / unweighted width")
ax.set_xscale("log")
ax.set_xlabel("Declared activity-discrepancy scale (dex; diagnostic only)")
ax.set_ylabel("Fraction or relative width")
ax.set_title("LHS 1140: activity weighting does not contract the historical interval")
ax.legend(frameon=False); ax.grid(alpha=0.2); fig.tight_layout()
save(fig, "figure_5_lhs1140_information_sensitivity")
