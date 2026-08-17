from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr

from sclh.history_information import (
    normalized_log_likelihood_weights,
    effective_sample_size,
    weighted_quantile,
    log10_spread,
)

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "E5C6"
OUT.mkdir(parents=True, exist_ok=True)

lhs = pd.read_csv(ROOT / "results/E5C3/e5c3_rotation_compatible_histories.csv")
toi = pd.read_csv(ROOT / "results/E5C2/e5c2_mors_rotation_conditioned_histories.csv")

obs_lx = 1.34e26
sigmas = [0.02, 0.05, 0.10, 0.20, 0.36, 0.50, 1.00]

# Baseline empirical distributions, not posterior distributions.
base_q_euv = np.quantile(lhs.cumulative_euv_fluence_erg_cm2, [0.05, 0.5, 0.95])
base_q_wind = np.quantile(lhs.cumulative_wind_column_kg_m2, [0.05, 0.5, 0.95])
base_width_euv = base_q_euv[2] - base_q_euv[0]
base_width_wind = base_q_wind[2] - base_q_wind[0]

rows = []
for s in sigmas:
    w = normalized_log_likelihood_weights(lhs.current_x_luminosity_erg_s, obs_lx, s)
    qe = weighted_quantile(lhs.cumulative_euv_fluence_erg_cm2, w)
    qw = weighted_quantile(lhs.cumulative_wind_column_kg_m2, w)
    rows.append({
        "sigma_dex_diagnostic_only": s,
        "effective_sample_size": effective_sample_size(w),
        "ess_fraction": effective_sample_size(w) / len(lhs),
        "max_normalized_weight": float(w.max()),
        "euv_q05": qe[0], "euv_q50": qe[1], "euv_q95": qe[2],
        "euv_width_ratio_to_unweighted": float((qe[2]-qe[0]) / base_width_euv),
        "wind_q05": qw[0], "wind_q50": qw[1], "wind_q95": qw[2],
        "wind_width_ratio_to_unweighted": float((qw[2]-qw[0]) / base_width_wind),
    })
weights_df = pd.DataFrame(rows)
weights_df.to_csv(OUT / "e5c6_lhs1140_weight_sensitivity.csv", index=False)

# Information/adequacy audit table. 'history_discrimination' refers to the
# ability to rank histories *within the already rotation-compatible set*.
px_euv = pearsonr(lhs.current_x_luminosity_erg_s, lhs.cumulative_euv_fluence_erg_cm2)
sx_euv = spearmanr(lhs.current_x_luminosity_erg_s, lhs.cumulative_euv_fluence_erg_cm2)
px_wind = pearsonr(lhs.current_x_luminosity_erg_s, lhs.cumulative_wind_column_kg_m2)
sx_wind = spearmanr(lhs.current_x_luminosity_erg_s, lhs.cumulative_wind_column_kg_m2)

info_rows = [
    {
        "system": "TOI-700 d", "observable": "present rotation",
        "role": "history admissibility gate", "history_discrimination": "partial",
        "current_state_log10_spread_dex": log10_spread(toi.current_euv_flux_erg_s_cm2),
        "cumulative_euv_log10_spread_dex": log10_spread(toi.cumulative_euv_fluence_erg_cm2),
        "cumulative_wind_log10_spread_dex": log10_spread(toi.cumulative_wind_column_kg_m2),
        "note": "194 rotation-compatible histories retain large cumulative-forcing spread."
    },
    {
        "system": "LHS 1140 b", "observable": "present rotation",
        "role": "history admissibility gate", "history_discrimination": "partial",
        "current_state_log10_spread_dex": log10_spread(lhs.current_euv_flux_erg_s_cm2),
        "cumulative_euv_log10_spread_dex": log10_spread(lhs.cumulative_euv_fluence_erg_cm2),
        "cumulative_wind_log10_spread_dex": log10_spread(lhs.cumulative_wind_column_kg_m2),
        "note": "97 rotation-compatible histories retain cumulative-forcing spread."
    },
    {
        "system": "LHS 1140 b", "observable": "present X-ray luminosity",
        "role": "model-adequacy / observation-mapping audit",
        "history_discrimination": "weak within rotation-compatible set",
        "current_state_log10_spread_dex": log10_spread(lhs.current_x_luminosity_erg_s),
        "cumulative_euv_log10_spread_dex": log10_spread(lhs.cumulative_euv_fluence_erg_cm2),
        "cumulative_wind_log10_spread_dex": log10_spread(lhs.cumulative_wind_column_kg_m2),
        "note": "MORS Lx is nearly common across histories; X-ray bandpass/model discrepancy must be homogenized before a physical rejection gate."
    },
    {
        "system": "LHS 1140 b", "observable": "age > 5 Gyr",
        "role": "support check", "history_discrimination": "none for retained set",
        "current_state_log10_spread_dex": np.nan,
        "cumulative_euv_log10_spread_dex": log10_spread(lhs.cumulative_euv_fluence_erg_cm2),
        "cumulative_wind_log10_spread_dex": log10_spread(lhs.cumulative_wind_column_kg_m2),
        "note": "All retained MORS crossing ages are >10.5 Gyr; the adopted lower bound removes no history."
    },
]
pd.DataFrame(info_rows).to_csv(OUT / "e5c6_observable_information_audit.csv", index=False)

summary = {
    "experiment": "E5-C6 multi-observable history-information and observation-operator audit",
    "central_result": "Present-day observables can audit the forward model without materially identifying the cumulative forcing history when they are conditionally redundant within a current-state-matched history set.",
    "theorem": {
        "statement": "If Y is conditionally independent of latent history H given the matched present state C and a fixed observation model m, then p(H|C,Y,m)=p(H|C,m).",
        "interpretation": "A present-day activity datum that is a memoryless function of the already-conditioned current state is an adequacy check, not additional historical information."
    },
    "bandpass_audit": {
        "mors_xray_definition": "approximately 0.1-2.4 keV (Johnstone et al. 2021)",
        "spinelli_xmm_table_definition": "0.2-2.4 keV; authors convert measured luminosities when comparing in a ROSAT-band rotation-activity diagram",
        "decision": "DOWNGRADE_E5C4_HARD_REJECTION_UNTIL_OBSERVATION_OPERATOR_HOMOGENIZED",
        "affected_claim": "physical model rejection based on hard Lx interval",
        "unaffected_claims": ["TOI-700 current-state/history degeneracy", "LHS 1140 current-state/history degeneracy", "need for explicit activity-model class and discrepancy"]
    },
    "lhs1140": {
        "n_rotation_compatible_histories": int(len(lhs)),
        "current_x_luminosity_max_over_min": float(lhs.current_x_luminosity_erg_s.max()/lhs.current_x_luminosity_erg_s.min()),
        "current_x_luminosity_log10_spread_dex": log10_spread(lhs.current_x_luminosity_erg_s),
        "cumulative_euv_max_over_min": float(lhs.cumulative_euv_fluence_erg_cm2.max()/lhs.cumulative_euv_fluence_erg_cm2.min()),
        "cumulative_wind_max_over_min": float(lhs.cumulative_wind_column_kg_m2.max()/lhs.cumulative_wind_column_kg_m2.min()),
        "pearson_current_lx_vs_cumulative_euv": {"r": float(px_euv.statistic), "p": float(px_euv.pvalue)},
        "spearman_current_lx_vs_cumulative_euv": {"rho": float(sx_euv.statistic), "p": float(sx_euv.pvalue)},
        "pearson_current_lx_vs_cumulative_wind": {"r": float(px_wind.statistic), "p": float(px_wind.pvalue)},
        "spearman_current_lx_vs_cumulative_wind": {"rho": float(sx_wind.statistic), "p": float(sx_wind.pvalue)},
        "age_lower_bound_retained_fraction": float(np.mean(lhs.rotation_match_age_myr >= 5000.0)),
        "diagnostic_weighting": {
            "sigma_grid_dex": sigmas,
            "warning": "These normalized weights are a sensitivity diagnostic, not a calibrated posterior or fitted discrepancy prior.",
            "result": "For sigma >=0.05 dex, ESS remains >95% of the 97-history set and the 5-95% cumulative-EUV interval does not contract; present Lx provides negligible history discrimination inside the rotation-matched set."
        }
    },
    "toi700": {
        "n_rotation_compatible_histories": int(len(toi)),
        "current_euv_max_over_min": float(toi.current_euv_flux_erg_s_cm2.max()/toi.current_euv_flux_erg_s_cm2.min()),
        "cumulative_euv_max_over_min": float(toi.cumulative_euv_fluence_erg_cm2.max()/toi.cumulative_euv_fluence_erg_cm2.min()),
        "cumulative_wind_max_over_min": float(toi.cumulative_wind_column_kg_m2.max()/toi.cumulative_wind_column_kg_m2.min()),
    },
    "status": "HISTORY_PARTIALLY_IDENTIFIED_PRESENT_ACTIVITY_MAINLY_ADEQUACY_INFORMATION",
    "next_model": "p(H_star,m_rot,m_struct,m_act,delta_act | P_rot,Lx,age,Mstar,...), with bandpass-homogenized observation operators and no global historical XUV rescaling from one present-day datum.",
}
with open(OUT / "e5c6_summary.json", "w", encoding="utf-8") as f:
    json.dump(summary, f, ensure_ascii=False, indent=2)

# Figures (separate plots; Matplotlib defaults, no custom style/color).
import matplotlib.pyplot as plt
fig, ax = plt.subplots(figsize=(5.6, 4.2))
ax.scatter(lhs.current_x_luminosity_erg_s/1e26, lhs.cumulative_euv_fluence_erg_cm2/1e19, s=22)
ax.axvspan(1.13, 1.53, alpha=0.15)
ax.set_xlabel(r"Present $L_X$ ($10^{26}$ erg s$^{-1}$)")
ax.set_ylabel(r"Cumulative EUV fluence ($10^{19}$ erg cm$^{-2}$)")
ax.set_title("LHS 1140 b: present activity vs cumulative history")
fig.tight_layout()
fig.savefig(OUT / "e5c6_lhs1140_present_lx_vs_cumulative_euv.png", dpi=220, bbox_inches="tight")
plt.close(fig)

fig, ax = plt.subplots(figsize=(5.6, 4.2))
ax.plot(weights_df["sigma_dex_diagnostic_only"], weights_df["ess_fraction"], marker="o", label="ESS fraction")
ax.plot(weights_df["sigma_dex_diagnostic_only"], weights_df["euv_width_ratio_to_unweighted"], marker="s", label="EUV 90% width ratio")
ax.axhline(1.0, linewidth=1)
ax.set_xscale("log")
ax.set_xlabel(r"Declared diagnostic discrepancy scale $\sigma_\delta$ (dex)")
ax.set_ylabel("Relative diagnostic")
ax.set_title("Present-day weighting does not contract cumulative history")
ax.legend(frameon=False)
fig.tight_layout()
fig.savefig(OUT / "e5c6_weight_sensitivity.png", dpi=220, bbox_inches="tight")
plt.close(fig)

print(json.dumps(summary, indent=2, ensure_ascii=False))
