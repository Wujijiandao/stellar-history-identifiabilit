"""E5-C7: value-of-information audit for independent stellar-age constraints.

This analysis uses only the frozen E5-C2/E5-C3 derived history tables.  It asks
how strongly a hypothetical independent age interval would contract the
model-conditional identified set for cumulative EUV and wind forcing.

The diagnostic is deliberately worst-case over the unknown interval location.
It is not a posterior and does not assert that a given age precision is currently
attainable for either M dwarf.
"""
from pathlib import Path
import json
import os
import sys

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from sclh.history_information import worst_case_age_window_width

DATA = ROOT / "data" / "derived"
OUT = ROOT / "results" / "E5C7"
OUT.mkdir(parents=True, exist_ok=True)

toi = pd.read_csv(DATA / "toi700_rotation_conditioned_histories.csv")
lhs = pd.read_csv(DATA / "lhs1140_rotation_compatible_histories.csv")

slices = [
    ("TOI-700 combined", toi),
    ("TOI-700 0.40 Msun", toi[toi.stellar_mass_msun == 0.40].copy()),
    ("TOI-700 0.45 Msun", toi[toi.stellar_mass_msun == 0.45].copy()),
    ("LHS 1140 0.20 Msun", lhs.copy()),
]

# Dense enough to show the discrete finite-support envelope while preserving
# exact paper-friendly checkpoints at 0, 25, 50, 100, 200, 250, 300, 500,
# 750, 1000 and 1500 Myr.
base = np.arange(0.0, 1600.0 + 1e-9, 10.0)
checkpoints = np.array([0, 25, 50, 100, 200, 250, 300, 500, 750, 1000, 1500], dtype=float)
windows = np.unique(np.concatenate([base, checkpoints]))

rows = []
for label, df in slices:
    for window in windows:
        rows.append({
            "model_slice": label,
            "age_window_width_myr": float(window),
            "n_histories_total": int(len(df)),
            "worst_case_cumulative_euv_width": worst_case_age_window_width(
                df.rotation_match_age_myr,
                df.cumulative_euv_fluence_erg_cm2,
                window,
            ),
            "worst_case_cumulative_wind_width": worst_case_age_window_width(
                df.rotation_match_age_myr,
                df.cumulative_wind_column_kg_m2,
                window,
            ),
        })

curve = pd.DataFrame(rows)
curve.to_csv(OUT / "e5c7_age_information_leverage.csv", index=False)

# Thresholds used in the manuscript.  For a target multiplicative width,
# the first pair of histories whose Q ratio exceeds the target determines the
# age-window scale at which the worst-case envelope must exceed that target.
def first_violating_age_separation(df, column, target):
    age = df.rotation_match_age_myr.to_numpy(dtype=float)
    q = df[column].to_numpy(dtype=float)
    best = None
    for i in range(len(df)):
        for j in range(i + 1, len(df)):
            ratio = max(q[i], q[j]) / min(q[i], q[j])
            if ratio > target:
                sep = abs(age[i] - age[j])
                if best is None or sep < best:
                    best = float(sep)
    return best

threshold_rows = []
for label, df in slices:
    for target in (1.10, 1.25, 1.50):
        threshold_rows.append({
            "model_slice": label,
            "target_max_width": target,
            "first_euv_violating_age_separation_myr": first_violating_age_separation(
                df, "cumulative_euv_fluence_erg_cm2", target
            ),
            "first_wind_violating_age_separation_myr": first_violating_age_separation(
                df, "cumulative_wind_column_kg_m2", target
            ),
        })
thresholds = pd.DataFrame(threshold_rows)
thresholds.to_csv(OUT / "e5c7_age_information_thresholds.csv", index=False)

summary = {
    "interpretation": (
        "Worst-case finite-support contraction under a hypothetical independent age interval; "
        "not a posterior and not a claim of observationally attainable age precision."
    ),
    "checkpoints": {},
}
for label, df in slices:
    summary["checkpoints"][label] = {}
    for window in (100.0, 300.0, 500.0):
        summary["checkpoints"][label][f"{int(window)}_Myr"] = {
            "euv_width": worst_case_age_window_width(
                df.rotation_match_age_myr, df.cumulative_euv_fluence_erg_cm2, window
            ),
            "wind_width": worst_case_age_window_width(
                df.rotation_match_age_myr, df.cumulative_wind_column_kg_m2, window
            ),
        }

(OUT / "e5c7_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
print("E5-C7 age-information leverage audit complete.")
