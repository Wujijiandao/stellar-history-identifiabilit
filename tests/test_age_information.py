from pathlib import Path
import numpy as np
import pandas as pd

from sclh.history_information import worst_case_age_window_width

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data" / "derived"


def test_age_window_width_is_monotone_non_decreasing():
    df = pd.read_csv(DATA / "lhs1140_rotation_compatible_histories.csv")
    widths = [
        worst_case_age_window_width(
            df.rotation_match_age_myr,
            df.cumulative_euv_fluence_erg_cm2,
            w,
        )
        for w in [0, 50, 100, 200, 300, 500, 750, 1000, 1500]
    ]
    assert np.all(np.diff(widths) >= -1e-12)


def test_lhs1140_100myr_age_window_contracts_euv_history():
    df = pd.read_csv(DATA / "lhs1140_rotation_compatible_histories.csv")
    width = worst_case_age_window_width(
        df.rotation_match_age_myr,
        df.cumulative_euv_fluence_erg_cm2,
        100.0,
    )
    assert np.isclose(width, 1.064458302685109, rtol=1e-12)
    assert width < 1.07


def test_toi700_fixed_mass_100myr_age_window_contracts_euv_history():
    df = pd.read_csv(DATA / "toi700_rotation_conditioned_histories.csv")
    for mass, expected in [(0.40, 1.119137325368272), (0.45, 1.1152167217539684)]:
        g = df[df.stellar_mass_msun == mass]
        width = worst_case_age_window_width(
            g.rotation_match_age_myr,
            g.cumulative_euv_fluence_erg_cm2,
            100.0,
        )
        assert np.isclose(width, expected, rtol=1e-12)
        assert width < 1.12


def test_pooled_toi700_age_alone_retains_model_slice_ambiguity():
    df = pd.read_csv(DATA / "toi700_rotation_conditioned_histories.csv")
    pooled = worst_case_age_window_width(
        df.rotation_match_age_myr,
        df.cumulative_euv_fluence_erg_cm2,
        100.0,
    )
    fixed_040 = worst_case_age_window_width(
        df[df.stellar_mass_msun == 0.40].rotation_match_age_myr,
        df[df.stellar_mass_msun == 0.40].cumulative_euv_fluence_erg_cm2,
        100.0,
    )
    assert pooled > 1.7
    assert fixed_040 < 1.12
