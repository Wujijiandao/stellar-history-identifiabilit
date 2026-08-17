from pathlib import Path
import json
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data" / "derived"

def test_toi700_headline_spreads():
    df = pd.read_csv(DATA / "toi700_rotation_conditioned_histories.csv")
    assert len(df) == 194
    assert np.isclose(df.current_euv_flux_erg_s_cm2.max()/df.current_euv_flux_erg_s_cm2.min(), 1.0770298566, rtol=2e-4)
    assert np.isclose(df.cumulative_euv_fluence_erg_cm2.max()/df.cumulative_euv_fluence_erg_cm2.min(), 2.652, rtol=5e-3)
    assert np.isclose(df.cumulative_wind_column_kg_m2.max()/df.cumulative_wind_column_kg_m2.min(), 2.73, rtol=5e-3)

def test_toi700_fixed_mass_control():
    df = pd.read_csv(DATA / "toi700_rotation_conditioned_histories.csv")
    for mass, expected_euv, expected_wind in [
        (0.40, 2.1962715253, 2.2530017419),
        (0.45, 2.1841394909, 2.4116923881),
    ]:
        g = df[df.stellar_mass_msun == mass]
        assert len(g) == 97
        assert np.isclose(g.cumulative_euv_fluence_erg_cm2.max()/g.cumulative_euv_fluence_erg_cm2.min(), expected_euv, rtol=1e-10)
        assert np.isclose(g.cumulative_wind_column_kg_m2.max()/g.cumulative_wind_column_kg_m2.min(), expected_wind, rtol=1e-10)

def test_lhs1140_headline_spreads():
    df = pd.read_csv(DATA / "lhs1140_rotation_compatible_histories.csv")
    assert len(df) == 97
    assert np.isclose(df.current_euv_flux_erg_s_cm2.max()/df.current_euv_flux_erg_s_cm2.min(), 1.0084657, rtol=5e-4)
    assert np.isclose(df.cumulative_euv_fluence_erg_cm2.max()/df.cumulative_euv_fluence_erg_cm2.min(), 1.56747, rtol=5e-4)
    assert np.isclose(df.cumulative_wind_column_kg_m2.max()/df.cumulative_wind_column_kg_m2.min(), 1.54284, rtol=5e-4)

def test_public_release_governance():
    meta = json.loads((DATA / "paper_results.json").read_text(encoding="utf-8"))
    assert meta["xray_observation_operator"]["hard_physical_rejection_claimed"] is False
    assert meta["inference_scope"]["track_ensemble_interpretation"] == "deterministic finite model support, not an IID sample or population posterior"

def test_v102_age_information_release_metadata():
    meta = json.loads((DATA / "paper_results.json").read_text(encoding="utf-8"))
    assert meta["public_release_version"] == "1.0.2"
    age = meta["age_information_leverage"]["age_window_100_myr"]
    assert np.isclose(age["toi700_0.40_cumulative_euv_width"], 1.119137325368272, rtol=1e-12)
    assert np.isclose(age["toi700_0.45_cumulative_euv_width"], 1.1152167217539684, rtol=1e-12)
    assert np.isclose(age["lhs1140_0.20_cumulative_euv_width"], 1.064458302685109, rtol=1e-12)
    assert age["toi700_combined_cumulative_euv_width"] > 1.7
