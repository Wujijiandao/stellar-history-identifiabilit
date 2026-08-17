import pytest
from sclh.model_adequacy import compare_positive_intervals


def test_interval_adequacy_overlap():
    r=compare_positive_intervals(1.0,2.0,1.5,3.0)
    assert r.overlaps
    assert r.multiplicative_gap == 1.0


def test_interval_adequacy_disjoint_high_prediction():
    r=compare_positive_intervals(2.0,3.0,0.8,1.0)
    assert not r.overlaps
    assert abs(r.multiplicative_gap-2.0) < 1e-12


def test_interval_adequacy_rejects_nonpositive():
    with pytest.raises(ValueError):
        compare_positive_intervals(0.0,1.0,1.0,2.0)
