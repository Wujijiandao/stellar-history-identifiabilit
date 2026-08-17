import math
import pytest
from sclh.activity_mapping import wright2018_turnover_days, wright2018_fully_convective_lx


def test_wright2018_lhs1140_central_crosscheck():
    p=wright2018_fully_convective_lx(0.1844,0.216,3096,131.0)
    assert 110 < p.convective_turnover_days < 122
    assert 1.05 < p.rossby_number < 1.20
    assert p.regime == 'unsaturated'
    assert -5.25 < p.log10_lx_over_lbol < -5.00
    assert 0.9e26 < p.lx_erg_s < 1.3e26


def test_wright2018_turnover_domain_guard():
    with pytest.raises(ValueError):
        wright2018_turnover_days(0.05)
