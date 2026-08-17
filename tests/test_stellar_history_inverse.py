import numpy as np
from sclh.stellar_history_inverse import rotation_period_days, infer_age_for_rotation_period, integrate_history_to_age, pressure_support_window


def test_rotation_inverse_linear_bracket():
    age=np.array([1000.,2000.,3000.])
    # choose omega from desired periods 30, 50, 70 d
    om=2*np.pi/(np.array([30.,50.,70.])*86400.)/2.67e-6
    m=infer_age_for_rotation_period(age,om,60.,min_age_myr=1000)
    assert m.within_tolerance
    assert abs(m.age_myr-2500.) < 1e-9
    assert abs(m.period_days-60.) < 1e-9


def test_history_integral_constant():
    age=np.array([0.,10.,20.,30.])
    y=np.ones_like(age)*2.0
    got=integrate_history_to_age(age,y,30.,start_age_myr=10.)
    expect=2.0*20.0*1e6*365.25*86400.
    assert abs(got/expect-1) < 1e-12


def test_pressure_support_fraction():
    age=np.array([0.,100.,200.,300.])
    p=np.array([300.,200.,100.,0.])
    a,b,f=pressure_support_window(age,p,300.,50.,150.,start_age_myr=0.1)
    assert a is not None and b is not None
    assert 145 < a < 155
    assert 245 < b < 255
    assert 0.30 < f < 0.36
