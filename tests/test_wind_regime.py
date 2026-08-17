import numpy as np
from sclh.wind_regime import ram_pressure_from_mdot_npa, chi_history, first_upcrossing_age

def test_ram_pressure_linear_mdot_and_v():
    p=ram_pressure_from_mdot_npa(2e-14,400,1)
    assert 1.0 < p < 2.0
    assert abs(ram_pressure_from_mdot_npa(4e-14,400,1)/p-2)<1e-12
    assert abs(ram_pressure_from_mdot_npa(2e-14,800,1)/p-2)<1e-12

def test_chi_decreases_with_pressure():
    c=chi_history(100,np.array([1.,4.,16.]))
    assert np.all(np.diff(c)<0)

def test_crossing_interpolation():
    age=np.array([0.,10.,20.]); y=np.array([.5,.8,1.2])
    x=first_upcrossing_age(age,y,1)
    assert abs(x-15)<1e-12
