"""Minimal, auditable stellar-wind pressure/regime utilities for SCLH E3-C0.

E3-C0 maps MORS mass-loss histories into ram-pressure histories using an explicit
terminal-wind-speed sensitivity bracket.  It is a magnetic *regime* diagnostic,
not a complete Johnstone 1D wind solution and not an atmospheric escape model.
"""
from __future__ import annotations
import numpy as np

M_SUN_KG = 1.98847e30
YEAR_S = 365.25 * 86400.0
AU_M = 1.495978707e11


def ram_pressure_from_mdot_npa(mdot_msun_yr, v_wind_kms: float, a_au: float = 1.0):
    """Spherical steady-wind ram pressure mdot*v/(4*pi*r^2), in nPa."""
    mdot = np.asarray(mdot_msun_yr, dtype=float) * M_SUN_KG / YEAR_S
    v = float(v_wind_kms) * 1e3
    r = float(a_au) * AU_M
    if v <= 0 or r <= 0:
        raise ValueError("wind speed and orbital distance must be positive")
    return mdot * v / (4.0 * np.pi * r * r) * 1e9


def chi_history(b_planet_nT: float, p_ram_npa):
    """Egan-topology regime coordinate B/[50 nT sqrt(P/nPa)]."""
    p = np.asarray(p_ram_npa, dtype=float)
    if np.any(p <= 0):
        raise ValueError("ram pressure must be positive")
    return float(b_planet_nT) / (50.0 * np.sqrt(p))


def first_upcrossing_age(age_myr, y, threshold: float = 1.0):
    """Linear-interpolated first y crossing from below to >= threshold."""
    age = np.asarray(age_myr, dtype=float)
    val = np.asarray(y, dtype=float)
    if len(age) != len(val) or len(age) == 0:
        raise ValueError("age/y lengths must match and be nonempty")
    above = val >= float(threshold)
    if above[0]:
        return float(age[0])
    idx = np.flatnonzero(above)
    if len(idx) == 0:
        return None
    i = int(idx[0])
    x0, x1 = age[i-1], age[i]
    y0, y1 = val[i-1], val[i]
    if y1 == y0:
        return float(x1)
    f = (float(threshold)-y0)/(y1-y0)
    return float(x0 + f*(x1-x0))
