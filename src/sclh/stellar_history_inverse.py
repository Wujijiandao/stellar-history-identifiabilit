"""Inverse conditioning of stellar-evolution tracks on present-day observables.

The functions here deliberately separate two operations:
1. conditioning a model track on a present observable (e.g. rotation period), and
2. integrating history-dependent exposures along the compatible track.

No probability over tracks is implied unless an external prior is supplied.
"""
from __future__ import annotations
from dataclasses import dataclass, asdict
import numpy as np

OMEGA_SUN_RAD_S = 2.67e-6
SEC_PER_DAY = 86400.0
SEC_PER_MYR = 1e6 * 365.25 * SEC_PER_DAY


@dataclass(frozen=True)
class RotationMatch:
    age_myr: float
    period_days: float
    omega_sun: float
    within_tolerance: bool

    def to_dict(self) -> dict:
        return asdict(self)


def rotation_period_days(omega_surface_sun: np.ndarray | float) -> np.ndarray:
    omega = np.asarray(omega_surface_sun, dtype=float) * OMEGA_SUN_RAD_S
    if np.any(omega <= 0):
        raise ValueError("rotation rate must be positive")
    return 2.0 * np.pi / omega / SEC_PER_DAY


def infer_age_for_rotation_period(
    age_myr: np.ndarray,
    omega_surface_sun: np.ndarray,
    target_period_days: float,
    *,
    min_age_myr: float = 0.0,
    tolerance_days: float = 1.0,
) -> RotationMatch:
    """Infer the first late-time age at which a track matches a target period.

    Uses linear interpolation in period versus age across a sign-changing bracket.
    If no bracket exists, returns the nearest point only when within tolerance.
    This is a model-conditioned inverse, not a stellar-age posterior.
    """
    age = np.asarray(age_myr, dtype=float)
    period = rotation_period_days(omega_surface_sun)
    if age.shape != period.shape or age.ndim != 1:
        raise ValueError("age and omega must be 1D arrays with identical shape")
    if len(age) < 2 or np.any(np.diff(age) <= 0):
        raise ValueError("age grid must be strictly increasing")
    eligible = age >= float(min_age_myr)
    if not np.any(eligible):
        raise ValueError("no ages satisfy min_age_myr")
    idx = np.where(eligible)[0]
    for i in idx[:-1]:
        y1 = period[i] - target_period_days
        y2 = period[i + 1] - target_period_days
        if y1 == 0:
            return RotationMatch(float(age[i]), float(period[i]), float(np.asarray(omega_surface_sun)[i]), True)
        if y1 * y2 <= 0:
            frac = (target_period_days - period[i]) / (period[i + 1] - period[i])
            t = age[i] + frac * (age[i + 1] - age[i])
            om = np.interp(t, age, np.asarray(omega_surface_sun, dtype=float))
            return RotationMatch(float(t), float(target_period_days), float(om), True)
    j = idx[np.argmin(np.abs(period[idx] - target_period_days))]
    d = abs(period[j] - target_period_days)
    return RotationMatch(float(age[j]), float(period[j]), float(np.asarray(omega_surface_sun)[j]), bool(d <= tolerance_days))


def integrate_history_to_age(age_myr: np.ndarray, y: np.ndarray, end_age_myr: float, *, start_age_myr: float = 10.0) -> float:
    """Trapezoid integral of y(t) dt, with dt converted from Myr to seconds."""
    age = np.asarray(age_myr, dtype=float)
    val = np.asarray(y, dtype=float)
    if age.shape != val.shape or age.ndim != 1:
        raise ValueError("age and y must be 1D arrays with identical shape")
    if not (age[0] <= start_age_myr < end_age_myr <= age[-1]):
        raise ValueError("integration bounds outside track or reversed")
    mask = (age >= start_age_myr) & (age <= end_age_myr)
    ta = age[mask]
    yy = val[mask]
    if len(ta) == 0 or ta[0] > start_age_myr:
        ta = np.insert(ta, 0, start_age_myr)
        yy = np.insert(yy, 0, np.interp(start_age_myr, age, val))
    if ta[-1] < end_age_myr:
        ta = np.append(ta, end_age_myr)
        yy = np.append(yy, np.interp(end_age_myr, age, val))
    return float(np.trapezoid(yy, ta * SEC_PER_MYR))


def pressure_support_window(
    age_myr: np.ndarray,
    pressure_npa: np.ndarray,
    end_age_myr: float,
    p_low_npa: float,
    p_high_npa: float,
    *,
    start_age_myr: float = 10.0,
) -> tuple[float | None, float | None, float]:
    """Return first entry/exit ages and fraction of history inside a pressure support interval.

    This is purely a support-overlap diagnostic. It does not interpolate an escape rate.
    """
    age = np.asarray(age_myr, dtype=float)
    p = np.asarray(pressure_npa, dtype=float)
    if p_low_npa >= p_high_npa:
        raise ValueError("p_low_npa must be < p_high_npa")
    if end_age_myr <= start_age_myr:
        raise ValueError("end_age_myr must exceed start_age_myr")
    # Dense interpolation is sufficient for a support-duration diagnostic.
    grid = np.linspace(start_age_myr, end_age_myr, max(1000, int(end_age_myr - start_age_myr) + 1))
    pg = np.interp(grid, age, p)
    inside = (pg >= p_low_npa) & (pg <= p_high_npa)
    if not np.any(inside):
        return None, None, 0.0
    g = grid[inside]
    frac = float(np.trapezoid(inside.astype(float), grid) / (end_age_myr - start_age_myr))
    return float(g[0]), float(g[-1]), frac
