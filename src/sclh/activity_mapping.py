"""Independent present-day rotation--activity mappings for SCLH adequacy audits.

The functions in this module are deliberately narrow: they provide a transparent
external cross-check of a *present-day* X-ray observation model.  They are not a
replacement for age-dependent MORS XUV tracks and must not be extrapolated into a
Gyr history without an independently validated time-dependent model.
"""
from __future__ import annotations
from dataclasses import dataclass, asdict
import math

L_SUN_ERG_S = 3.828e33
T_SUN_K = 5772.0


@dataclass(frozen=True)
class Wright2018ActivityPrediction:
    mass_msun: float
    radius_rsun: float
    teff_k: float
    rotation_period_days: float
    convective_turnover_days: float
    rossby_number: float
    log10_lbol_erg_s: float
    log10_lx_over_lbol: float
    lx_erg_s: float
    regime: str

    def to_dict(self) -> dict:
        return asdict(self)


def stellar_lbol_erg_s(radius_rsun: float, teff_k: float) -> float:
    if radius_rsun <= 0 or teff_k <= 0:
        raise ValueError("stellar radius and effective temperature must be positive")
    return L_SUN_ERG_S * radius_rsun**2 * (teff_k / T_SUN_K) ** 4


def wright2018_turnover_days(mass_msun: float) -> float:
    """Central empirical mass--turnover-time fit from Wright et al. (2018), Eq. 6.

    Valid in the published calibration range 0.08 < M/Msun < 1.36.
    """
    if not (0.08 < mass_msun < 1.36):
        raise ValueError("Wright2018 turnover relation is calibrated for 0.08 < M/Msun < 1.36")
    log10_tau = 2.33 - 1.50 * mass_msun + 0.31 * mass_msun**2
    return 10.0**log10_tau


def wright2018_fully_convective_lx(
    mass_msun: float,
    radius_rsun: float,
    teff_k: float,
    rotation_period_days: float,
    *,
    beta: float = -2.3,
    rossby_sat: float = 0.14,
    log10_lx_lbol_sat: float = -3.05,
) -> Wright2018ActivityPrediction:
    """Central fully-convective rotation--activity prediction from Wright+2018.

    Uses the paper's central beta, Ro_sat and saturated fractional luminosity,
    plus its empirical mass--convective-turnover relation.  This function emits
    a deterministic central cross-check only.  Wright et al. explicitly fitted
    an extra scatter term; therefore callers must not interpret this central
    value as a zero-scatter likelihood or a hard support interval.
    """
    if rotation_period_days <= 0:
        raise ValueError("rotation period must be positive")
    tau = wright2018_turnover_days(mass_msun)
    ro = rotation_period_days / tau
    if ro <= rossby_sat:
        log_frac = log10_lx_lbol_sat
        regime = "saturated"
    else:
        log_frac = log10_lx_lbol_sat + beta * math.log10(ro / rossby_sat)
        regime = "unsaturated"
    lbol = stellar_lbol_erg_s(radius_rsun, teff_k)
    lx = lbol * 10.0**log_frac
    return Wright2018ActivityPrediction(
        mass_msun=float(mass_msun), radius_rsun=float(radius_rsun), teff_k=float(teff_k),
        rotation_period_days=float(rotation_period_days), convective_turnover_days=float(tau),
        rossby_number=float(ro), log10_lbol_erg_s=float(math.log10(lbol)),
        log10_lx_over_lbol=float(log_frac), lx_erg_s=float(lx), regime=regime,
    )
