"""Diagnostics for how present-day observables constrain latent stellar histories.

These utilities deliberately distinguish *relative-likelihood sensitivity* from a
calibrated posterior. They are intended for model-adequacy and information audits
when observational/model-discrepancy scales are not yet identified.
"""
from __future__ import annotations

import numpy as np


def normalized_log_likelihood_weights(model_values, observed_value, sigma_dex):
    """Return normalized Gaussian weights in log10 space.

    Parameters
    ----------
    model_values : array-like
        Positive model predictions.
    observed_value : float
        Positive observed central value.
    sigma_dex : float
        Declared diagnostic scale in dex. This is not inferred here and must not
        be interpreted as a fitted model-discrepancy prior.
    """
    x = np.asarray(model_values, dtype=float)
    if np.any(x <= 0) or observed_value <= 0 or sigma_dex <= 0:
        raise ValueError("values and sigma_dex must be positive")
    residual = np.log10(x / observed_value)
    logw = -0.5 * (residual / float(sigma_dex)) ** 2
    logw -= np.max(logw)
    w = np.exp(logw)
    w /= np.sum(w)
    return w


def effective_sample_size(weights):
    w = np.asarray(weights, dtype=float)
    s = np.sum(w)
    if s <= 0:
        raise ValueError("weights must have positive total")
    w = w / s
    return float(1.0 / np.sum(w * w))


def weighted_quantile(values, weights, quantiles=(0.05, 0.5, 0.95)):
    x = np.asarray(values, dtype=float)
    w = np.asarray(weights, dtype=float)
    q = np.asarray(quantiles, dtype=float)
    if len(x) != len(w):
        raise ValueError("values and weights must have equal length")
    if np.any(q < 0) or np.any(q > 1):
        raise ValueError("quantiles must lie in [0,1]")
    idx = np.argsort(x)
    x = x[idx]
    w = w[idx]
    w = w / np.sum(w)
    cdf = np.cumsum(w)
    return np.interp(q, cdf, x)


def log10_spread(values):
    x = np.asarray(values, dtype=float)
    if np.any(x <= 0):
        raise ValueError("values must be positive")
    return float(np.log10(np.max(x) / np.min(x)))


def multiplicative_width(values):
    """Return max(values) / min(values) for positive finite support."""
    x = np.asarray(values, dtype=float)
    if x.size == 0:
        raise ValueError("values must be non-empty")
    if np.any(~np.isfinite(x)) or np.any(x <= 0):
        raise ValueError("values must be positive and finite")
    return float(np.max(x) / np.min(x))


def worst_case_age_window_width(match_age_myr, values, window_myr):
    """Worst residual multiplicative width after an independent age-window constraint.

    The input histories are already conditioned on the paper's present-state
    constraints.  Suppose an additional, independent age measurement restricts
    the model-conditional matching age to *some* interval of total width
    ``window_myr``.  Because the interval centre is not fixed here, this function
    returns the largest multiplicative width that can remain in any such interval.

    This is a deterministic finite-support sensitivity diagnostic, not a posterior
    interval and not a claim that the assumed age precision is observationally
    attainable.
    """
    age = np.asarray(match_age_myr, dtype=float)
    x = np.asarray(values, dtype=float)
    if age.size == 0 or x.size == 0 or len(age) != len(x):
        raise ValueError("match_age_myr and values must be non-empty and equal length")
    if np.any(~np.isfinite(age)) or np.any(~np.isfinite(x)) or np.any(x <= 0):
        raise ValueError("inputs must be finite and values positive")
    if window_myr < 0:
        raise ValueError("window_myr must be non-negative")

    width = 1.0
    tol = 1e-9
    # It is sufficient to anchor the interval at an observed support age: the
    # extrema of any finite-support interval occur at support points.
    for left in age:
        mask = (age >= left - tol) & (age <= left + float(window_myr) + tol)
        if np.any(mask):
            width = max(width, multiplicative_width(x[mask]))
    return float(width)
