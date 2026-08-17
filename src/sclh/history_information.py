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
