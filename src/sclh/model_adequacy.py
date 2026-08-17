"""Small, auditable utilities for model-support and present-observable adequacy checks.

These functions do not assign probabilities.  They compare interval-valued model
predictions with independent observational intervals so that a history inversion
cannot be accepted merely because one observable is matched.
"""
from __future__ import annotations
from dataclasses import dataclass, asdict


@dataclass(frozen=True)
class IntervalAdequacy:
    predicted_low: float
    predicted_high: float
    observed_low: float
    observed_high: float
    overlaps: bool
    multiplicative_gap: float

    def to_dict(self) -> dict:
        return asdict(self)


def compare_positive_intervals(
    predicted_low: float,
    predicted_high: float,
    observed_low: float,
    observed_high: float,
) -> IntervalAdequacy:
    """Compare two positive intervals without inventing a likelihood.

    ``multiplicative_gap`` is 1 when intervals overlap.  Otherwise it is the
    factor by which the nearer upper endpoint of the lower interval must be
    multiplied to reach the nearer lower endpoint of the higher interval.
    """
    vals = [predicted_low, predicted_high, observed_low, observed_high]
    if any(v <= 0 for v in vals):
        raise ValueError("interval endpoints must be positive")
    if predicted_low > predicted_high or observed_low > observed_high:
        raise ValueError("low endpoint must not exceed high endpoint")
    overlaps = not (predicted_high < observed_low or observed_high < predicted_low)
    if overlaps:
        gap = 1.0
    elif predicted_high < observed_low:
        gap = observed_low / predicted_high
    else:
        gap = predicted_low / observed_high
    return IntervalAdequacy(
        float(predicted_low), float(predicted_high), float(observed_low), float(observed_high),
        bool(overlaps), float(gap)
    )
