import numpy as np
from sclh.history_information import (
    normalized_log_likelihood_weights,
    effective_sample_size,
    weighted_quantile,
    log10_spread,
)


def test_equal_predictions_give_uniform_weights():
    w = normalized_log_likelihood_weights([2.0, 2.0, 2.0], 1.0, 0.2)
    assert np.allclose(w, np.ones(3) / 3)
    assert np.isclose(effective_sample_size(w), 3.0)


def test_weighted_quantile_is_order_invariant():
    q1 = weighted_quantile([1, 2, 3], [1, 1, 1], [0.5])[0]
    q2 = weighted_quantile([3, 1, 2], [1, 1, 1], [0.5])[0]
    assert np.isclose(q1, q2)


def test_log10_spread():
    assert np.isclose(log10_spread([1, 10]), 1.0)
