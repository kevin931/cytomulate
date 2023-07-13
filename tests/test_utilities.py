import pytest
import numpy as np
from cytomulate.utilities import spline_function, \
    polynomial_function, \
    brownian_bridge_function, \
    trajectories, \
    univariate_noise_model, \
    estimate_noise_model


@pytest.mark.parametrize("x, y, smoothing_factor, t, expected", [
    (np.linspace(0, 1, 10), np.zeros(10), 0.5, 0.0, 0.0),
    (np.linspace(0, 1, 10), np.zeros(10), 0.5, 1.0, 0.0),
])
def test_spline_function(x, y, smoothing_factor, t, expected):
    f = spline_function(x , y, smoothing_factor)
    assert f(t) == expected


@pytest.mark.parametrize("coefficients, end_value, t, expected", [
    ([1, 2, 3], 1.0, 0.0, 0.0),
    ([1, 2, 3], 1.0, 1.0, 1.0),
])
def test_polynomial_function(coefficients, end_value, t, expected):
    f = polynomial_function(coefficients, end_value)
    assert np.around(f(t)) == expected


@pytest.mark.parametrize("end_value, N, lb, ub, t, expected", [
    (5.0, 5, 0, 1, 0.0, 0.0),
    (5.0, 5, 0, 1, 1.0, 5.0),
])
def test_brownian_bridge_function(end_value, N, lb, ub, t, expected):
    f = brownian_bridge_function(end_value, N, lb, ub)
    assert np.around(f(t)) == expected


@pytest.mark.parametrize("end_values, coefficients, x, y, t, expected", [
    ([1, 2, 3, 4], None, None, None, 0, [0, 0, 0, 0]),
    ([1, 2, 3, 4], None, None, None, 1, [1, 2, 3, 4]),
    ([1, 2, 3, 4], [1, 2, 3], None, None, 0, [0, 0, 0, 0]),
    ([1, 2, 3, 4], [1, 2, 3], None, None, 1, [1, 2, 3, 4]),
    (None, None, np.linspace(0, 1, 10), np.zeros(10), 0, 4),
    (None, None, np.linspace(0, 1, 10), np.zeros(10), 1, 4),
    (None, None, None, None, 1, 4),
])
def test_trajectories(end_values, coefficients, x, y, t, expected):
    try:
        if end_values is not None:
            f = trajectories(end_values, coefficients, x, y)
            results = [np.around(f[i](t)) for i in range(4)]
            assert results == expected
        else:
            f = trajectories(end_values, coefficients, x, y)
            assert isinstance(f, list)
    except ValueError:
        assert True


@pytest.mark.parametrize("kwargs, size, expected", [
    ({"noise_distribution":"normal", "loc":0, "scale":1}, 5, (5, )),
    ({"noise_distribution":"half_normal", "loc":0, "scale":1}, (5, 3), (5, 3)),
    ({"noise_distribution":"uniform", "low":0, "high":1}, 5, (5, )),
    ({"noise_distribution":"uniform", "low":0, "high":1}, (5, 3), (5, 3)),
    ({"noise_distribution":"gamma"}, (5, 3), (5, 3)),
])
def test_univariate_noise_model(kwargs, size, expected):
    try:
        f = univariate_noise_model(**kwargs)
        assert f(size).shape == expected
    except ValueError:
        assert True
        

def test_univariate_noise_model_warning():
    with pytest.warns(UserWarning) as record:
        f = univariate_noise_model()
        assert str(record[0].message) == "The default `noise_distribution` is now changed from `normal` to `uniform` as of v0.2.0. Please see the release notes for details."
    

@pytest.mark.parametrize("data, noise_distribution, size, expected", [
    (-np.abs(np.random.normal(size=(5,5), loc=0, scale=1)), "half_normal", 5, (5, )),
    (np.random.uniform(size=(5,5), low=-1, high=0), "uniform", 5, (5, )),
])
def test_estimate_noise_model(data, noise_distribution, size, expected):
    f = estimate_noise_model(data=data,
                             noise_distribution=noise_distribution)
    assert f(size).shape == expected
