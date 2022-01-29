# Math computation
import numpy as np
from numpy import random as rd
from scipy.interpolate import Akima1DInterpolator


def linear_function(end_values):
    def line_segments(t):
        return t * end_values

    return line_segments


def smooth_brownian_bridge(end_values, N = 5, sigma2 = 1):

    if sigma2 > 0:
        t_interval: "np.ndarray" = np.linspace(0, 1, num=N, endpoint=True)
        delta: float = 1 / (N - 1)
        # We first generate a Wiener process
        wiener_process: "np.ndarray" = rd.normal(0, np.sqrt(sigma2), N - 1) * np.sqrt(delta)
        wiener_process = np.cumsum(wiener_process)
        wiener_process = np.concatenate(([0], wiener_process))
        # Then we can construct a Brownian bridge
        brownian_bridge: "np.ndarray" = np.array([wiener_process[i] - \
                                                  t_interval[i] * (wiener_process[N - 1] - end_values) \
                                                  for i in range(N)])
        # Akima spline to interpolate
        spline_function: Akima1DInterpolator = Akima1DInterpolator(t_interval, brownian_bridge)
        return spline_function
    else:
        return linear_function(end_values)

