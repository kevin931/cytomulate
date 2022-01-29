# Math computation
import numpy as np
from numpy import random as rd
from scipy.interpolate import Akima1DInterpolator


def linear_function(end_value):
    def line_segment(t):
        return np.array(t) * end_value

    return line_segment


def smooth_brownian_bridge(end_values, N = 5, function_type = "linear", lb = 0, ub = 1):
    end_values = end_values.reshape((-1,1))
    n_markers = end_values.shape[0]
    if function_type == "nonlinear":
        sigma2 = np.abs(end_values) * np.random.uniform(lb, ub, n_markers).reshape((-1,1))
        t_interval: "np.ndarray" = np.linspace(0, 1, num=N, endpoint=True)
        delta: float = 1 / (N - 1)

        # We first generate a Wiener process
        wiener_process = np.zeros((n_markers, N))
        for i in range(n_markers):
            wiener_process[i, 1:] = rd.normal(0, np.sqrt(sigma2[i,0]), N - 1) * np.sqrt(delta)

        wiener_process = np.cumsum(wiener_process, axis=1)
        brownian_bridge = np.zeros((n_markers, N))
        for i in range(N):
            brownian_bridge[:,[i]] = wiener_process[:,[i]] - t_interval[i] * (wiener_process[:,[N - 1]] - end_values)

        # Akima spline to interpolate
        spline_functions = []
        for i in range(n_markers):
            spline_functions.append(Akima1DInterpolator(t_interval, brownian_bridge[i,:]))
        return spline_functions
    elif function_type == "linear":
        linear_functions = []
        for i in range(n_markers):
            linear_functions.append(linear_function(end_values[i,0]))
        return linear_functions
    else:
        raise ValueError('Unknown type')

