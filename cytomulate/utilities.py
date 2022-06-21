# Math computation
import numpy as np
from numpy import random as rd

# Polynomials and spline functions
from numpy.polynomial import polynomial
from scipy.interpolate import Akima1DInterpolator
from scipy.interpolate import UnivariateSpline

# Typing
from typing import Union, Optional, List, Tuple, Callable


def spline_function(x: np.ndarray,
                    y: np.ndarray,
                    smoothing_factor: float = 0.5) -> Callable:
    """Generate a smoothing spline function
    This is mainly used for generating temporal effects
    Parameters
    ----------
    x: np.ndarray
        An array of x values
    y: np.ndarray
        An array of y values
    smoothing_factor: float
        The smoothing factor used in spline fitting

    Returns
    -------
    spline_values: Callable
        A spline function that can be evaluated at point t
    """
    # We first normalize the x values to the unit interval [0,1]
    x = (x - np.min(x))/(np.max(x) - np.min(x))
    # Then we fit a smoothing spline
    spl = UnivariateSpline(x, y)
    spl.set_smoothing_factor(smoothing_factor)

    def spline_values(t):
        # We subtract the function value at 0 to ensure that the resulting spline always starts at 0
        return (spl(t) - spl(0))[()]
    return spline_values


def polynomial_function(coefficients: Union[list, np.ndarray],
                        end_value: float) -> Callable:
    """Generate a polynomial on [0,1]
    This can be used to generate temporal effect and generate differentiation path
    
    Parameters
    ----------
    coefficients: list or np.ndarray
        An array of the polynomial coefficients.
        The resulting polynomial will almost surely not have the same coefficients.
        They are used to specify the rough "shape" of the polynomial
    end_value: float
        The desired end value

    Returns
    -------
    polynomial_values: Callable
        A polynomial function that can be evaluated at point t
    """
    # We first use the provided coefficient to generate the base polynomial
    # However, the resulting polynomial does not guarantee that it will start at 0
    # and end at end_value
    base_polynomial = polynomial.Polynomial(coefficients, domain=[0, 1])
    base_polynomial_at_end_points = polynomial.polyval([0, 1], base_polynomial.coef)
    # We use a linear function to adjust the base polynomial
    adjust_polynomial = polynomial.Polynomial([base_polynomial_at_end_points[0],
                                               base_polynomial_at_end_points[1]-end_value - base_polynomial_at_end_points[0]],
                                              domain=[0, 1])
    adjusted_polynomial_coefficients = polynomial.polysub(base_polynomial.coef, adjust_polynomial.coef)
    adjusted_polynomial = polynomial.Polynomial(adjusted_polynomial_coefficients, domain=[0, 1])

    def polynomial_values(t):
        return (polynomial.polyval(t, adjusted_polynomial.coef))[()]

    return polynomial_values


def brownian_bridge_function(end_value: float,
                             N: int = 5,
                             lb: float = 0,
                             ub: float = 0.1) -> Callable:
    """Generate a random function that starts at 0 and ends at end_value
    This can be used to generate temporal effect and generate differentiation path
    
    Parameters
    ----------
    end_value: float
        The desired end value
    N: int
        Number of steps when simulating a brownian bridge
    lb: float
        The lower bound of the variance scaling factor
    ub: float
        The upper bound of the variance scaling factor

    Returns
    -------
    spline_values: Callable
        A function that can be evaluated at time t
    """
    # We first generate the variance and the time steps
    sigma2 = np.abs(end_value) * np.random.uniform(lb, ub, 1).reshape((-1,1))
    t_interval: "np.ndarray" = np.linspace(0, 1, num=N, endpoint=True)
    delta: float = 1 / (N - 1)

    # We then generate a Wiener process
    wiener_process = np.zeros((1, N))

    wiener_process[0, 1:] = rd.normal(0, np.sqrt(sigma2[0,0]), N - 1) * np.sqrt(delta)

    wiener_process = np.cumsum(wiener_process, axis=1)

    # Then we adjust the process to make sure it starts at 0 and ends at end_value
    brownian_bridge = np.zeros((1, N))
    for i in range(N):
        brownian_bridge[:,[i]] = wiener_process[:,[i]] - t_interval[i] * (wiener_process[:,[N - 1]] - end_value)
    # To smooth the brownian bridge, we use the Akima spline to interpolate
    spl = Akima1DInterpolator(t_interval, brownian_bridge[0,:])

    def spline_values(t):
        return (spl(t))[()]
    return spline_values


def trajectories(end_values: Optional[Union[list, np.ndarray]] = None,
                 coefficients: Optional[Union[list, np.ndarray]] = None,
                 x: Optional[np.ndarray] = None,
                 y: Optional[np.ndarray] = None,
                 **kwargs) -> List[Callable]:
    """Vectorize the spline function or the polynomial function or the brownian bridge function
    
    Parameters
    ----------
    end_values: list or np.ndarray
        An array of the end values
    coefficients: list or np.ndarray
        If polynomial is desired, an array of the polynomial coefficients
    x: np.ndarray
        If spline is sought after, the x values
    y: np.ndarray
        If spline is sought after, the y values
    kwargs: 
        Other arguments for the brownian bridge method or the spline function

    Returns
    -------
    trajectories_functions: List[Callable]
        A list of functions
    """
    trajectories_functions = []
    if end_values is not None:
        # If end_values is supplied, then we know it's either polynomial or brownian bridge
        end_values = np.array(end_values).reshape((-1, 1))
        n_markers = end_values.shape[0]
        if coefficients is None:
            # If coefficients is not supplied, then we know it's brownian bridge
            for i in range(n_markers):
                trajectories_functions.append(brownian_bridge_function(end_value=end_values[i,0], **kwargs))
        else:
            # Otherwise, it's polynomials
            for i in range(n_markers):
                trajectories_functions.append(polynomial_function(coefficients=coefficients,
                                                                  end_value=end_values[i,0]))
    elif (x is not None) and (y is not None):
        # For spline, we need both x and y
        trajectories_functions.append(spline_function(x, y, **kwargs))
    else:
        raise ValueError('Unknown type')

    return trajectories_functions


def univariate_noise_model(noise_distribution: str = "normal",
                           **kwargs) -> Callable:
    """Generate a noise distribution
    This is mainly used to generate background noise in the cytof_data object
    
    Parameters
    ----------
    noise_distribution: str
        Either "normal" or "uniform"
    kwargs:
        extra parameters needed for numpy.random.normal or numpy.random.uniform

    Returns
    -------
    model: Callable
        A RV generator that only takes size as its input
    """
    if noise_distribution == "normal":
        def model(size):
            return rd.normal(**kwargs, size=size)
    elif noise_distribution == "uniform":
        def model(size):
            return rd.uniform(**kwargs, size=size)
    else:
        raise ValueError('Unknown noise distribution')
    return model

