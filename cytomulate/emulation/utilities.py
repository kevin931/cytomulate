# Math computation
import numpy as np
from numpy import random as rd
from collections import Counter

# Polynomials and spline functions
from numpy.polynomial import polynomial
from scipy.interpolate import Akima1DInterpolator
from scipy.interpolate import UnivariateSpline
from scipy.linalg import eigh


def spline_function(x, y, smoothing_factor = 0.5):
    """
    Generate a smoothing spling function
    Parameters
    ----------
    x: an array of x values
    y: an array of y values
    smoothing_factor: the smoothing factor used in spline fitting

    Returns
    -------
    A spline function that can be evaulated at point t
    """
    # We first normalize the x values to the unit interval [0,1]
    x = (x - np.min(x))/(np.max(x) - np.min(x))
    # Then we fit a smoothing spline
    spl = UnivariateSpline(x, y)
    spl.set_smoothing_factor(smoothing_factor)

    def spline_values(t):
        # We subtract the function value at 0 to ensure that the resulting spline always starts at 0
        return spl(t) - spl(0)
    return spline_values


def polynomial_function(coefficients, end_value):
    """
    Generate a polynomial on [0,1]
    Parameters
    ----------
    coefficients: an array of the polynomial coefficients
    end_value: the desired end value

    Returns
    -------
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
        return polynomial.polyval(t, adjusted_polynomial.coef)

    return polynomial_values


def brownian_bridge_function(end_value, N = 5, lb = 0, ub = 1):
    """
    Generate a random function that starts at 0 and ends at end_value
    Parameters
    ----------
    end_value: the desired end value
    N: number of steps when simulating a brownian bridge
    lb: the lower bound of the variance scaling factor
    ub: the upper bound of the variance scaling factor

    Returns
    -------
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
        return spl(t)
    return spline_values


def trajectories(end_values=None, coefficients=None, x=None, y=None, **kwargs):
    """
    Vectorize the spline function or the polynomial function or the brownian bridge function
    Parameters
    ----------
    end_values: an array of the end values
    coefficients: if polynomial is desired, an array of the polynomial coefficients
    x: if spline is sought after, the x values
    y: if spline is sought after, the y values
    kwargs: other arguments for target functions

    Returns
    -------
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

def find_psm(matrix):
    eigen_result = eigh(matrix)
    eigen_vals = eigen_result[0]
    eigen_vals = np.clip(eigen_vals, a_min=0, a_max=None)
    eigen_vals = np.diag(eigen_vals)

    eigen_vecs = eigen_result[1]

    return eigen_vecs @ eigen_vals @ eigen_vecs.T

def univariate_noise_model(type="normal", **kwargs):
    if type == "normal":
        def model(size):
            return rd.normal(**kwargs, size=size)
    elif type == "uniform":
        def model(size):
            return rd.uniform(**kwargs, size=size)
    else:
        raise ValueError('Unknown type')
    return model

def KLdivergence(x, y):
  """Compute the Kullback-Leibler divergence between two multivariate samples.
  Parameters
  ----------
  x : 2D array (n,d)
    Samples from distribution P, which typically represents the true
    distribution.
  y : 2D array (m,d)
    Samples from distribution Q, which typically represents the approximate
    distribution.
  Returns
  -------
  out : float
    The estimated Kullback-Leibler divergence D(P||Q).
  References
  ----------
  Pérez-Cruz, F. Kullback-Leibler divergence estimation of
continuous distributions IEEE International Symposium on Information
Theory, 2008.
  """
  from scipy.spatial import cKDTree as KDTree

  # Check the dimensions are consistent
  x = np.atleast_2d(x)
  y = np.atleast_2d(y)

  n,d = x.shape
  m,dy = y.shape

  assert(d == dy)


  # Build a KD tree representation of the samples and find the nearest neighbour
  # of each point in x.
  xtree = KDTree(x)
  ytree = KDTree(y)

  # Get the first two nearest neighbours for x, since the closest one is the
  # sample itself.
  r = xtree.query(x, k=2, eps=.01, p=2)[0][:,1]
  s = ytree.query(x, k=1, eps=.01, p=2)[0]

  # There is a mistake in the paper. In Eq. 14, the right side misses a negative sign
  # on the first term of the right hand side.
  return -np.log(r/s).sum() * d / n + np.log(m / (n - 1.))


def cell_type_discrepancy(observed_matrix,
                          simulated_matrix,
                          observed_cell_types,
                          simulated_cell_types,
                          cell_type,
                          mean_ord = 2,
                          cov_ord = 2):
    """
    Calculate the discrepancy between the simulated data and the observed data given a cell type
    Parameters
    ----------
    observed_matrix: the observed expression matrix
    simulated_matrix: the simulated expression matrix
    observed_cell_types: an array of the observed cell type labels
    simulated_cell_types: an array of the simulated cell type labels
    cell_type: the desired cell type
    mean_ord: the type of norm for mean
    cov_ord: the type of norm for covariance

    Returns
    -------
    The discrepancies based on mean, covariance, and KL divergence
    """
    # We first extract the corresponding portions
    observed_index = np.where(observed_cell_types == cell_type)[0]
    simulated_index = np.where(simulated_cell_types == cell_type)[0]
    observed_y = observed_matrix[observed_index, :]
    simulated_y = simulated_matrix[simulated_index, :]

    # To find the discrepancy based on mean
    m_obs = np.mean(observed_y, axis=0)
    m_simu = np.mean(simulated_y, axis=0)
    mean_discrepancy = np.linalg.norm(m_obs - m_simu, ord=mean_ord)

    # To find the discrepancy based on covariance
    cov_obs = np.cov(observed_y, rowvar=False)
    cov_simu = np.cov(simulated_y, rowvar=False)
    cov_discrepancy = np.linalg.norm(cov_obs - cov_simu, ord=cov_ord)

    # Finally, the KL divergence
    kl_divergence = KLdivergence(observed_y, simulated_y)

    return mean_discrepancy, cov_discrepancy, kl_divergence


def combined_cell_type_discrepancy(observed_matrix,
                          simulated_matrix,
                          observed_cell_types,
                          simulated_cell_types,
                          mean_ord = 2,
                          cov_ord = 2):
    """
    Calculated a weighted sum of 3 discrepancy measures
    Parameters
    ----------
    observed_matrix: the observed expression matrix
    simulated_matrix: the simulated expression matrix
    observed_cell_types: an array of the observed cell type labels
    simulated_cell_types: an array of the simulated cell type labels
    mean_ord: the type of norm for mean
    cov_ord: the type of norm for covariance

    Returns
    -------
    The discrepancies based on mean, covariance, and KL divergence
    """
    # We will weight the discrepancies by the cell abundances
    cell_types = np.unique(simulated_cell_types)
    n = simulated_cell_types.shape[0]
    n_cell_types = dict(Counter(simulated_cell_types))
    mean_d = 0
    cov_d = 0
    KL_d = 0

    for cell_type in cell_types:
        m, c, k = cell_type_discrepancy(observed_matrix,
                          simulated_matrix,
                          observed_cell_types,
                          simulated_cell_types,
                          cell_type,
                          mean_ord,
                          cov_ord)
        mean_d += n_cell_types[cell_type]/n * m
        cov_d += n_cell_types[cell_type] / n * c
        KL_d += n_cell_types[cell_type] / n * k

    return mean_d, cov_d, KL_d