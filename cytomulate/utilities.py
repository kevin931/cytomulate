# Math computation
import numpy as np
from numpy import random as rd
from scipy.interpolate import Akima1DInterpolator


def linear_function(end_value):
    def line_segment(t):
        return np.array(t) * end_value

    return line_segment


def smooth_brownian_bridge(end_values, N = 5, function_type = "linear", lb = 0, ub = 1):
    end_values = np.array(end_values).reshape((-1,1))
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
