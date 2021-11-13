# Package for math compuatation
import numpy as np
from numpy import random as rd
from scipy.interpolate import Akima1DInterpolator


def smoothBrownianBridge(start_value=0, end_value=0, \
                         N=5, sigma2=1):
    """
    Simulate a spline-smoothed brownian bridge
    :param start_value: the starting value of the brownian bridge
    :param end_value: the ending value of the brownian bridge
    :param N: number of steps
    :param sigma2: variance
    :return: a spline function defined on interval [0,1]
    """
    t_interval = np.linspace(0, 1, num=N, endpoint=True)
    delta = 1 / (N - 1)
    # We first generate a Wiener process
    wiener_process = rd.normal(0, np.sqrt(sigma2), N - 1) * np.sqrt(delta)
    wiener_process = np.cumsum(wiener_process)
    wiener_process = np.concatenate(([0], wiener_process))
    # Then we can construct a Brownian bridge
    brownian_bridge = np.array([start_value + wiener_process[i] - \
                                t_interval[i] * (wiener_process[N-1] - end_value + start_value) \
                                for i in range(N)])
    # Akima spline to interpolate
    splineFunction = Akima1DInterpolator(t_interval, brownian_bridge)
    return splineFunction

def generatePruferSequence(node_ids):
    """
    Generate a Prufer sequence
    :param node_ids: an list or an array of IDs of nodes
    :return: a Prufer sequence in ascending order
    """
    S = rd.choice(node_ids, size=len(node_ids)-2, replace=True)
    S.sort()
    return S


