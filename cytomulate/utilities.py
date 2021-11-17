# Math computation
import numpy as np
from numpy import random as rd
from scipy.interpolate import Akima1DInterpolator

# List manipulation
from copy import deepcopy


def linear_function(start_value, end_value):
    """ Generate a linear function
    
    :param start_value: the starting value
    :param end_value: the ending value
    :return: a function that interpolate the two points
    """
    def line_segment(t):
        return start_value + t * (end_value - start_value)
    return line_segment


def smooth_brownian_bridge(start_value=0, end_value=0, \
                           N=5, sigma2=1):
    """Simulate a spline-smoothed brownian bridge
    
    :param start_value: the starting value of the brownian bridge
    :param end_value: the ending value of the brownian bridge
    :param N: number of steps
    :param sigma2: variance
    :return: a spline function defined on interval [0,1]
    """
    if sigma2 > 0:
        t_interval = np.linspace(0, 1, num=N, endpoint=True)
        delta = 1 / (N - 1)
        # We first generate a Wiener process
        wiener_process = rd.normal(0, np.sqrt(sigma2), N - 1) * np.sqrt(delta)
        wiener_process = np.cumsum(wiener_process)
        wiener_process = np.concatenate(([0], wiener_process))
        # Then we can construct a Brownian bridge
        brownian_bridge = np.array([start_value + wiener_process[i] - \
                                    t_interval[i] * (wiener_process[N - 1] - end_value + start_value) \
                                    for i in range(N)])
        # Akima spline to interpolate
        spline_function = Akima1DInterpolator(t_interval, brownian_bridge)
        return spline_function
    else:
        return linear_function(start_value, end_value)


def generate_prufer_sequence(node_ids):
    """Generate a Prufer sequence
    
    :param node_ids: an list or an array of IDs of nodes
    :return: a Prufer sequence
    """
    S = rd.choice(node_ids, size=len(node_ids) - 2, replace=True)
    return S


def generate_random_tree(node_ids=[], S=[]):
    """Generate a random tree given the nodes and a Prufer sequence
    
    :param node_ids: IDs of the nodes
    :param S: a Prufer sequence
    :return: a nested list whose elements are pairs of ids in ascending order
    """
    nodes = deepcopy(node_ids)
    seq = deepcopy(S)
    nodes.sort()
    edges = []
    while len(seq) > 0:
        # We find the smallest element in nodes that is
        # not in the Prufer sequence
        # Since it's already sorted, we simply iterate thru the nodes
        counter = 0
        while counter < len(nodes):
            if nodes[counter] not in seq:
                break
            counter += 1
        temp = [nodes[counter], seq[0]]
        temp.sort()
        edges.append(temp)
        nodes = np.delete(nodes, counter)
        seq = np.delete(seq, 0)
    edges.append([nodes[0], nodes[1]])
    return edges


