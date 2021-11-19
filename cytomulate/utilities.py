# Math computation
import numpy as np
from numpy import random as rd
from scipy.interpolate import Akima1DInterpolator

# List manipulation
from copy import deepcopy

# File IO
import _csv
import csv
import os

from typing import Union, Optional, Any, List


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


def generate_random_tree(node_ids):
    """Generate a random tree given the nodes and a Prufer sequence
    :param node_ids: IDs of the nodes
    :return: a nested list whose elements are pairs of ids in ascending order
    """
    S = generate_prufer_sequence(node_ids)
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


class FileIO():
    
    @staticmethod
    def load_data(file: str,
                  col_names: bool = True,
                  drop_columns: Optional[Union[int, List[int]]]=None,
                  delim: str = "\t",
                  dtype = float
                  ) -> List["np.ndarray"]:
        
        """Load CyTOF data into a list of arrays
        
        :param file: Full file path
        :type file: str
        :param col_names: Whether the first row is column names, defaults to True
        :type col_names: bool
        :param drop_columns: Indicies of columns to drop (starts at 0), defaults to None
        :type drop_columns: bool, optional
        :param delim: File delimiter, defaults to "\t"
        :type delim: str
        :param dtype: Expression matrix data type (not including the col_names if applicable), defaults to float.
        :type stype: float

        :return: A list of two arrays with column names and expression matrix
        :rtype: List[np.ndarray]
        """
    
        return_files: List["np.ndarray"] = []    
        skiprows: int=0
        
        if col_names:
            names: "np.ndarray" = np.loadtxt(fname=file, dtype ="str", max_rows=1, delimiter=delim)
            if drop_columns is not None:
                names = np.delete(names, drop_columns)
            return_files.append(names)
            skiprows = 1
        else:
            return_files.append(np.array(None))

        # Load Data
        f: "np.ndarray" = np.loadtxt(fname=file, dtype=dtype, skiprows=skiprows, delimiter=delim)
        if drop_columns is not None:
            f = np.delete(f, drop_columns, axis=1)
 
        return_files.append(f)
            
        return return_files
    
    
    @staticmethod
    def save_2d_list_to_csv(data: List[List[Any]], path: str):
        """Save a nested list to a CSV file.

        :param data: The nested list to be written to disk
        :type data: List[List[Any]]
        :param path: Path to save the CSV file
        :type path: str
        """
        
        i: int
        j: int   
        
        with open(path, "w") as f:      
            w: "_csv._writer" = csv.writer(f)
            for i in range(len(data[0])):
                row: List[Any] = []
                for j in range(len(data)):
                    row.append(data[j][i])
                w.writerow(row)
            
            
    @staticmethod
    def save_np_array(array: "np.ndarray",
                      path: str,
                      col_names: Optional["np.ndarray"]=None,
                      dtype: str="%.18e") -> None:
        """Save a NumPy array to a plain text file

        :param array: The NumPy array to be saved
        :type array: np.ndarray
        :param file: Path to save the plain text file
        :type file: str
        :param col_names: Column names to be save as the first row, defaults to None
        :type col_names: np.ndarray, optional
        :param dtype: NumPy data type, defaults to "%.18e"
        :type dtype: str, optional
        """
        with open(path, "w") as f:
            if col_names is not None:
                f.write("\t".join(list(map(str, col_names))))
                f.write("\n")
            np.savetxt(f, array, delimiter="\t", fmt=dtype)
            
    
    @staticmethod
    def make_dir(dir_path: str, add_number_if_dir_exists: bool = False, _counter: int=0) -> str:
        """Create a new directory

        :param dir_path: Path to the new directory to be created
        :type dir_path: str
        :param add_number_if_dir_exists: If the directory already exists, append a number to the
            name until creation of directory is successful, defaults to True
        :type add_number_if_dir_exists: bool, optional
        :return: The path to the new directory
        :rtype: str
        
        .. Warning:: The ``add_number_if_dir_exists`` can be dangerously when this method is run inside
            of a loop. This behavior may be removed in the future.
        """
        dir_path = dir_path.rstrip("/")
        if _counter==0:
            new_dir_path = dir_path
        else:
            new_dir_path = dir_path + str(_counter)
        
        try:
            os.makedirs(new_dir_path)
        except FileExistsError:
            if add_number_if_dir_exists: 
                new_dir_path = FileIO.make_dir(dir_path, _counter = _counter+1)
            else:
                raise
            
        return new_dir_path