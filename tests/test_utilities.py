import cytomulate
import numpy as np

import pytest
import os
import shutil
import csv
import _csv
import math

from scipy.interpolate import Akima1DInterpolator
from typing import List, Any, Union, Callable
from numpy.typing import ArrayLike


@pytest.mark.parametrize("start_value,end_value,t, expected_value",
            [(1, 2, 0.5, 1.5),
            (1.1, 2.2, 0.5, 1.65),
            (1, 2.2, 0, 1),
            (1.1, 2, 1, 2)]
            )
def test_linear_function(start_value: Union[int, float], end_value: Union[int, float], t: Union[int, float], expected_value: Union[int, float]):
    func: Callable[[Union[int, float]], Union[int, float]] = cytomulate.utilities.linear_function(start_value, end_value)
    interpolation: Union[int, float] = func(t)
    assert callable(func)
    assert isinstance(interpolation, int) or isinstance(interpolation, float)
    assert math.isclose(interpolation, expected_value)
    

# TODO: Resolve the issue of t input type
# @pytest.mark.parametrize("t",
#                          [-0.1, 1.1])
# def test_linear_function_exception(t):
#     func: Callable[[Union[int, float]], Union[int, float]] = cytomulate.utilities.linear_function(1, 2)
#     try:
#         func(t)
#     except ValueError:
#         assert True
#     else:
#         assert False


@pytest.mark.parametrize("t,expected_dim",
                         [(0.3, ()),
                          ([0.3, 0.4], (2,)),
                          (np.array([0.3, 0.4]), (2,))])
def test_smooth_brownian_bridge(t: ArrayLike, expected_dim: tuple):
    func: Akima1DInterpolator = cytomulate.utilities.smooth_brownian_bridge(1, 2, sigma2=0.5) # type: ignore
    interpolation: "np.ndarray" = func(t)
    assert callable(func)
    assert isinstance(func, Akima1DInterpolator)
    assert isinstance(interpolation, np.ndarray)
    assert interpolation.shape == expected_dim


def test_smooth_brownian_bridge_linear():
    func: Callable[[Union[int, float]], Union[int, float]] = cytomulate.utilities.smooth_brownian_bridge(1, 2, sigma2=0) #type: ignore
    interpolation: Union[int, float] = func(0.5)
    assert callable(func)
    assert isinstance(interpolation, int) or isinstance(interpolation, float)
    assert math.isclose(interpolation, 1.5)


def test_generate_prufer_sequence():
    nodes: List[int] = [1, 2, 3, 4, 5]
    prufer: "np.ndarray" = cytomulate.utilities.generate_prufer_sequence(nodes)
    assert isinstance(prufer, np.ndarray)
    assert len(prufer) == len(nodes) - 2
    for element in prufer:
        assert element in nodes
        
        
def test_generate_prufer_sequence_exception():
    nodes: List[int] = [1]
    try:
        cytomulate.utilities.generate_prufer_sequence(nodes)
    except ValueError:
        assert True
    else:
        assert False
    
        
def test_generate_random_tree():
    nodes: List[int] = [1, 2, 3, 4, 5]
    tree: List[List[int]] = cytomulate.utilities.generate_random_tree(nodes)
    assert isinstance(tree, list)
    assert isinstance(tree[0], list)
    assert len(tree[0]) == 2
     

class TestFileIO():
    
    @classmethod
    def setup_class(cls):
        os.mkdir("./tmp_pytest/")
        
        with open("./tmp_pytest/file_read_test_csv.txt", "w") as f:
            tsv_writer: "_csv._writer" = csv.writer(f, delimiter=",")
            tsv_writer.writerow(["col1", "col2", "col3"])
            tsv_writer.writerow([1, 2, 3])
            tsv_writer.writerow([4, 5, 6])
            
        with open("./tmp_pytest/file_read_test_tsv.txt", "w") as f:
            tsv_writer: "_csv._writer" = csv.writer(f, delimiter="\t")
            tsv_writer.writerow([1.1, 2.2, 3.3])
            tsv_writer.writerow([4.4, 5.5, 6.6])
        
    
    @pytest.mark.parametrize("path,delim,col_names,dtype",
                [("./tmp_pytest/file_read_test_csv.txt", ",", True, int),
                 ("./tmp_pytest/file_read_test_tsv.txt", "\t", False, float)]
                )
    def test_load_data_filetype(self, path: str, delim: str, col_names: bool, dtype):
        out_file: List["np.ndarray"] = cytomulate.utilities.FileIO.load_data(path, col_names, delim=delim, dtype = dtype)
        assert isinstance(out_file, list)
        assert isinstance(out_file[0], np.ndarray)
        assert isinstance(out_file[1], np.ndarray)
        assert len(out_file) == 2
        assert out_file[1].shape == (2, 3)
        
        if path == "./tmp_pytest/file_read_test_csv.txt":
            assert "col1" in out_file[0]
            assert out_file[1].dtype == np.dtype("int64")
        else:
            assert None in out_file[0]
            assert out_file[1].dtype == np.dtype("float64")
            
    
    @pytest.mark.parametrize("drop_cols,expected_shape",
            [([0, 1], (2, 1)), (1, (2,2))]
            )      
    def test_load_data_drop_col(self, drop_cols, expected_shape):
        path: str = "./tmp_pytest/file_read_test_csv.txt"
        out_file: List["np.ndarray"] = cytomulate.utilities.FileIO.load_data(path, True, drop_columns=drop_cols, delim=",", dtype = int)
        assert out_file[1].shape == expected_shape
        
        
    @pytest.mark.parametrize("save_list,path",
                    [([[1, 2, 3], [1, 2, 3]], "./tmp_pytest/2d_list_1.csv"),
                    ([[1.4, 2.2], ["a", "b"]], "./tmp_pytest/2d_list_2.csv"),
                    ([[1.4, 1], [1, "a"]], "./tmp_pytest/2d_list_3.csv")]
                    )   
    def test_save_2d_list_to_csv(self, save_list: List[List[Any]], path: str):
        cytomulate.utilities.FileIO.save_2d_list_to_csv(save_list, path)
        assert os.path.exists(path)


    def test_save_2d_list_to_csv_exception(self):
        path: str = "./tmp_pytest/2d_list_3.csv"
        save_list: List[List[Any]] = [[1.4, 1], [1, "a"]]
        try:
            cytomulate.utilities.FileIO.save_2d_list_to_csv(save_list, path)
        except FileExistsError:
            assert True
        else:
            assert False


    def test_save_np_array(self):
        arr: "np.ndarray" = np.array([[2.1, 2.2, 2.3], [1.1, 2.2, 3.3]])
        path: str = "./tmp_pytest/nparr.txt"
        cytomulate.utilities.FileIO.save_np_array(arr, path)
        assert os.path.exists("./tmp_pytest/nparr.txt")
        
        
    def test_save_np_array_colnames(self): 
        arr: "np.ndarray" = np.array([[2.1, 2.2, 2.3], [1.1, 2.2, 3.3]])
        header: "np.ndarray" = np.array(["A", "B"])
        path: str = "./tmp_pytest/nparr_colnames.txt"
        cytomulate.utilities.FileIO.save_np_array(arr, path, col_names=header)
        assert os.path.exists("./tmp_pytest/nparr_colnames.txt")
        
    
    def test_save_np_array_exception(self):
        arr: "np.ndarray" = np.array([[2.1, 2.2, 2.3], [1.1, 2.2, 3.3]])
        path: str = "./tmp_pytest/nparr.txt"
        try:
            cytomulate.utilities.FileIO.save_np_array(arr, path)
        except FileExistsError:
            assert True
        else:
            assert False
        
        
    def test_dir_create(self):
        cytomulate.utilities.FileIO.make_dir("./tmp_pytest/test_create")
        assert os.path.exists("./tmp_pytest/test_create")
        
     
    @pytest.mark.parametrize("path_in,counter,path_out",
                        [("./tmp_pytest/test_recur", 0, "./tmp_pytest/test_recur1"),
                        ("./tmp_pytest/test_recur_counter", 2, "./tmp_pytest/test_recur_counter2")]
                        )   
    def test_dir_create_recursive(self, path_in: str, counter: int, path_out: str):
        os.mkdir(path_in)
        cytomulate.utilities.FileIO.make_dir(path_in, add_number_if_dir_exists=True, _counter=counter)
        assert os.path.exists(path_out)
    
    
    @pytest.mark.parametrize("path_in,recursive,expected",
                    [("./tmp_pytest/test_recur_out", False, "./tmp_pytest/test_recur_out"),
                    ("./tmp_pytest/test_recur_out", True, "./tmp_pytest/test_recur_out1")]
                    )  
    def test_dir_create_return(self, path_in: str, recursive: bool, expected: str):
        path_out: str = cytomulate.utilities.FileIO.make_dir(path_in, add_number_if_dir_exists=recursive)
        assert path_out == expected
    
    
    def test_dir_create_exception(self):
        try:
            cytomulate.utilities.FileIO.make_dir("./tmp_pytest/test_create")
        except FileExistsError:
            assert True
        else:
            assert False
    
    
    @classmethod
    def teardown_class(cls):
        shutil.rmtree("./tmp_pytest/")
