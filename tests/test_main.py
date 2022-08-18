from cytomulate import __main__
import cytomulate
import pytest
import numpy as np

import sys
import os
import shutil
from io import StringIO

from typing import Dict, Union

OPT_PCK: Dict[str, bool] = {"PyCytoData": True}
try:
    from PyCytoData import PyCytoData, FileIO
except ImportError:
    OPT_PCK["PyCytoData"] = False


def test_parser_version():
    screen_stdout = sys.stdout
    string_stdout = StringIO()
    sys.stdout = string_stdout
    
    try:
        __main__.parser.parse_args(["--version"])
    except SystemExit:
        output = string_stdout.getvalue()
        expected = cytomulate.__version__ + "\n"
        assert output == expected
        sys.stdout = screen_stdout
    else:
        sys.stdout = screen_stdout
        assert False
        

@pytest.mark.parametrize("cmdarg", ["--help", "-h"])
def test_parser_help(cmdarg):
    screen_stdout = sys.stdout
    string_stdout = StringIO()
    sys.stdout = string_stdout
    
    try:
        __main__.parser.parse_args([cmdarg])
    except SystemExit:
        output = string_stdout.getvalue()
        assert "usage: " in output
        sys.stdout = screen_stdout
    else:
        sys.stdout = screen_stdout
        assert False
        

if OPT_PCK["PyCytoData"]:
    class TestMain():
        
        @classmethod
        def setup_class(cls):
            os.mkdir("./tmp_pytest/")
            exprs: np.ndarray = np.random.rand(100, 10)
            cell_types: np.ndarray = np.repeat(["a", "b"], 50)
            
            np.savetxt("./tmp_pytest/exprs.txt", exprs, delimiter="\t")
            np.savetxt("./tmp_pytest/cell_types.txt", cell_types, delimiter="\t", fmt="%s")
            
        
        @pytest.mark.parametrize("n_batches,n_markers,n_types,n_trees,trajectory,batch_effect,batch_effect_var,temporal_effect,temporal_effect_var,count",
                                [(1, 20, 10, 2, False, False, 0.1, False, 0.1, 1),
                                (2, 30, 5, 1, False, True, 0.1, False, 0.1, 2),
                                (1, 40, 20, 3, True, False, 0.1, True, 0.1, 3),
                                (1, 40, 20, 3, True, True, 0.1, True, 0.1, 4)])
        def test_creation_mode(self, mocker, n_batches: int, n_markers: int, n_types: int, n_trees: int, trajectory: bool,
                               batch_effect: bool, batch_effect_var: float,temporal_effect: bool, temporal_effect_var: float, count):
            os.mkdir("./tmp_pytest/creation{}/".format(count))
            cmdargs: mocker.MagicMock = mocker.MagicMock()
            cmdargs.creation = True
            cmdargs.emulation = False
            cmdargs.out_dir = "./tmp_pytest/creation{}/".format(count)
            cmdargs.n_cells = 1000
            cmdargs.n_batches = n_batches
            cmdargs.n_markers = n_markers
            cmdargs.n_types = n_types
            cmdargs.n_trees = n_trees
            cmdargs.make_new_dir = False
            cmdargs.trajectory = trajectory
            cmdargs.temporal_effect = temporal_effect
            cmdargs.temporal_effect_var = temporal_effect_var
            cmdargs.batch_effect = batch_effect
            cmdargs.batch_effect_var = batch_effect_var
            
            try:
                __main__.main(cmdargs)
            except SystemExit:
                assert os.path.exists("./tmp_pytest/creation{}/exprs.txt".format(count))
                assert os.path.exists("./tmp_pytest/creation{}/cell_types.txt".format(count))
                assert os.path.exists("./tmp_pytest/creation{}/sample_index.txt".format(count))
                
                exprs: PyCytoData = FileIO.load_expression("./tmp_pytest/creation{}/exprs.txt".format(count), col_names=True)
                cell_types: Union[np.ndarray, tuple] = FileIO.load_delim("./tmp_pytest/creation{}/cell_types.txt".format(count), dtype=str)
                sample_index: Union[np.ndarray, tuple] = FileIO.load_delim("./tmp_pytest/creation{}/sample_index.txt".format(count), dtype=str)
                assert isinstance(sample_index, np.ndarray)
                assert isinstance(cell_types, np.ndarray)
                exprs.cell_types = cell_types
                exprs.sample_index = sample_index
                
                
                assert exprs.n_cells == 1000*n_batches
                assert exprs.n_channels == n_markers
                assert exprs.n_cell_types <= n_types
                assert exprs.n_samples == n_batches
                
        
        @pytest.mark.parametrize("colnames,n_batches,trajectory,batch_effect,batch_effect_var,temporal_effect,temporal_effect_var,count",
                                [(False,1, False, False, 0.1, False, 0.1, 1),
                                (False,2, False, False, 0.1, True, 0.1, 2),
                                (True,1, True, True, 0.1, True, 0.2, 3),
                                (False,1, False, True, 0.2, False, 0.1, 4)])  
        def test_emulation_mode(self, mocker, colnames: bool, n_batches: int, trajectory: bool,
                                batch_effect: bool, batch_effect_var: float,temporal_effect: bool, temporal_effect_var: float, count: int):
            os.mkdir("./tmp_pytest/emulation{}/".format(count))
            cmdargs: mocker.MagicMock = mocker.MagicMock()
            cmdargs.creation = False
            cmdargs.emulation = True
            cmdargs.out_dir = "./tmp_pytest/emulation{}/".format(count)
            cmdargs.exprs = "./tmp_pytest/exprs.txt"
            cmdargs.exprs_colnames = colnames
            cmdargs.exprs_delim = "\t"
            cmdargs.cell_types = "./tmp_pytest/cell_types.txt"
            cmdargs.cell_types_colnames = colnames
            cmdargs.cell_types_delim = "\t"
            cmdargs.n_cells = 1000
            cmdargs.n_batches = n_batches
            cmdargs.make_new_dir = False
            cmdargs.trajectory = trajectory
            cmdargs.temporal_effect = temporal_effect
            cmdargs.temporal_effect_var = temporal_effect_var
            cmdargs.batch_effect = batch_effect
            cmdargs.batch_effect_var = batch_effect_var
            
            try:
                __main__.main(cmdargs)
            except SystemExit:
                assert os.path.exists("./tmp_pytest/emulation{}/exprs.txt".format(count))
                assert os.path.exists("./tmp_pytest/emulation{}/cell_types.txt".format(count))
                assert os.path.exists("./tmp_pytest/emulation{}/sample_index.txt".format(count))
                
                exprs: PyCytoData = FileIO.load_expression("./tmp_pytest/emulation{}/exprs.txt".format(count), col_names=True)
                cell_types: Union[np.ndarray, tuple] = FileIO.load_delim("./tmp_pytest/emulation{}/cell_types.txt".format(count), dtype=str)
                sample_index: Union[np.ndarray, tuple] = FileIO.load_delim("./tmp_pytest/emulation{}/sample_index.txt".format(count), dtype=str)
                assert isinstance(sample_index, np.ndarray)
                assert isinstance(cell_types, np.ndarray)
                exprs.cell_types = cell_types
                exprs.sample_index = sample_index
                
                
                assert exprs.n_cells == 1000*n_batches
                assert exprs.n_channels == 10
                assert exprs.n_cell_types <= 2
                assert exprs.n_samples == n_batches
                
                
        def test_emulation_creation_error(self, mocker):
            cmdargs: mocker.MagicMock = mocker.MagicMock()
            cmdargs.creation = True
            cmdargs.emulation = True
            
            try:
                __main__.main(cmdargs)
            except ValueError as e:
                assert "Cannot run 'Emulation Mode' and  'Creation Mode' simultaneously." in str(e)
                
                
        def test_make_new_dir(self, mocker):
            cmdargs: mocker.MagicMock = mocker.MagicMock()
            cmdargs.creation = True
            cmdargs.emulation = False
            cmdargs.out_dir = "./tmp_pytest/creation5/"
            cmdargs.n_cells = 1000
            cmdargs.n_batches = 1
            cmdargs.n_markers = 20
            cmdargs.n_types = 5
            cmdargs.n_trees = 2
            cmdargs.make_new_dir = True
            cmdargs.trajectory = False
            cmdargs.temporal_effect = False
            cmdargs.batch_effect = False
            try:
                __main__.main(cmdargs)
            except SystemExit:
                assert os.path.exists("./tmp_pytest/creation5/exprs.txt")
                assert os.path.exists("./tmp_pytest/creation5/cell_types.txt")
                assert os.path.exists("./tmp_pytest/creation5/sample_index.txt")
            else:
                assert False
            
            
        @classmethod
        def teardown_class(cls):
            shutil.rmtree("./tmp_pytest/")