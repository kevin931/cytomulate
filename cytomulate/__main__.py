"""The Command-Line Interface (CLI) of Cytomulate

For convenience, users can access cytomulate without the need of writing a script. Consistent
with python CLI syntax, this ``__main__.py`` can be accessed through ``python -m cytomulate``. 
The CLI has the advantages of convenience while also supporting most functionalities of the
packages.

:Example:

    Check version:

    .. code-block:: bash

        python -m cytomulate --version
        
    Get help:
    
    .. code-block:: bash

        python -m cytomulate -h
        
**Advantage**:
    * Shell-scripting friendly
    * Easy-to-use
    * No python knowledge required
    
**Drawbacks**:
    * Not as configurable as a python program
    * Top level functionality and integration only (i.e. No partial package usage)
"""


from cytomulate import CreationCytofData, EmulationCytofData
from cytomulate import __version__
import numpy as np

import os
import sys
import argparse
from typing import Dict, Union

OPT_PCK: Dict[str, bool] = {"PyCytoData": True}

try:
    from PyCytoData import PyCytoData, FileIO
except ImportError:
    OPT_PCK["PyCytoData"] = False


parser: "argparse.ArgumentParser" = argparse.ArgumentParser(description="cytomulate: CyTOF Simulation")
parser.add_argument("--version", action="version", version=__version__)

# Creation Mode
parser.add_argument("--creation", help="Creation mode of Cytomulate.", action="store_true")
parser.add_argument("--n_cells", type=int, help="The number of cells to simulate")
parser.add_argument("--n_batches", type=int, help="The number of batches", default=1)
parser.add_argument("--n_types", type=int, help="The number of cell types to simulate.", default=10)
parser.add_argument("--n_markers", type=int, help="The number of markers to simulate.", default=20)
parser.add_argument("--n_trees", type=int, help="The number of differentiation trees to simulate.", default=2)


# Emulation Mode
parser.add_argument("--emulation", help="Creation mode of Cytomulate.", action="store_true")
parser.add_argument("--exprs", help="The path to existing expression matrix.")
parser.add_argument("--cell_types", help="The path to existing cell types.")
parser.add_argument("--exprs_colnames", help="Whether the first row of the existing expression matrix is column names.", action="store_true")
parser.add_argument("--cell_types_colnames", help="Whether the first row of the existing cell types is column names.", action="store_true")
parser.add_argument("--exprs_delim", help="The delimiter of existing expression matrix.", type=str, default="\t")
parser.add_argument("--cell_types_delim", help="The delimiter of existing cell types.", type=str, default="\t")

# Complex Simulation
parser.add_argument("--trajectory", help="Add cell differentiation trajectory.", action="store_true")
parser.add_argument("--temporal_effect", help="Add temporal effect to the dataset. Only the Brownian Bridge method is supported.", action="store_true")
parser.add_argument("--temporal_effect_var", help="Add cell differentiation trajectory.", type=float, default=0.1)
parser.add_argument("--batch_effect", help="Add batch effect to the dataset.", action="store_true")
parser.add_argument("--batch_effect_var", help="Add cell differentiation trajectory.", type=float, default=0.1)

# Output
parser.add_argument("-o", "--out_dir", help="Directory to save simulation data files", type=str)
parser.add_argument("--make_new_dir", help="Directory to save simulation data files", action="store_true")


def main(cmdargs: argparse.Namespace):
    """The main method for cytomulate.

    This is the command-line driver method: it supports both creation
    and emulation mode. To get help, run:
    
    .. code-block::

        python -m cytomulate -h
    

    :param cmdargs: The command line arguments and flags.
    :type cmdargs: argparse.Namespace
    """
    
    if not OPT_PCK["PyCytoData"]:
        print("`PyCytoData` is necessary for CLI. To install, run 'pip install PyCytoData'.")
        sys.exit(1)
            
    if cmdargs.emulation and cmdargs.creation:
        raise ValueError("Cannot run 'Emulation Mode' and  'Creation Mode' simultaneously.")
    
    if cmdargs.make_new_dir:
        os.makedirs(cmdargs.out_dir)
    
    cmdargs.out_dir.rstrip("/")
    cmdargs.out_dir.rstrip("\\")
    
    exprs_path: str = cmdargs.out_dir + "/exprs.txt"
    cell_types_path: str = cmdargs.out_dir + "/cell_types.txt"
    sample_index_path: str = cmdargs.out_dir + "/sample_index.txt"
    
    cytof_data: Union[CreationCytofData, EmulationCytofData]
    
    if cmdargs.creation:
        cytof_data = CreationCytofData(n_batches=cmdargs.n_batches, n_types=cmdargs.n_types, n_markers=cmdargs.n_markers, n_trees=cmdargs.n_trees)
        cytof_data.initialize_cell_types()
        if cmdargs.temporal_effect:
            cytof_data.generate_temporal_effects(variance=cmdargs.temporal_effect_var)
        if cmdargs.batch_effect:
            cytof_data.generate_overall_batch_effects(variance=cmdargs.batch_effect_var)
            cytof_data.generate_local_batch_effects(variance=cmdargs.batch_effect_var)
        if cmdargs.trajectory:
            cytof_data.generate_cell_graph()
        exprs: PyCytoData = cytof_data.sample_to_pycytodata(n_samples = cmdargs.n_cells)
        FileIO.save_np_array(exprs.expression_matrix, path=exprs_path, col_names=exprs.channels)
        FileIO.save_np_array(exprs.cell_types, path=cell_types_path, dtype="%s")
        FileIO.save_np_array(exprs.sample_index, path=sample_index_path, dtype="%s")
          
    if cmdargs.emulation:
        exprs_skiprow: int = 1 if cmdargs.exprs_colnames else 0
        cell_types_skiprow: int = 1 if cmdargs.cell_types_colnames else 0
            
        expression_matrix: Union[np.ndarray, tuple] = FileIO.load_delim(cmdargs.exprs, exprs_skiprow, delim=cmdargs.exprs_delim)
        cell_types: Union[np.ndarray, tuple] = FileIO.load_delim(cmdargs.cell_types, cell_types_skiprow, delim=cmdargs.cell_types_delim, dtype=str)
        assert isinstance(expression_matrix, np.ndarray)
        assert isinstance(cell_types, np.ndarray)
        
        cytof_data = EmulationCytofData(n_batches=cmdargs.n_batches)
        cytof_data.initialize_cell_types(expression_matrix=expression_matrix,
                                         labels=cell_types)
        if cmdargs.temporal_effect:
            cytof_data.generate_temporal_effects(variance=cmdargs.temporal_effect_var)
        if cmdargs.batch_effect:
            cytof_data.generate_overall_batch_effects(variance=cmdargs.batch_effect_var)
            cytof_data.generate_local_batch_effects(variance=cmdargs.batch_effect_var)
        if cmdargs.trajectory:
            cytof_data.generate_cell_graph()
        exprs: PyCytoData = cytof_data.sample_to_pycytodata(n_samples = cmdargs.n_cells)
        FileIO.save_np_array(exprs.expression_matrix, path=exprs_path, col_names=exprs.channels)
        FileIO.save_np_array(exprs.cell_types, path=cell_types_path, dtype="%s")
        FileIO.save_np_array(exprs.sample_index, path=sample_index_path, dtype="%s")
        
    sys.exit(0)


if __name__ == "__main__":
    cmdargs: argparse.Namespace = parser.parse_args()
    main(cmdargs=cmdargs)