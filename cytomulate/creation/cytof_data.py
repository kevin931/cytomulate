# Progress bar
from tqdm import tqdm

# Math
import numpy as np
from scipy.stats import truncnorm

# Classes to use
from cytomulate.creation.cell_type import CreationCellType
from cytomulate.creation.cell_graph import CreationCellGraph

# Superclass
from cytomulate.cytof_data_general import GeneralCytofData

# Typing
from typing import Optional, Callable, Union


class CreationCytofData(GeneralCytofData):
    def __init__(self,
                 n_batches: int = 1,
                 n_types: int = 10,
                 n_markers: int = 20,
                 n_trees: int = 2,
                 background_noise_model: Optional[Union[Callable, dict]] = None) -> None:
        """Initialize the CreationCytofData object

        Parameters
        ----------
        n_batches: int
            Number of batches to be simulated
        n_types: int
            Number of cell types to be simulated
        n_markers: int
            Number of markers (columns) to be used
        n_trees: int
            Number of trees in the cell graph
        background_noise_model: Callable or dict
            The function used to generate background noise. It should only take one input: size
        """
        super().__init__(n_batches, background_noise_model)

        self.n_markers = n_markers

        self.cell_graph = CreationCellGraph()

        labels = np.arange(n_types)

        # We first initialize all the cell types
        cell_id = 0
        for c_type in labels:
            self.cell_type_labels_to_ids[c_type] = cell_id
            self.cell_type_ids_to_labels[cell_id] = c_type
            self.cell_types[c_type] = CreationCellType(label=c_type, cell_id=cell_id, n_markers=n_markers)
            cell_id += 1

        # Since the simulation of markers and expressions depends on the cell graph
        # we initialize the cell graph here
        self.cell_graph.initialize_graph(self.cell_types, n_trees)

    def initialize_cell_types(self,
                              L: int = 4,
                              scale: float = 0.5,
                              n_components: int = 1,
                              variance_mode: float = 0.01) -> None:
        """Initialize cell type objects

        Parameters
        ----------
        L: int
            Number of levels of expressions
        scale: float
            The scale parameter used in generating expression levels
        n_components: int
            Number of components in a GMM
        variance_mode: float
            The mode of the variance of the inverse wishart distribution
        """
        # We first generate high expression levels and low expression levels
        # Truncated normals are used to ensure the ordering
        high_expressions = np.cumsum(truncnorm.rvs(a=0, b=np.inf, loc=0, scale=scale, size = L))
        low_expressions = np.cumsum(truncnorm.rvs(a=0, b=np.inf, loc=0, scale=scale/10, size=L-1))
        low_expressions = np.append(0, low_expressions)

        for c_type in tqdm(self.cell_graph.serialized_graph):
            self.cell_types[c_type].generate_marker_expression_patterns(self.cell_types, self.cell_graph.graph)
            self.cell_types[c_type].generate_marker_expressions(self.cell_types, self.cell_graph.graph,
                                                                high_expressions, low_expressions)
            self.cell_types[c_type].generate_model(n_components, variance_mode)

    def generate_cell_graph(self, **kwargs) -> None:
        """Generate cell differentiation paths

        Parameters
        ----------
        kwargs:
            Parameters used for trajectory generation
        """
        self.cell_graph.generate_trajectories(self.cell_types, **kwargs)
