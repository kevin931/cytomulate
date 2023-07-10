# Progress bar
from tqdm import tqdm

# List manipulation
from copy import deepcopy
import numpy as np
from collections import Counter

# Classes to include
from cytomulate.emulation.cell_type import EmulationCellType
from cytomulate.emulation.cell_graph import EmulationCellGraph

# Superclass
from cytomulate.cytof_data_general import GeneralCytofData

# Typing
from typing import Union, Optional, Callable, Tuple, List


class EmulationCytofData(GeneralCytofData):
    def __init__(self,
                 n_batches: int = 1,
                 background_noise_model: Optional[Union[Callable, dict]] = None,
                 bead_label: Optional[Union[str, int]] = None) -> None:
        """The Emulation Mode object for Cytomulate.
        
        This class serves as a starting point for the Emulation Mode of Cytomulate. The constructor
        defines the key parameters of the simulation, including the number of batches. Unlike the
        Creation mode, other parameters such as the number of protein markers are fixed from the
        dataset rather than user-soecified. The number of cells is defined later at a sampling step.

        Parameters
        ----------
        n_batches: int
            The number of batches to be simulated
        background_noise_model: Callable or dict
            The model used to generate random values. It should have only one input: size
        bead_label: str or int
            The label for beads
        """
        super().__init__(n_batches, background_noise_model)

        self.bead_label = bead_label

        self.observed_cell_abundances = {}

        self.cell_graph = EmulationCellGraph()

    def initialize_cell_types(self,
                              expression_matrix: np.ndarray,
                              labels: np.ndarray,
                              max_components: int = 9,
                              min_components: int = 1,
                              covariance_types: Union[List[str], Tuple[str]] = ("full", "tied", "diag", "spherical")) -> None:
        """Initialize cell type models by fitting Gaussian mixtures
        
        This method fits the GMM models for each cell type. Namely, a Gaussian Mixture Model
        is generated for each cell type at this stage according to the parameters specified.
        An extensive model selection procedure based on the Bayesian Information Criterion (BIC)
        is performed when multiple possibilities of components and covariance types are
        specified. See details in `max_components` and `covariance_types`.
        
        Parameters
        ----------
        expression_matrix: np.ndarray
            A matrix containing the expression levels of cell events
        labels: np.ndarray
            A vector of cell type labels
        max_components: int
            The maximal number of components for a Gaussian mixture. Used for Gaussian mixture model selection.
            This must be smaller or equal to the `max_components`. If `max_components` equals `min_components`,
            the exact number will be used for fitting. Otherwise, a model selection procedure will ensue using
            Bayesian Information Criterion.
        min_components: int
            The minimal number of components for a Gaussian mxitrue. Used for Gaussian mixture model selection.
            This must be smaller or equal to the `max_components`. See `max_components` for details on model
            selection.
        covariance_types: list or tuple
            The candidate types of covariances used for Gaussian mixture model selection. If only one is specified,
            no model selection will be performed based on the covariance structure.
        """
        self.n_markers = np.shape(expression_matrix)[1]

        unique_labels = np.unique(labels)

        abundances = Counter(labels)

        cell_id = 0
        for c_type in tqdm(unique_labels):
            self.observed_cell_abundances[c_type] = abundances[c_type]/len(labels)

            self.cell_type_labels_to_ids[c_type] = cell_id
            self.cell_type_ids_to_labels[cell_id] = c_type

            self.cell_types[c_type] = EmulationCellType(label=c_type, cell_id=cell_id, n_markers=self.n_markers)

            ind = np.where(labels == c_type)[0]
            D = expression_matrix[ind, :]

            self.cell_types[c_type].fit(data=D,
                                        max_components=max_components,
                                        min_components=min_components,
                                        covariance_types=covariance_types)

            cell_id += 1


    def generate_cell_graph(self,
                            graph_topology: str = "forest",
                            **kwargs) -> None:
        """Generate a cell graph as well as differentiation paths
        
        This method is part of complex simulation's cellular trajectory simulation. It
        generates differentiation paths, which will be used at the sampling stage.

        Parameters
        ----------
        graph_topology: str
            Type of graph to be generated
        kwargs:
            Other parameters used for trajectory generation
        """
        self.cell_graph.initialize_graph(self.cell_types, bead_label=self.bead_label)
        self.cell_graph.prune_graph(graph_topology)
        self.cell_graph.generate_trajectories(self.cell_types, **kwargs)


    def generate_cell_abundances(self,
                                 use_observed: bool = True,
                                 is_random: bool = True) -> None:
        """Generate cell abundances

        Generate the cell abundances for all cell types: namely, the amount
        of cells in each cell type. This method supports either data-based
        cell abundance or randomly-generated cell abundance. In the latter
        case, each cell type's probability can be further randomized.

        Parameters
        ----------
        use_observed: bool
            Whether the cell abundances should use the observed ones
        is_random: bool
            In the case that `user_obsersed` is `False`, whether the cell abundances'
            probability should be randomly generated. If `True`, the abundance of each
            cell type is sampled from a dirichlet distribution. If `False`, then all cell
            types an have equal probability.
            
        Note
        -----
        If you wish to use the default observed cell abundance from the data,
        it is not necessary to call this method. Otherwise, you should always
        set ``used_observed`` to ``False``.
        """
        if use_observed:
            for b in range(self.n_batches):
                self.cell_abundances[b] = deepcopy(self.observed_cell_abundances)
        else:
            super().generate_cell_abundances(is_random)
