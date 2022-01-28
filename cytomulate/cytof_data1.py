# Math computation
import numpy as np
from numpy import random as rd

from sklearn import cluster

from cytomulate.cell_type import CellType
from cytomulate.forest import Forest
from cytomulate.utilities import smooth_brownian_bridge

from numpy.typing import ArrayLike
from typing import Optional, Sequence


class CytofData:
    def __init__(self, n_batches: int, n_trees: int, n_cells_per_batch: Sequence,\
                 n_cell_types: int, n_markers: int, expression_matrix: Optional["np.ndarray"] = None, \
                 cell_type_indicator: Optional[ArrayLike] = None):
        """Class for CyTOF simulation data.
        
        This is the main entry point for CyTOF simulation and the main class used.

        :param n_batches: The number of batches to simulate
        :type n_batches: int
        :param n_trees: The number of cell-differentiation trees
        :type n_trees: int
        :param n_cells_per_batch: The number of cells per batch
        :type n_cells_per_batch: Sequence
        :param n_markers: The number of markers
        :type n_markers: int
        :param expression_matrix: An existing expression matrix for data-based simulation, defaults to None
        :type expression_matrix: np.ndarray, optional
        :param cell_type_indicator: The names of cell types present or to be simulated
        :type cell_type_indicator: ArrayLike, optional
        """

        if expression_matrix is not None:
            self.simulation_mode = "Data"
        else:
            self.simulation_mode = "Model"

        self.background_noise_level = None

        self.expression_matrix = expression_matrix
        self.cell_type_indicator = cell_type_indicator
        self.cell_types = {}


        self.n_markers = n_markers
        self.n_cell_types = n_cell_types

        self.n_batches = n_batches
        self.batch_effect_overall = np.zeros(self.n_batches)

        self.batch_effect_by_cell_types_markers = {}
        for b in range(self.n_batches):
            self.batch_effect_by_cell_types_markers["batch" + str(b)] = np.zeros((self.n_cell_types, self.n_markers))

        self.temporal_effect = []

        self.n_cells_per_batch = n_cells_per_batch

        self.n_trees = n_trees

        self.cell_type_proportions = np.zeros((self.n_batches, n_cell_types))

        self.cytof_data = {}
        for b in range(self.n_batches):
            self.cytof_data["batch" + str(b)] = {}
            self.cytof_data["batch" + str(b)]["cell_type_indices"] = np.zeros(self.n_cells_per_batch[b])
            self.cytof_data["batch" + str(b)]["expression_matrix"] = np.zeros((self.n_cells_per_batch[b], self.n_markers))

        self.forest = None


    def initialize_cell_types(self, **kwargs):
        """ Initialize the cell types for the object based on simulation mode.
        
        If the ``simulation_mode`` is "Data" and if cell types are not provided, this function
        performs a clustering on the data and fills in the cell type indicator to initialize
        the cell types. If cell types are provided, or the simulation mode is "Model",
        this function simply initializes cell type objects. The resulting cell types are stored
        in the ``cell_types`` attribute of the object.
        
        :param **kwargs: Keyword-only arguments passed into ``sklearn.cluster.KMeans`` for clustering.
        """
        if (self.simulation_mode == "Data") and (self.cell_type_indicator is None):
            assert self.expression_matrix is not None
            cluster_model = cluster.KMeans(**kwargs).fit(self.expression_matrix)
            self.n_cell_types = np.max(cluster_model.labels_) + 1
            for id in range(self.n_cell_types):
                self.cell_types[id] = CellType(id = id, name= id, n_markers=self.n_markers)
                self.cell_types[id].fit_model(dat = self.expression_matrix[cluster_model.labels_ == id, :])

        elif (self.simulation_mode == "Data") and (self.cell_type_indicator is not None):
            assert self.expression_matrix is not None
            cell_type_names = np.unique(self.cell_type_indicator)
            self.n_cell_types = len(cell_type_names)
            for id in range(len(cell_type_names)):
                self.cell_types[id] = CellType(id=id, name=cell_type_names[id], n_markers=self.n_markers)
                self.cell_types[id].fit_model(dat=self.expression_matrix[self.cell_type_indicator == cell_type_names[id], :])

        else:
            for id in range(self.n_cell_types):
                self.cell_types[id] = CellType(id = id, name= id, n_markers=self.n_markers)


    def get_markers_pattern_and_background_noise_level(self):
        """
        If the simulation_mode is Data then this function
        uses the results of the clustering to estimate the
        background noise level as well as the markers pattern
        for each cell type
        Otherwise, the user needs to provide one
        :return:
        """
        pass


    def initialize_cell_type_models(self):
        """
        This function uses the markers pattern information
        and the clustering result to fit a Gaussian Mixture
        to each cell type
        The results will be biased. We will adjust them once we
        have the tree structure
        :return:
        """
        pass


    def grow_forest(self):
        """Generate a forest consisting of cell-differentiation trees.
        
        The forest uses ``n_trees``, ``n_cell_types``, and ``n_markers`` 
        to grow a skeleton structure of the dataset.
        """
        self.forest = Forest(self.n_trees, self.n_cell_types, self.n_markers)
        self.forest.assign_cell_types()
        self.forest.sketch_trees()
        self.forest.grow_trees()


    def generate_batch_effect(self):
        """Generate batch effect.
        """
        pass


    def generate_temporal_effect(self):
        """Generate temporal effect.
        
        This method generates a tempral effect uses a smooth brownian bridge for each batch. 
        For each batch, all functions are stored in the ``temporal_effect`` attribute.
        """
        for b in range(self.n_batches):
            temporal_function = smooth_brownian_bridge(0, rd.normal(0, 0.1), N=5, sigma2=0.1)
            self.temporal_effect.append(temporal_function)


    def generate_cell_type_proportions(self):
        """Generate Cell Type Proportions.
        
        This method uses a Dirichlet Distribution to randomly generate cell type
        proportions for each batch. 
        """
        self.cell_type_proportions = rd.dirichlet(np.ones(self.n_cell_types),\
                                                  self.n_batches)


    def grow_leaves(self):
        """Grow leaf nodes from cell-differentiation trees of the forest.
        
        For each batch, leaves are assigned and computed to the trees. The
        results are stored in the ``cytof_data`` dictionary of the object.
        """
        assert self.forest is not None
        for b in range(self.n_batches):
            cell_type_indices = rd.choice(self.n_cell_types, self.n_cells_per_batch[b],\
                                                replace=True, p=self.cell_type_proportions[b,:])
            self.cytof_data["batch" + str(b)]["cell_type_indices"] = cell_type_indices
            expression_matrix = np.zeros((self.n_cells_per_batch[b], self.n_markers))
            for event in range(self.n_cells_per_batch[b]):
                cell = self.forest.find_cell_type_by_id(cell_type_indices[event])
                child_cell = rd.choice(cell.children, 1)
                for marker_id in range(self.n_markers):
                    mu = cell.expression_level[marker_id]
                    sig2 = cell.variance_level[marker_id]
                    pseudo_time = rd.beta(a=0.4, b=1, size=1)

                    mu += child_cell[1](pseudo_time)

                    expression_matrix[event, marker_id] = rd.normal(loc = mu, scale = np.sqrt(sig2), size = 1)

            self.cytof_data["batch" + str(b)]["expression_matrix"] = expression_matrix

