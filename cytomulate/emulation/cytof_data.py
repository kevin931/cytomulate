
from tqdm import tqdm

import numpy as np
from collections import Counter
from copy import deepcopy
from emulation.cell_type import EmulationCellType
from emulation.cell_graph import EmulationCellGraph

# Superclass
from cytof_data_general import GeneralCytofData


class EmulationCytofData(GeneralCytofData):
    def __init__(self, n_batches=1, background_noise_model=None, bead_label=None):
        super().__init__(n_batches, background_noise_model)
        """Initializer of CytofData
        Parameters
        ----------
        n_batches: number of batches
        """

        self.bead_label = bead_label

        self.observed_cell_abundances = {}

        self.cell_graph = EmulationCellGraph()

    def initialize_cell_types(self, expression_matrix,
                              labels,
                              max_components=9,
                              min_components=1,
                              covariance_types=("full", "tied", "diag", "spherical")):
        """Initialize cell type models by fitting Gaussian mixtures
        
        Parameters
        ----------
        expression_matrix: A matrix containing the expression levels of cell events
        labels: A vector of cell type labels
        max_components: Used for Gaussian mixture model selection. The maximal number of components for a Gaussian mixture
        min_components: Used for Gaussian mixture model selection. The minimal number of components for a Gaussian mxitrue
        covariance_types: Used for Gaussian mixture model selection. The candidate types of covariances

        Returns
        -------

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

    def generate_cell_graph(self, graph_topology = "forest", **kwargs):
        self.cell_graph.initialize_graph(self.cell_types, bead_label=self.bead_label)
        self.cell_graph.prune_graph(graph_topology)
        self.cell_graph.generate_trajectories(self.cell_types, **kwargs)

    def generate_cell_abundances(self, use_observed=True, is_random=True):
        if use_observed:
            for b in range(self.n_batches):
                self.cell_abundances[b] = deepcopy(self.observed_cell_abundances)
        else:
            super().generate_cell_abundances(is_random)
