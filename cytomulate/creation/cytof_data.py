from tqdm import tqdm

import numpy as np
from scipy.stats import truncnorm

from creation.cell_type import CreationCellType
from creation.cell_graph import CreationCellGraph

# Superclass
from cytof_data_general import GeneralCytofData


class CreationCytofData(GeneralCytofData):
    def __init__(self, n_batches=1, n_types=10, n_markers=20,
                 n_trees=2, background_noise_model=None):
        super().__init__(n_batches, background_noise_model)

        self.n_markers = n_markers

        self.cell_graph = CreationCellGraph()

        labels = np.arange(n_types)

        cell_id = 0
        for c_type in labels:
            self.cell_type_labels_to_ids[c_type] = cell_id
            self.cell_type_ids_to_labels[cell_id] = c_type
            self.cell_types[c_type] = CreationCellType(label=c_type, cell_id=cell_id, n_markers=n_markers)
            cell_id += 1

        self.cell_graph.initialize_graph(self.cell_types, n_trees)

    def initialize_cell_types(self, L=4, scale=0.5, n_components = 5):
        high_expressions = np.cumsum(truncnorm.rvs(a=0, b=np.inf, loc=0, scale=scale, size = L))
        low_expressions = np.cumsum(truncnorm.rvs(a=0, b=np.inf, loc=0, scale=scale/10, size=L-1))
        low_expressions = np.append(0, low_expressions)

        for c_type in tqdm(self.cell_graph.serialized_graph):
            self.cell_types[c_type].generate_marker_expression_patterns(self.cell_types, self.cell_graph.graph)
            self.cell_types[c_type].generate_marker_expressions(self.cell_types, self.cell_graph.graph,
                                                                high_expressions, low_expressions)
            self.cell_types[c_type].generate_model(n_components)

    def generate_cell_graph(self, **kwargs):
        self.cell_graph.generate_trajectories(self.cell_types, **kwargs)
