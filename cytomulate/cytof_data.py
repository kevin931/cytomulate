# Math computation
import numpy as np
from numpy import random as rd

from cytomulate.forest import Forest
from cytomulate.utilities import smooth_brownian_bridge


class CytofData:
    def __init__(self, n_batches, n_trees, n_cells_per_batch,\
                 n_cell_types, n_markers):
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
            self.cytof_data["batch" + str(b)] = np.zeros((self.n_cells_per_batch[b], self.n_markers))

        self.forest = Forest(self.n_trees, self.n_cell_types, self.n_markers)

    def generate_batch_effect(self):
        pass

    def generate_temporal_effect(self):
        for b in range(self.n_batches):
            temporal_function = smooth_brownian_bridge(0, rd.normal(0, 0.1), N=5, sigma2=0.1)
            self.temporal_effect.append(temporal_function)

    def generate_cell_type_proportions(self):
        # We use a Dirichlet distribution to generate
        # cell type proportions for each batch
        self.cell_type_proportions = rd.dirichlet(np.ones(self.n_cell_types),\
                                                  self.n_batches)

    def grow_forest(self):
        self.forest.assign_cell_types()
        self.forest.sketch_trees()
        self.forest.grow_trees()

    def grow_leaves(self):
        pass