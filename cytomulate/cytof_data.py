# Math computation
import numpy as np
from numpy import random as rd

from cytomulate.forest import Forest


class CytofData:
    def __init__(self, n_batches, n_trees, n_cells_per_batch,\
                 n_cell_types, n_markers):
        self.n_batches = n_batches
        self.n_trees = n_trees
        self.n_cells = n_cells_per_batch
        self.n_cell_types = n_cell_types
        self.n_markers = n_markers

        # We use a Dirichlet distribution to generate
        # cell type proportions for each batch
        self.cell_type_proportions = rd.dirichlet(np.ones(self.n_cell_types),\
                                                  self.n_batches)

        self.cytof_data = {}
        for b in range(self.n_batches):
            self.cytof_data["batch" + str(b)] = np.zeros((self.n_cells[b], self.n_markers))

        self.forest = Forest(self.n_trees, n_cell_types)
