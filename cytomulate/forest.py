# Math computation
import numpy as np
from numpy import random as rd
from scipy.stats import truncnorm

from cytomulate.cell_type import CellType
from cytomulate.tree import Tree


class Forest:
    def __init__(self, n_trees, n_cell_types, n_markers):
        self.n_cell_types = n_cell_types
        self.n_trees = n_trees
        self.n_markers = n_markers

        self.cell_types = [CellType(i, i, n_markers) for i in range(self.n_cell_types)]

        self.L = 4
        self.high_expression = np.cumsum(truncnorm.rvs(0, np.Inf, scale = 0.5, size = self.L))
        self.low_expression = np.cumsum(truncnorm.rvs(0, np.Inf, scale = 0.05, size = self.L-1))
        self.low_expression = np.concatenate(([0], self.low_expression))

        self.trees = []

    def __contains__(self, tree_id):
        pass

    def __str__(self):
        return "Forest"

    def find_cell_type_by_id(self, id):
        for c in self.cell_types:
            if c.id == id:
                return c


    def generate_root_marker_patterns(self):
        p = 0.4
        is_valid = False
        while not is_valid:
            # Each column of the matrix will be the marker
            # pattern of a cell type
            marker_pattern_matrix = rd.binomial(1, p, self.n_markers * self.n_trees).reshape((self.n_markers,\
                                                                                              self.n_trees))
            col_sum = np.sum(marker_pattern_matrix, axis = 0)
            row_sum = np.sum(marker_pattern_matrix, axis = 1)
            temp = np.unique(marker_pattern_matrix, axis = 1)
            if np.any(col_sum > 0) and np.any(row_sum > 0) and (temp.shape[1] == marker_pattern_matrix.shape[1]):
                is_valid = True
        return marker_pattern_matrix

    def assign_cell_types(self):
        marker_pattern_matrix = self.generate_root_marker_patterns()
        nodes = rd.permutation(self.cell_types)
        dividing_points = [self.n_cell_types]
        if self.n_trees > 1:
            temp = rd.choice(self.n_cell_types - 1, self.n_trees - 1, replace=False)
            temp.sort()
            dividing_points = np.concatenate((temp, dividing_points))
        dividing_points = np.concatenate(([0], dividing_points))
        for t in range(self.n_trees):
            cell_type_list = nodes[(dividing_points[t]):(dividing_points[t + 1])]
            root = rd.choice(cell_type_list, 1)
            root.variance_level = np.zeros(self.n_markers)
            root.marker_pattern = marker_pattern_matrix[:,t]
            root.gating_markers = list(rd.choice(self.n_markers, 2, replace = False))
            for marker_id in range(self.n_markers):
                if root.marker_pattern[marker_id] == 0:
                    level = rd.choice(self.low_expression, 1)
                else:
                    level = rd.choice(self.high_expression, 1)

                if level == 0:
                    root.expression_level[marker_id] = 0
                else:
                    root.expression_level[marker_id] = truncnorm.rvs(0, np.Inf, loc = level, scale = 0.01, size = 1)
                    root.variance_level[marker_id] = 1 / rd.gamma(100, 1 / 10, size=1)

            tree = Tree(t, t, cell_type_list, root, self.high_expression, self.low_expression)
            self.trees.append(tree)

    def sketch_trees(self):
        for t in range(self.n_trees):
            self.trees[t].sketch_tree()
            self.trees[t].sketch_family_tree()

    def grow_trees(self):
        for t in range(self.n_trees):
            self.trees[t].grow_tree()

    def visualize_forest(self):
        pass

