# Math computation
import numpy as np
from numpy import random as rd

from cytomulate.cell_type import CellType
from cytomulate.tree import Tree


class Forest:
    def __init__(self, n_trees, n_cell_types, n_markers):
        self.n_cell_types = n_cell_types
        self.n_trees = n_trees
        self.n_markers = n_markers

        self.cell_types = [CellType(i, i, n_markers) for i in range(self.n_cell_types)]

        self.trees = []

    def __contains__(self, tree_id):
        pass

    def __str__(self):
        return "Forest"

    def assign_cell_types(self):
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
            tree = Tree(t, t, cell_type_list, root)
            self.trees.append(tree)

    def sketch_trees(self):
        for t in range(self.n_trees):
            self.trees[t].sketch_tree()

    def grow_trees(self):
        for t in range(self.n_trees):
            self.trees[t].grow_tree()

    def visualize_forest(self):
        pass

