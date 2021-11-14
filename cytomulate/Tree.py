# Math computation
from numpy import random as rd

# Tree and path generation
from utilities import generate_random_tree
from utilities import smooth_brownian_bridge


class Tree:
    def __init__(self, id, name, cell_type_list):
        self.id = id
        self.name = name
        self.cell_types = cell_type_list
        self.size = len(self.cell_types)
        self.root = rd.choice(self.cell_types, 1)
        self.edges = []
        self.differentiation_paths = []

    def construct_tree(self):
        """
        Construct a tree structure out of the cell types
        :return: a nested list
        """
        node_ids = [x.id for x in self.cell_types]
        self.edges = generate_random_tree(node_ids)

    def construct_paths(self):
        pass
