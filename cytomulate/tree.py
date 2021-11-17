# Math computation
from numpy import random as rd

# Tree and path generation
from cytomulate.utilities import generate_random_tree
from cytomulate.utilities import smooth_brownian_bridge


class Tree:
    def __init__(self, id, name, cell_type_list, root):
        self.id = id
        self.name = name
        self.cell_types = cell_type_list
        self.size = len(self.cell_types)
        self.root = root
        self.edges = []
        self.differentiation_paths = []

    def find_cell_type_by_id(self, id):
        for c in self.cell_types:
            if c.id == id:
                return c

    def sketch_tree(self):
        """
        Construct the edges of the tree

        :return: a nested list
        """
        node_ids = [x.id for x in self.cell_types]
        self.edges = generate_random_tree(node_ids)

    def sketch_family_tree(self):
        """
        Use the edges and the root to get parents and children

        """
        # We then use breadth first search
        doing_list = [self.root]
        while len(doing_list) > 0:
            parent_cell = doing_list.pop(0)
            for e in self.edges:
                if parent_cell.id in e:
                    child_cell_id = set(e) - {parent_cell.id}
                    if not (parent_cell.parent is None):
                        child_cell_id = child_cell_id - {parent_cell.parent.id}
                    if len(child_cell_id) == 1:
                        child_cell = self.find_cell_type_by_id(child_cell_id.pop())
                        child_cell.parent = parent_cell
                        doing_list.append(child_cell)
                        parent_cell.children.append(child_cell)

    def grow_tree(self):
        pass

    def visualize_tree(self):
        pass
