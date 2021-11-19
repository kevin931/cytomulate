# Math computation
from numpy import random as rd
from scipy.stats import truncnorm

# List manipulation
from copy import deepcopy

# Tree and path generation
from cytomulate.utilities import generate_random_tree
from cytomulate.utilities import smooth_brownian_bridge


class Tree:
    def __init__(self, id, name, cell_type_list, root, high_expression, low_expression):
        self.id = id
        self.name = name
        self.cell_types = cell_type_list
        self.size = len(self.cell_types)
        self.root = root
        self.n_markers = len(self.root.marker_pattern)
        self.high_expression = high_expression
        self.low_expression = low_expression
        self.edges = []
        # self.differentiation_paths = []

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
                        parent_cell.children.append([child_cell, []])

    def inherit_marker_patterns(self):
        # We will again use BFS to grow the tree
        p = 0.2
        doing_list = [self.root]
        while len(doing_list) > 0:
            parent_cell = doing_list.pop(0)
            for child_list in parent_cell.children:
                child_cell = child_list[0]
                child_cell.marker_pattern = deepcopy(parent_cell.marker_patter)
                child_cell.gating_markers = deepcopy(parent_cell.gating_markers)

                fickle_markers = set(range(self.n_markers)) - set(child_cell.gating_markers)
                flip_markers = set(rd.choice(list(fickle_markers), int(p * len(fickle_markers))))

                counter = 1
                while len(fickle_markers) > 0:
                    marker_id = fickle_markers.pop()
                    if marker_id in flip_markers:
                        child_cell.marker_pattern[marker_id] = -1 * child_cell.marker_pattern[marker_id] + 1
                        if counter <= 2:
                            child_cell.gating_markers.append(marker_id)
                            counter += 1
                doing_list.append(child_cell)

    def inherit_expression_variance_levels(self):
        doing_list = [self.root]
        while len(doing_list) > 0:
            parent_cell = doing_list.pop(0)
            for child_list in parent_cell.children:
                child_cell = child_list[0]
                child_cell.expression_level = deepcopy(parent_cell.expression_level)
                child_cell.variance_level = deepcopy(parent_cell.variance_level)
                for marker_id in range(self.n_markers):

                    if marker_id in parent_cell.gating_markers:
                        continue

                    if child_cell.marker_pattern[marker_id] == 0:
                        level = rd.choice(self.low_expression, 1)
                    else:
                        level = rd.choice(self.high_expression, 1)

                    if level == 0:
                        child_cell.expression_level[marker_id] = 0
                        child_cell.variance_level[marker_id] = 0
                    else:
                        child_cell.expression_level[marker_id] = truncnorm.rvs(0, np.Inf, loc=level, scale=0.01, size=1)
                        child_cell.variance_level[marker_id] = 1 / rd.gamma(100, 1 / 10, size=1)
                doing_list.append(child_cell)

    def grow_branches(self):
        doing_list = [self.root]
        while len(doing_list) > 0:
            parent_cell = doing_list.pop(0)
            for child_list in parent_cell.children:
                child_cell = child_list[0]
                for marker_id in range(self.n_markers):
                    if parent_cell.expression_level[marker_id] == child_cell.expression_level[marker_id]:
                        child_list[1].append(None)
                    else:
                        child_list[1].append(smooth_brownian_bridge(0,\
                                               child_cell.expression_level[marker_id] - \
                                               parent_cell.expression_level[marker_id],N = 5,\
                                               sigma2= 0.1))
                doing_list.append(child_cell)

    def grow_tree(self):
        self.inherit_marker_patterns()
        self.inherit_expression_variance_levels()
        self.grow_branches()

    def visualize_tree(self):
        pass
