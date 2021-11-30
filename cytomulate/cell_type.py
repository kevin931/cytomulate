# Math computation
import numpy as np


class CellType:
    def __init__(self, id, name, n_markers):
        self.id = id
        self.name = name


        self.parent = None
        self.children = []
        self.marker_pattern = np.zeros(n_markers)
        self.expression_level = np.zeros(n_markers)
        self.variance_level = np.zeros(n_markers)

        # parent_cell_type will be a dictionary on size 1
        # whose key is the id of the parent
        # and whose value is the actual CellType
        # object of the parent
        self.parent_cell_type = {}

        # children_cell_types will be a dictionary
        # whose keys are the ids of the children
        # and whose values are the actual
        # CellType objects of the children
        self.children_cell_types = {}

        # differential_paths_to_children will be
        # a dictionary whose keys are the ids of
        # the children and whose values are the
        # collections of functions of
        # generated differential paths
        self.differential_paths_to_children = {}

        # We expect the markers_pattern to be a list of 1's and 0's
        self.markers_pattern = np.zeros(n_markers)

        # We will use Gaussian Mixture to model or learn
        # the expression of this CellType
        # Consequently, there will be no need for
        # us to store expression levels and
        # variance levels separately
        self.model_for_expressed_markers = []

        self.gating_markers = []

    def generate_from_paths(self, child_id):
        pass 

    def generate_from_model(self):
        pass

    def generate_expressions(self, differentiate = True):
        pass