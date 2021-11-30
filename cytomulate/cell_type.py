# Math computation
import numpy as np
from numpy import random as rd

# List manipulation
from copy import deepcopy

# Gaussian Mixtures
from sklearn.mixture import GaussianMixture


class CellType:
    def __init__(self, id, name, n_markers):
        self.id = id
        self.name = name
        self.n_markers = n_markers

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
        self.markers_pattern = np.zeros((1, n_markers))

        # We will use Gaussian Mixture to model or learn
        # the expression of this CellType
        # Consequently, there will be no need for
        # us to store expression levels and
        # variance levels separately
        # We expect the class of this field
        # to be sklearn.mixture._gaussian_mixture.GaussianMixture
        self.model_for_expressed_markers = None

        # Gating markers will be a set of indices
        self.gating_markers = {}

    def inherit_markers_pattern(self, mutation_probability = 0.2, \
                                n_additional_gating_markers = 2):
        if len(self.parent_cell_type) == 0:
            # if this is the root
            raise Exception("This is the root.")
        # Get the parent
        parent = list(self.parent_cell_type.values())[0]
        # simply copy the information needed
        self.markers_pattern = deepcopy(parent.markers_pattern)
        self.gating_markers = deepcopy(parent.gating_markers)

        # fickle markers are indices of the markers that can change
        # expression levels or even from 1 to 0
        fickle_markers = set(range(self.n_markers)) - self.gating_markers
        # flip markers are fickle markers that will flip from 0 to 1 and vice versa
        flip_markers = set(rd.choice(list(fickle_markers), int(mutation_probability * len(fickle_markers))))
        for m in flip_markers:
            self.markers_pattern[0, m] = -1 * self.markers_pattern[0, m] + 1
        new_gating_markers = set(rd.choice(list(fickle_markers), n_additional_gating_markers))
        self.gating_markers = self.gating_markers.union(new_gating_markers)

    def inherit_model(self):
        if len(self.parent_cell_type) == 0:
            # if this is the root
            raise Exception("This is the root.")
        # Get the parent
        parent = list(self.parent_cell_type.values())[0]

        self.model_for_expressed_markers = deepcopy(parent.model_for_expressed_markers)

    def fit_model(self, dat = None):
        self.model_for_expressed_markers = GaussianMixture(n_components = 4, random_state = 0).fit(dat)

    def generate_from_paths(self, child_id, alpha = 0.4, beta = 1):
        """
        Generates an array of differential points that can be later on added to the true expression
        :param child_id: ID of the child to differentiate to
        :param alpha: alpha parameter in the Beta distribution
        :param beta: beta parameter in the Beta distribution
        :return: an array of differential points
        """
        differential_paths = self.differential_paths_to_children[child_id]
        differential_points = np.zeros((1, self.n_markers))
        for m in range(self.n_markers):
            path = differential_paths[m]
            if path is not None:
                # If there is no differential path, we keep it at 0
                pseudo_time = rd.beta(a=alpha, b=beta, size=1)
                differential_points[0, m] = path(pseudo_time)
        return differential_points

    def generate_initial_expressions(self):
        """
        Generates
        :return:
        """
        initial_expressions = np.zeros((1, self.n_markers))
        expressed = self.model_for_expressed_markers.sample(1)[0][0]
        counter = 0
        for m in range(self.n_markers):
            if self.markers_pattern[0, m] != 0:
                initial_expressions[0, m] = expressed[counter]
                counter += 1
        return initial_expressions

    def generate_final_expressions(self, differentiate = True, alpha = 0.4, beta = 1):
        """
        Generates the final expression
        :param differentiate: A boolean variable indicating if differentiation is being considered
        :param alpha: alpha parameter in the Beta distribution
        :param beta: beta parameter in the Beta distribution
        :return: An array of the final expressions
        """
        final_expressions = self.generate_initial_expressions()
        if differentiate:
            children_ids = list(self.children_cell_types)
            child_id = rd.choice(children_ids, 1)[0]
            differential_points = self.differential_paths_to_children(child_id, alpha, beta)
            final_expressions += differential_points
        return final_expressions
