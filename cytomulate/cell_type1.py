# Math computation
import numpy as np
from numpy import random as rd

# List manipulation
from copy import deepcopy

# Gaussian Mixtures
from sklearn.mixture import GaussianMixture

from typing import Union, Optional


class CellType:
    """The Cell Type Class

    This class the basis for all cells in our simulation. It stores all the relavant information for
    the cell and computes the expressions for the cell.
    """
    def __init__(self, id: int, name: Union[str, int], n_markers: int):
        """The CellType Class constructor

        This constructor sets up the bare-minimum of a cell.

        :param id: The ID of the cell
        :type id: int
        :param name: The name of the cell
        :type name: Union[str, int]
        :param n_markers: The number of markers
        :type n_markers: int
        """
        self.id = id
        self.name = name
        self.n_markers = n_markers

        self.parent = None
        self.children = []
        self.marker_pattern = np.zeros(n_markers)
        self.expression_level = np.zeros(n_markers)
        self.variance_level = np.zeros(n_markers)

        # parent_cell_type will be a dictionary of size 1
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
        self.expressed_markers = None
        self.unexpressed_markers = None
        self.true_zero_indices = {}
        # We will use Gaussian Mixture to model or learn
        # the expression of this CellType
        # Consequently, there will be no need for
        # us to store expression levels and
        # variance levels separately
        # We expect the class of this field
        # to be sklearn.mixture._gaussian_mixture.GaussianMixture
        self.model_for_expressed_markers = None
        self.model_for_unexpressed_markers = None
        # Gating markers will be a set of indices
        self.gating_markers = set()

        self.observed_mean = np.zeros((1, n_markers))
        self.observed_cov = np.zeros((n_markers, n_markers))
        self.overall_mean = np.zeros((1, n_markers))
        self.overall_cov = np.zeros((n_markers, n_markers))
        self.background_noise_level = None

    def sample_cell(self, n_samples = 1):
        result = np.zeros((n_samples, self.n_markers))
        result[:, self.expressed_markers], _ = self.model_for_expressed_markers["all"].sample(n_samples)
        for m in self.unexpressed_markers:
            result[:, [m]], _ = self.model_for_unexpressed_markers[m].sample(n_samples)
        return result

    def inherit_markers_pattern(self, mutation_probability = 0.2, \
                                n_additional_gating_markers = 2):
        """Inherit markers pattern from parent

        If the cell is not a root, this method inherits the marker pattern
        from the parent cell.

        :param mutation_probability: The mutation probability of markers by flipping the marker pattern, defaults to 0.2
        :type mutation_probability: float, optional
        :raises Exception: This is the root cell without a parent
        :raises ValueError: Mutation probability is not between 0 and 1
        :raises ValueError: Marker number exceeds the number of remaining markers (not counting gating markers)
        :raises ValueError: Marker numbers are not positive
        """
        
        if len(self.parent_cell_type) == 0:
            raise Exception("This is the root.")
        if mutation_probability > 1 or mutation_probability < 0:
            raise ValueError("Mutation probability must be between 0 and 1.")
        if n_additional_gating_markers > self.n_markers - len(self.gating_markers):
            raise ValueError("Maker number exceeds the number of remaining markers.")
        if n_additional_gating_markers < 0:
            raise ValueError("Marker number has to be positive.")
        
        # Get the parent
        parent = list(self.parent_cell_type.values())[0]
        # simply copy the information needed
        self.markers_pattern = deepcopy(parent.markers_pattern)
        self.gating_markers = deepcopy(parent.gating_markers)

        # fickle markers are indices of the markers that can change
        # expression levels or even from 1 to 0
        fickle_markers = set(range(self.n_markers)) - self.gating_markers
        # flip markers are fickle markers that will flip from 0 to 1 and vice versa
        flip_markers = set(rd.choice(list(fickle_markers), int(mutation_probability * len(fickle_markers)), False))
        for m in flip_markers:
            self.markers_pattern[0, m] = -1 * self.markers_pattern[0, m] + 1
        new_gating_markers = set(rd.choice(list(fickle_markers), n_additional_gating_markers, False))
        self.gating_markers = self.gating_markers.union(new_gating_markers)


    def inherit_model(self):
        """Inherit the model from the parent cell

        If the cell is not the root, i.e. there is a parent cell, this model inherits
        the Gaussian Mixture Model from the parent cell.

        :raises Exception: This is the root cell without parent.
        """
        if len(self.parent_cell_type) == 0:
            # if this is the root
            raise Exception("This is the root.")
        # Get the parent
        parent = list(self.parent_cell_type.values())[0]
        self.model_for_expressed_markers = deepcopy(parent.model_for_expressed_markers)


    def fit_model(self, dat: Optional[np.ndarray] = None, max_n_components: int = 9, min_n_components: int = 1):
        """Fit the Gaussian Mixture Model of the cell

        For each cell, a Gaussian Mixture Model is used to act as the base for initial expression levels. 
        This method fits a Gaussian Mixture Model with a range of components and the best model based on
        BIC.

        :param dat: The expression matrix for data-based simulation, defaults to None
        :type dat: Optional[np.ndarray], optional
        :param max_n_components: Maximum number of components, defaults to 9
        :type max_n_components: int, optional
        :param min_n_components: Minimum number of components, defaults to 1
        :type min_n_components: int, optional
        """
        best_model = GaussianMixture(n_components = max_n_components, covariance_type = "diag").fit(dat)
        best_bic = best_model.bic(dat)
        for n_comp in range(min_n_components, max_n_components):
            temp_model = GaussianMixture(n_components = n_comp, covariance_type = "diag").fit(dat)
            temp_bic = temp_model.bic(dat)
            if temp_bic < best_bic:
                best_model = temp_model
                best_bic = temp_bic
        self.model_for_expressed_markers = best_model

        self.overall_mean = self.model_for_expressed_markers.weights_.dot(self.model_for_expressed_markers.means_)
        self.overall_var = None


    def generate_from_paths(self, child_id: Union[str, int], alpha: float = 0.4, beta: float = 1) -> np.ndarray:
        """Generates an array of differential points
        
        This array of differential points are used to simulation the differentiation path between
        cell populations. These values are to be added to the original expressions to simulate the
        final expressions.
        
        :param child_id: ID of the child to differentiate to
        :type child_id:
        :param alpha: alpha parameter in the Beta distribution, defaults to 0.4
        :type alpha: float
        :param beta: beta parameter in the Beta distribution, defaults to to 1
        :type float: float
        :return: an array of differential points
        :rtype: np.ndarray
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


    def generate_initial_expressions(self) -> np.ndarray:
        """ Generates the initial expressions of the cell
        
        This method samples from the Gaussian Mixture Model stored in
        ``model_for_expressed_markers``. The resulting expressions are
        a 1-dimensional array with a value for each marker.
        
        :return: The initial expressions of the cell
        :rtype: np.ndarray
        """
        assert self.model_for_expressed_markers is not None
        initial_expressions = np.zeros((1, self.n_markers))
        initial_expressions = self.model_for_expressed_markers.sample(1)[0][0]
        # counter = 0
        # for m in range(self.n_markers):
        #     if self.markers_pattern[0, m] != 0:
        #         initial_expressions[0, m] = expressed[counter]
        #         counter += 1
        return initial_expressions


    def generate_final_expressions(self, differentiate: bool = True, alpha: float = 0.4, beta: float = 1) -> np.ndarray:
        """ Generates the final expression of the cell
        
        The final expression of the cell takes differentiation into consideration. The values are generated
        using the initial expressions as a starting point and add the differentiation path, which is
        in turn generated by a Brownian Bridge or a linear function connecting the cell and its child node.
        
        :param differentiate: A boolean variable indicating if differentiation is being considered,
            defaults to True
        :type differentiate: bool
        :param alpha: alpha parameter in the Beta distribution, defaults to 0.4
        :type alpha: float
        :param beta: beta parameter in the Beta distribution, defaults to 1
        :type beta: float
        :return: An array of the final expressions
        :rtype: np.ndarray
        """
        final_expressions = self.generate_initial_expressions()
        if differentiate:
            children_ids = list(self.children_cell_types)
            child_id = rd.choice(children_ids, 1)[0]
            differential_points = self.generate_from_paths(child_id, alpha, beta)
            final_expressions += differential_points.reshape(self.n_markers, )
        return final_expressions
