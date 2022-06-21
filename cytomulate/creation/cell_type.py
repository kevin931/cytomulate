# Math computation
import numpy as np
from copy import deepcopy

# Statistical models
from scipy.stats import truncnorm
from scipy.stats import invwishart
from sklearn.mixture import GaussianMixture

# Typing
from typing import Union, Optional, Any, List, Callable
from cytomulate.creation.cell_graph import CreationCellGraph

# Superclass
from cytomulate.cell_type_general import GeneralCellType


class CreationCellType(GeneralCellType):
    def __init__(self,
                 label: Union[str, int],
                 cell_id: int,
                 n_markers: int) -> None:
        """Initialize the GeneralCellType object

        Parameters
        ----------
        label: str or int
            The label (name) for the cell type
        cell_id: int
            The id number assigned to the cell type
        n_markers: int
            Number of markers used in the experiment
        """
        super().__init__(label, cell_id, n_markers)
        self.highly_expressed_markers = None
        self.lowly_expressed_markers = None
        # Gating markers are the markers that will not change during differentiation
        self.gating_markers = None
        # The probability of a marker being highly expressed
        self.p = 0.4


    def generate_marker_expression_patterns(self,
                                            cell_types: dict,
                                            cell_graph: CreationCellGraph) -> None:
        """Generate marker patterns for the cell types

        Parameters
        ----------
        cell_types: dict
            A dictionary of CreationCellType objects
        cell_graph: CreationCellGraph
            A cell graph object
        """
        predecessor_cell_labels = list(cell_graph.predecessors(self.label))
        if len(predecessor_cell_labels) == 0:
            # This means this is a root
            indicator = np.random.binomial(n=1, p=self.p, size=self.n_markers)
            self.highly_expressed_markers = self.markers[np.where(indicator == 1)]
            self.lowly_expressed_markers = self.markers[np.where(indicator == 0)]
            self.gating_markers = np.random.choice(self.markers, size=2, replace=False)
        elif len(predecessor_cell_labels) == 1:
            # First inherit from the parent lists
            self.highly_expressed_markers = deepcopy(cell_types[predecessor_cell_labels[0]].highly_expressed_markers)
            self.lowly_expressed_markers = deepcopy(cell_types[predecessor_cell_labels[0]].lowly_expressed_markers)
            self.gating_markers = deepcopy(cell_types[predecessor_cell_labels[0]].gating_markers)

            # Then for the markers that can be changed
            # we randomly choose a subset to flip their signs
            changeable_markers = np.array(list(set(self.markers) - set(self.gating_markers)))
            indicator = np.random.binomial(n=1, p=self.p/2, size=len(changeable_markers))
            markers_to_change = changeable_markers[np.where(indicator == 1)]
            gating_counter = 0
            for m in markers_to_change:
                if m in self.highly_expressed_markers:
                    self.highly_expressed_markers = np.delete(self.highly_expressed_markers,
                                                              np.where(self.highly_expressed_markers == m))
                    self.lowly_expressed_markers = np.append(self.lowly_expressed_markers, m)
                else:
                    self.lowly_expressed_markers = np.delete(self.lowly_expressed_markers,
                                                              np.where(self.lowly_expressed_markers == m))
                    self.highly_expressed_markers = np.append(self.highly_expressed_markers, m)
                if gating_counter <= 1:
                    self.gating_markers = np.append(self.gating_markers, m)
                    gating_counter += 1
        else:
            raise ValueError('Cell graph is not a tree or collection of non-overlapping trees')


    def generate_marker_expressions(self,
                                    cell_types: dict,
                                    cell_graph: CreationCellGraph,
                                    high_expressions: np.ndarray,
                                    low_expressions: np.ndarray) -> None:
        """Generate the actual expressions

        Parameters
        ----------
        cell_types: dict
            A dictionary of CreationCellType objects
        cell_graph: CreationCellGraph
            A cell graph object
        high_expressions: np.ndarray
            An array of high expression levels
        low_expressions: np.ndarray
            An array of low expression levels
        """
        predecessor_cell_labels = list(cell_graph.predecessors(self.label))
        # We first generate the means of the cell expression menas
        if len(predecessor_cell_labels) == 0:
            # This means this is a root
            for m in self.highly_expressed_markers:
                self.cell_mean[m] = np.random.choice(high_expressions, size=1)[0]
            for m in self.lowly_expressed_markers:
                self.cell_mean[m] = np.random.choice(low_expressions, size=1)[0]
        elif len(predecessor_cell_labels) == 1:
            predecessor_gating_markers = cell_types[predecessor_cell_labels[0]].gating_markers
            for m in self.markers:
                if m in predecessor_gating_markers:
                    self.cell_mean[m] = cell_types[predecessor_cell_labels[0]].cell_mean[m]
                else:
                    if m in self.highly_expressed_markers:
                        self.cell_mean[m] = np.random.choice(high_expressions, size=1)[0]
                    if m in self.lowly_expressed_markers:
                        self.cell_mean[m] = np.random.choice(low_expressions, size=1)[0]
        else:
            raise ValueError('Cell graph is not a tree or collection of non-overlapping trees')

        # Now we use the means to generate the actual expressions
        for m in self.markers:
            self.cell_mean[m] = truncnorm.rvs(a=0, b=np.inf, loc=self.cell_mean[m], scale=0.01, size=1)[0]


    def generate_model(self,
                       n_components: int,
                       variance_mode: float = 0.01) -> None:
        """Randomly generate models for cell types

        Parameters
        ----------
        n_components: int
            Number of components in a GMM
        variance_mode: float
            The mode of the variance of the inverse wishart distribution
        """
        X = np.random.multivariate_normal(mean=self.cell_mean, cov=np.eye(self.n_markers),
                                          size=np.max([self.n_markers**2, n_components]))
        self.model = GaussianMixture(n_components=n_components).fit(X)
        for n in range(n_components):
            self.model.means_[n, :] = self.cell_mean

        base_correlation = invwishart.rvs(df=self.n_markers + 2, scale=np.eye(self.n_markers),
                                                 size=1).reshape((1, self.n_markers, self.n_markers))
        base_var = np.diag(base_correlation[0,:,:])
        base_var = np.diag(np.sqrt(1/base_var))
        base_correlation = base_var @ base_correlation @ base_var

        scaled_matrices = []
        variance_var = 0.01
        shape_par = (variance_mode**2)/variance_var + 2
        scale_par = variance_mode * (shape_par - 1)
        for n in range(n_components):
            scale_matrix = np.diag(1/np.random.gamma(shape_par, 1/scale_par, size=self.n_markers))
            scaled_matrices.append(scale_matrix @ base_correlation[0,:,:] @ scale_matrix)
        scaled_matrices = np.array(scaled_matrices)

        self.model.covariances_ = scaled_matrices

        # self.model.covariances_ = invwishart.rvs(df=self.n_markers + 2, scale=np.eye(self.n_markers) * variance_mode,
        #                                          size=n_components).reshape((n_components, self.n_markers, self.n_markers))
        self.model.weights_ = np.random.dirichlet(np.ones(n_components), size=1).reshape(-1)
