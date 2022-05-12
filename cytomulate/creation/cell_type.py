# Math computation
import numpy as np
from copy import deepcopy

# Statistical models
from scipy.stats import truncnorm
from scipy.stats import invwishart
from sklearn.mixture import GaussianMixture

# Typing
from typing import Union, Optional, Any, List, Callable
from creation.cell_graph import CreationCellGraph

# Superclass
from cell_type_general import GeneralCellType


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

    def generate_marker_expressions(self, cell_types, cell_graph, high_expressions, low_expressions):
        predecessor_cell_labels = list(cell_graph.predecessors(self.label))
        if len(predecessor_cell_labels) == 0:
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

        for m in self.markers:
            self.cell_mean[m] = truncnorm.rvs(a=0, b=np.inf, loc=self.cell_mean[m], scale=0.01, size=1)[0]

    def generate_model(self, n_components):
        X = np.random.multivariate_normal(mean=self.cell_mean, cov=np.eye(self.n_markers),
                                          size=np.max([self.n_markers**2, n_components]))
        self.model = GaussianMixture(n_components=n_components).fit(X)
        for n in range(n_components):
            self.model.means_[n, :] = self.cell_mean
        self.model.covariances_ = invwishart.rvs(df=self.n_markers + 2, scale=np.eye(self.n_markers)/100,
                                                 size=n_components)
        self.model.weights_ = np.random.dirichlet(np.ones(n_components), size=1).reshape(-1)
