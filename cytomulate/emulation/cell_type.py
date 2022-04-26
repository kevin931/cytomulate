# Math computation
import numpy as np

# Statistical models
from sklearn.mixture import GaussianMixture

# Typing
from typing import Union, Optional, Any, List, Callable

# Superclass
from cell_type_general import GeneralCellType


class EmulationCellType(GeneralCellType):
    def __init__(self, label, cell_id, n_markers):
        super().__init__(label, cell_id, n_markers)
        self.observed_n = None
        self.observed_mean = None
        self.observed_covariance = None

    def fit(self, data, max_components,
            min_components, covariance_types):

        self.observed_n = data.shape[0]
        min_components = np.min([min_components, self.observed_n])
        max_components = np.min([max_components, self.observed_n])
        self.observed_mean = np.mean(data, axis=0)
        self.observed_covariance = np.cov(data, rowvar=False)

        self.markers = np.array(range(len(self.observed_mean)))

        smallest_bic = np.Inf
        current_bic = 0
        best_gm = None
        for n_components in range(min_components, max_components + 1):
            for cv_type in covariance_types:
                gm = GaussianMixture(n_components=n_components,
                                     covariance_type=cv_type).fit(data)
                current_bic = gm.bic(data)
                if current_bic < smallest_bic:
                    smallest_bic = current_bic
                    best_gm = gm

        self.model = best_gm
        if self.model.n_components == 1:
            self.model.weights_[0] = 1.

