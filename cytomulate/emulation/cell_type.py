# Math computation
import numpy as np

# Statistical models
from sklearn.mixture import GaussianMixture


class CellType:
    def __init__(self, label, cell_id):
        self.label = label
        self.id = cell_id
        self.markers = None
        self.model = None
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

    def sample_cell(self, n_samples, clip):
        n_markers = len(self.observed_mean)
        X = np.zeros((n_samples, n_markers))
        X[:, self.markers], _ = self.model.sample(n_samples)
        expressed_index = (X > 0)
        if clip:
            X = np.clip(X, a_min=0, a_max=None)
        return X, expressed_index
