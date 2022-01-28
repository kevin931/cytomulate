
import numpy as np
from sklearn.cluster import KMeans
from sklearn.mixture import GaussianMixture


class CellType:
    def __init__(self, label, id):
        self.label = label
        self.id = id
        self.highly_expressed_markers = None
        self.lowly_expressed_markers = None
        self.model_for_highly_expressed_markers = None
        self.model_for_lowly_expressed_markers = None
        self.observed_mean = None
        self.observed_covariance = None

    def fit(self, data, max_components, min_components, covariance_types):
        col_median = np.median(data, axis=0).reshape(-1,1)
        # We will use K-means with 2 groups to
        # group markers into two groups
        # expressed and unexpressed
        kmeans = KMeans(n_clusters=2).fit(col_median)

        group_0_mean = np.mean(col_median[kmeans.labels_ == 0])
        group_1_mean = np.mean(col_median[kmeans.labels_ == 1])
        if group_1_mean > group_0_mean:
            highly_expressed_group = 1
            lowly_expressed_group = 0
        else:
            highly_expressed_group = 0
            lowly_expressed_group = 1

        self.highly_expressed_markers = np.where(kmeans.labels_ == highly_expressed_group)[0]
        self.lowly_expressed_markers = np.where(kmeans.labels_ == lowly_expressed_group)[0]

        # For unexpressed markers
        # we fit gaussian mixture model with 2 components to each marginal
        # as it seems reasonable to assume unexpressed markers are independent
        # one with background noise
        # one with lowly expressed protein
        # The main purpose is to speed up the algorithm
        # since the majority of the markers are un/lowly expressed
        lowly_expressed_data = data[:, self.lowly_expressed_markers]
        highly_expressed_data = data[:, self.highly_expressed_markers]

        self.model_for_highly_expressed_markers = {}
        self.model_for_lowly_expressed_markers = {}

        counter = 0
        for m in self.lowly_expressed_markers:
            self.model_for_lowly_expressed_markers[m] = GaussianMixture(n_components=2).fit(lowly_expressed_data[:, counter].reshape(-1, 1))
            counter += 1

        # Since we expect only a few markers to be highly expressed
        # we can probably afford fitting the entire data
        # with multivariate GMM as well as some model selection
        smallest_bic = np.Inf
        current_bic = 0
        best_gm = None
        for n_components in range(min_components, max_components + 1):
            for cv_type in covariance_types:
                gm = GaussianMixture(n_components=n_components,
                                     covariance_type=cv_type).fit(highly_expressed_data)
                current_bic = gm.bic(highly_expressed_data)
                if current_bic < smallest_bic:
                    smallest_bic = current_bic
                    best_gm = gm

        self.model_for_highly_expressed_markers["all"] = best_gm

    def adjust_models(self):
        pass

    def sample_cell(self, n_samples):
        n_markers = len(self.lowly_expressed_markers) + len(self.highly_expressed_markers)
        result = np.zeros((n_samples, n_markers))
        result[:, self.highly_expressed_markers], _ = self.model_for_highly_expressed_markers["all"].sample(n_samples)
        for m in self.lowly_expressed_markers:
            result[:, [m]], _ = self.model_for_lowly_expressed_markers[m].sample(n_samples)
        return result
