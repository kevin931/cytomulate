
import numpy as np
from sklearn.cluster import KMeans
from sklearn.mixture import GaussianMixture
from scipy.stats import t


class CellType:
    def __init__(self, label, cell_id):
        self.label = label
        self.id = cell_id
        self.highly_expressed_markers = None
        self.lowly_expressed_markers = None
        self.model_for_highly_expressed_markers = None
        self.model_for_lowly_expressed_markers = None
        self.observed_n = None
        self.observed_mean = None
        self.observed_covariance = None

    def fit(self, data, max_components,
            min_components, covariance_types,
            is_bead=False, bead_channels=None):

        self.observed_n = data.shape[0]
        self.observed_mean = np.mean(data, axis=0)
        self.observed_covariance = np.cov(data, rowvar=False)

        if is_bead:
            self.highly_expressed_markers = bead_channels
            self.lowly_expressed_markers = list(set(range(len(self.observed_mean))) - set(self.highly_expressed_markers))
        else:
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


        if is_bead:
            n_components = 1
        else:
            n_components = 1

        counter = 0
        for m in self.lowly_expressed_markers:
            self.model_for_lowly_expressed_markers[m] = GaussianMixture(n_components=1).fit(lowly_expressed_data[:, counter].reshape(-1, 1))
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

    def adjust_models(self, background_noise_variance):
        # We deal with expressed markers first
        n_components = self.model_for_highly_expressed_markers["all"].n_components
        if self.model_for_highly_expressed_markers["all"].covariance_type == "spherical":
            self.model_for_highly_expressed_markers["all"].covariances_ -= background_noise_variance
        elif self.model_for_highly_expressed_markers["all"].covariance_type == "tied":
            self.model_for_highly_expressed_markers["all"].covariances_ -= background_noise_variance * np.eye(len(self.highly_expressed_markers))
        elif self.model_for_highly_expressed_markers["all"].covariance_type == "diag":
            for c in range(n_components):
                self.model_for_highly_expressed_markers["all"].covariances_[c, :] -= background_noise_variance
        elif self.model_for_highly_expressed_markers["all"].covariance_type == "full":
            for c in range(n_components):
                self.model_for_highly_expressed_markers["all"].covariances_[c, :,:] -= background_noise_variance * np.eye(len(self.highly_expressed_markers))
        else:
            raise ValueError('Unknown covariance type')

        # Then we deal with lowly/ unexpressed markers
        for m in self.lowly_expressed_markers:
            n_components = self.model_for_lowly_expressed_markers[m].n_components
            for c in range(n_components):
                self.model_for_lowly_expressed_markers[m].covariances_[c, :, :] -= background_noise_variance
                if self.model_for_lowly_expressed_markers[m].covariances_[c, :, :] <= 0:
                    self.model_for_lowly_expressed_markers[m].means_[c] = 0
                    self.model_for_lowly_expressed_markers[m].covariances_[c, :, :] = 0
                else:
                    t_stat = (self.model_for_lowly_expressed_markers[m].means_[c]) / \
                             np.sqrt(self.model_for_lowly_expressed_markers[m].covariances_[c, :, :]/self.observed_n)
                    n095quantile = t.ppf(0.95, df=self.observed_n - 1)
                    if t_stat <= n095quantile:
                        self.model_for_lowly_expressed_markers[m].means_[c] = 0
                        self.model_for_lowly_expressed_markers[m].covariances_[c, :, :] = 0

    def sample_cell(self, n_samples):
        n_markers = len(self.observed_mean)
        result = np.zeros((n_samples, n_markers))
        result[:, self.highly_expressed_markers], _ = self.model_for_highly_expressed_markers["all"].sample(n_samples)
        for m in self.lowly_expressed_markers:
            result[:, [m]], _ = self.model_for_lowly_expressed_markers[m].sample(n_samples)
        return result
