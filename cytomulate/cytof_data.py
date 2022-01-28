
import numpy as np
from collections import Counter

from cell_type import CellType

class CytofData:
    def __init__(self):

        self.background_noise_variance = None
        self.cell_types = {}
        self.cell_types_data = {}
        self.cell_type_labels_and_ids = {}
        self.observed_cell_abundances = None
        self.n_markers = None

    def initialize_cell_types(self, expression_matrix, labels):
        self.n_markers = np.shape(expression_matrix)[1]
        self.cell_type_labels_and_ids = dict.fromkeys(np.unique(labels))
        abundances = Counter(labels)

        self.observed_cell_abundances = np.zeros(len(self.cell_type_labels_and_ids))
        id = 0
        for c_type in self.cell_type_labels_and_ids:
            self.observed_cell_abundances[id] = abundances[c_type]/len(labels)
            self.cell_type_labels_and_ids[c_type] = id
            self.cell_types[c_type] = CellType(label=c_type, id=id)

            ind = np.where(labels == c_type)[0]
            self.cell_types_data[c_type] = expression_matrix[ind, :]

            self.cell_types[c_type].observed_mean = np.mean(self.cell_types_data[c_type], axis=0)
            self.cell_types[c_type].observed_covariance = np.cov(self.cell_types_data[c_type], rowvar=False)

            id += 1

    def fit_cell_type_model(self, max_components=9,
                            min_components=1,
                            covariance_types=("full", "tied", "diag", "spherical")):
        for c_type in self.cell_types:
            self.cell_types[c_type].fit(data=self.cell_types_data[c_type],
                                        max_components=max_components,
                                        min_components=min_components,
                                        covariance_types=covariance_types)


    def estimate_background_noise_variance(self):
        self.background_noise_variance = np.Inf
        for c_type in self.cell_types:
            pass

    def sample(self, n_samples):
        pass