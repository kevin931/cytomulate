
import numpy as np
from collections import Counter

from cell_type import CellType


class CytofData:
    def __init__(self):

        self.background_noise_variance = None
        self.cell_types = {}
        self.cell_type_labels_and_ids = {}
        self.observed_cell_abundances = None
        self.n_markers = None

    def initialize_cell_types(self, expression_matrix,
                              labels,
                              max_components=9,
                              min_components=1,
                              covariance_types=("full", "tied", "diag", "spherical")):
        self.n_markers = np.shape(expression_matrix)[1]
        self.cell_type_labels_and_ids = dict.fromkeys(np.unique(labels))
        abundances = Counter(labels)
        self.background_noise_variance = np.Inf
        self.observed_cell_abundances = np.zeros(len(self.cell_type_labels_and_ids))
        cell_id = 0
        for c_type in self.cell_type_labels_and_ids:
            self.observed_cell_abundances[cell_id] = abundances[c_type]/len(labels)
            self.cell_type_labels_and_ids[c_type] = cell_id
            self.cell_types[c_type] = CellType(label=c_type, id=cell_id)
            self.cell_types[c_type].observed_n = abundances[c_type]
            ind = np.where(labels == c_type)[0]
            D = expression_matrix[ind, :]
            self.cell_types[c_type].observed_mean = np.mean(D, axis=0)
            self.cell_types[c_type].observed_covariance = np.cov(D, rowvar=False)
            self.cell_types[c_type].fit(data=D,
                                        max_components=max_components,
                                        min_components=min_components,
                                        covariance_types=covariance_types)

            for m in self.cell_types[c_type].lowly_expressed_markers:
                est_v = np.min(self.cell_types[c_type].model_for_lowly_expressed_markers[m].covariances_)
                if est_v <= self.background_noise_variance:
                    self.background_noise_variance = est_v

            cell_id += 1


    def sample(self, n_samples):
        pass