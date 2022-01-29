
import numpy as np
from collections import Counter

from cell_type import CellType


class CytofData:
    def __init__(self):

        self.n_markers = None
        self.background_noise_variance = None
        self.bead_label = None
        self.observed_cell_abundances = {}
        self.cell_type_labels_to_ids = {}
        self.cell_type_ids_to_labels = {}
        self.cell_types = {}

    def initialize_cell_types(self, expression_matrix,
                              labels,
                              bead_label=None,
                              bead_channels=None,
                              max_components=9,
                              min_components=1,
                              covariance_types=("full", "tied", "diag", "spherical")):
        self.n_markers = np.shape(expression_matrix)[1]

        unique_labels = np.unique(labels)

        abundances = Counter(labels)

        self.background_noise_variance = np.Inf
        cell_id = 0
        for c_type in unique_labels:
            self.observed_cell_abundances[c_type] = abundances[c_type]/len(labels)

            self.cell_type_labels_to_ids[c_type] = cell_id
            self.cell_type_ids_to_labels[cell_id] = c_type

            self.cell_types[c_type] = CellType(label=c_type, cell_id=cell_id)

            ind = np.where(labels == c_type)[0]
            D = expression_matrix[ind, :]
            self.cell_types[c_type].fit(data=D,
                                        max_components=max_components,
                                        min_components=min_components,
                                        covariance_types=covariance_types,
                                        is_bead=(c_type == bead_label),
                                        bead_channels=bead_channels)

            if (bead_label is None) or (c_type == bead_label):
                for m in self.cell_types[c_type].lowly_expressed_markers:
                    est_v = np.min(self.cell_types[c_type].model_for_lowly_expressed_markers[m].covariances_)
                    if est_v <= self.background_noise_variance:
                        self.background_noise_variance = est_v

            cell_id += 1

    def adjust_cell_types(self):
        for c_type in self.cell_types:
            self.cell_types[c_type].adjust_models(self.background_noise_variance)

    def sample(self, n_samples,
               cell_abundances = None,
               ordered_by = "ids"):
        cell_numbers = {}
        if cell_abundances is None:
            cell_numbers = self.observed_cell_abundances
            ordered_by = "labels"
        else:
            if type(cell_abundances) is dict:
                if ordered_by == "ids":
                    if cell_abundances.keys() != self.cell_type_ids_to_labels.keys():
                        raise ValueError('Keys do not match')
                    for cell_id in range(len(self.cell_types)):
                        cell_numbers[self.cell_type_ids_to_labels[cell_id]] = cell_abundances[cell_id]
                elif ordered_by == "labels":
                    if cell_abundances.keys() != self.cell_type_labels_to_ids.keys():
                        raise ValueError('Keys do not match')
                    cell_numbers = cell_abundances
                else:
                    raise ValueError('Unknown order type')
            else:
                if ordered_by == "ids":
                    for cell_id in range(len(self.cell_types)):
                        cell_numbers[self.cell_type_ids_to_labels[cell_id]] = cell_abundances[cell_id]
                elif set(ordered_by) == self.cell_type_labels_to_ids.keys():
                    for i in range(len(self.cell_types)):
                        cell_numbers[ordered_by[i]] = cell_abundances[i]
                else:
                    raise ValueError('Unknown order type')

        cell_probabilities = np.zeros(len(self.cell_types))
        for cell_id in range(len(self.cell_types)):
            cell_probabilities[cell_id] = cell_numbers[self.cell_type_ids_to_labels[cell_id]]
        cell_probabilities /= np.sum(cell_probabilities)

        expression_matrix = np.zeros((n_samples, self.n_markers))

        n_per_cell_type = np.random.multinomial(n_samples, cell_probabilities)
        ids = np.repeat(range(len(n_per_cell_type)), n_per_cell_type)

        start_n = 0
        end_n = 0
        for cell_id in range(len(n_per_cell_type)):
            c_type = self.cell_type_ids_to_labels[cell_id]
            n = n_per_cell_type[cell_id]
            if n == 0:
                continue
            end_n += n
            expression_matrix[start_n : end_n, :] = self.cell_types[c_type].sample_cell(n)
            start_n += n

        return expression_matrix, ids
