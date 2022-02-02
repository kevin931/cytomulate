
import numpy as np
from collections import Counter
from utilities import smooth_brownian_bridge
from cell_type import CellType
from cell_network import CellNetwork


class CytofData:
    def __init__(self, n_batches = 1):

        self.n_markers = None

        self.background_noise_variance = None
        self.bead_label = None

        self.n_batches = n_batches
        self.overall_batch_effects = {}
        self.local_batch_effects = {}
        self.temporal_effects = {}

        self.observed_cell_abundances = {}
        self.cell_type_labels_to_ids = {}
        self.cell_type_ids_to_labels = {}
        self.cell_types = {}

        self.cell_network = None

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

    def generate_cell_network(self, network_topology = "forest"):
        self.cell_network = CellNetwork()
        self.cell_network.initialize_network(self.cell_types, bead_label=self.bead_label)
        self.cell_network.prune_network(network_topology)

    def generate_cell_network_trajectories(self, cell_types, N = 5,
                              function_type = "linear", lb = 0, ub = 1):
        self.cell_network.generate_trajectories(self.cell_types, N, function_type, lb, ub)

    def generate_overall_batch_effects(self, variance=0.001):
        if self.n_batches == 1:
            self.overall_batch_effects[0] = 0
        else:
            batch_effects = np.zeros(self.n_batches)
            batch_effects[:(self.n_batches-1)] = np.random.normal(loc=0, scale=np.sqrt(variance),
                                             size=self.n_batches - 1)
            batch_effects[self.n_batches-1] = -np.sum(batch_effects)
            for b in range(self.n_batches):
                self.overall_batch_effects[b] = batch_effects[b]

    def generate_local_batch_effects(self, variance = 0.001):
        if self.n_batches == 1:
            self.local_batch_effects[0] = np.zeros((len(self.cell_types), self.n_markers))
        else:
            for b in range(self.n_batches):
                self.local_batch_effects[b] = np.zeros((len(self.cell_types), self.n_markers))
                self.local_batch_effects[b][np.ix_(range(len(self.cell_types)-1), range(self.n_markers-1))] = \
                    np.random.normal(loc=0, scale=np.sqrt(variance), size=(len(self.cell_types)-1)*(self.n_markers-1)).reshape((len(self.cell_types)-1, self.n_markers-1))
                self.local_batch_effects[b][:, self.n_markers-1] = -np.sum(self.local_batch_effects[b], axis=1)
                self.local_batch_effects[b][len(self.cell_types)-1, :] = -np.sum(self.local_batch_effects[b], axis=0)

    def generate_temporal_effects(self, variance= 0.001, N = 5,
                              function_type = "linear", lb = 0, ub = 1):
        for b in range(self.n_batches):
            self.temporal_effects[b] = smooth_brownian_bridge(np.random.normal(0,np.sqrt(variance),1),
                                                              N, function_type, lb, ub)

    def sample_one_batch(self, n_samples,
                         cell_abundances = None):
        if cell_abundances is None:
            cell_abundances = self.observed_cell_abundances
        # We record the order of the cell types
        cell_type_order = list(cell_abundances.keys())

        cell_probabilities = np.zeros(len(cell_abundances))
        counter = 0
        for c_type in cell_abundances:
            cell_probabilities[counter] = cell_abundances[c_type]
            counter += 1
        cell_probabilities /= np.sum(cell_probabilities)

        expression_matrix = np.zeros((n_samples, self.n_markers))
        pseudo_time = np.zeros((n_samples, self.n_markers))

        n_per_cell_type = np.random.multinomial(n_samples, cell_probabilities)
        labels = np.repeat(cell_type_order, n_per_cell_type)

        start_n = 0
        end_n = 0
        for cell_id in range(len(n_per_cell_type)):
            c_type = self.cell_type_ids_to_labels[cell_id]
            n = n_per_cell_type[cell_id]
            if n == 0:
                continue
            end_n += n
            X = self.cell_types[c_type].sample_cell(n)
            G, T, children_labels = self.cell_network.sample_network(n, c_type)
            E = np.random.normal(loc=0, scale=np.sqrt(self.background_noise_variance), size=(n, self.n_markers))
            expression_matrix[start_n : end_n, :] = X + G + E
            pseudo_time[start_n:end_n, :] = T
            start_n += n

        return expression_matrix, labels

    def sample(self, n_samples,
                     cell_abundances=None):
        expression_matrices = {}
        for b in range(self.n_batches):
            expression_matrices[b] = self.sample_one_batch(n_samples[b],
                                                           cell_abundances[b])

        return expression_matrices
