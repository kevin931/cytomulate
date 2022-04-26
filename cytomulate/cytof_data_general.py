
from copy import deepcopy

import numpy as np

from utilities import trajectories
from cell_graph_general import GeneralCellGraph


class GeneralCytofData:
    def __init__(self, n_batches=1, background_noise_model=None):

        self.n_markers = None

        self.background_noise_model = background_noise_model

        self.n_batches = n_batches
        self.overall_batch_effects = {}
        self.local_batch_effects = {}
        self.temporal_effects = {}
        self.cell_abundances = {}

        self.cell_type_labels_to_ids = {}
        self.cell_type_ids_to_labels = {}
        self.cell_types = {}

        self.cell_graph = GeneralCellGraph()

    def generate_cell_abundances(self, is_random=True):
        if is_random:
            for b in range(self.n_batches):
                self.cell_abundances[b] = {}
                abundance = np.random.dirichlet(np.ones(len(self.cell_types)), size=1).reshape(-1)
                counter = 0
                for c_type in self.cell_types:
                    self.cell_abundances[b][c_type] = abundance[counter]
                    counter += 1
        else:
            for b in range(self.n_batches):
                self.cell_abundances[b] = {}
                abundance = np.ones(len(self.cell_types))/len(self.cell_types)
                counter = 0
                for c_type in self.cell_types:
                    self.cell_abundances[b][c_type] = abundance[counter]
                    counter += 1

    def generate_overall_batch_effects(self, variance=0.001):
        if self.n_batches == 1:
            self.overall_batch_effects[0] = 0
        else:
            batch_effects = np.zeros(self.n_batches)
            batch_effects[:(self.n_batches - 1)] = np.random.normal(loc=0, scale=np.sqrt(variance),
                                                                    size=self.n_batches - 1)
            batch_effects[self.n_batches - 1] = -np.sum(batch_effects)
            for b in range(self.n_batches):
                self.overall_batch_effects[b] = batch_effects[b]

    def generate_local_batch_effects(self, variance=0.001):
        if self.n_batches == 1:
            self.local_batch_effects[0] = {}
            for c_type in self.cell_type_labels_to_ids:
                self.local_batch_effects[0][c_type] = np.zeros(self.n_markers)
        else:
            for b in range(self.n_batches):
                self.local_batch_effects[b] = {}
                temp = np.zeros((len(self.cell_types), self.n_markers))
                temp[np.ix_(range(len(self.cell_types) - 1), range(self.n_markers - 1))] = \
                    np.random.normal(loc=0, scale=np.sqrt(variance),
                                     size=(len(self.cell_types) - 1) * (self.n_markers - 1)).reshape(
                        (len(self.cell_types) - 1, self.n_markers - 1))
                temp[:, self.n_markers - 1] = -np.sum(temp, axis=1)
                temp[len(self.cell_types) - 1, :] = -np.sum(temp, axis=0)
                counter = 0
                for c_type in self.cell_type_labels_to_ids:
                    self.local_batch_effects[b][c_type] = temp[counter, :].reshape(-1)
                    counter += 1

    def generate_temporal_effects(self, variance=None, coefficients=None, x=None, y=None, **kwargs):
        for b in range(self.n_batches):
            if variance is not None:
                if coefficients is None:
                    self.temporal_effects[b] = trajectories(end_values=np.random.normal(0, np.sqrt(variance), 1),
                                                            **kwargs)
                else:
                    self.temporal_effects[b] = trajectories(end_values=np.random.normal(0, np.sqrt(variance), 1),
                                                            coefficients=coefficients[b],
                                                            **kwargs)
            else:
                self.temporal_effects[b] = trajectories(x=x[b], y=y[b], **kwargs)

    def sample_one_batch(self, n_samples,
                         cell_abundances=None,
                         batch=0,
                         clip=True):
        if cell_abundances is None:
            cell_abundances = self.cell_abundances[batch]

        # We record the order of the cell types
        cell_type_order = list(cell_abundances.keys())

        n_per_cell_type = np.zeros(len(cell_abundances), dtype=int)

        if np.sum(list(cell_abundances.values())) == n_samples:
            counter = 0
            for c_type in cell_type_order:
                n_per_cell_type[counter] = cell_abundances[c_type]
                counter += 1
        elif np.sum(list(cell_abundances.values())) <= 1.5:
            cell_probabilities = np.zeros(len(cell_abundances))
            counter = 0
            for c_type in cell_type_order:
                cell_probabilities[counter] = cell_abundances[c_type]
                counter += 1
            cell_probabilities /= np.sum(cell_probabilities)
            n_per_cell_type = np.random.multinomial(n_samples, cell_probabilities)
        else:
            raise ValueError('Unknown cell abundances type')

        expression_matrix = np.zeros((n_samples, self.n_markers))
        expressed_index_matrix = np.zeros((n_samples, self.n_markers))
        pseudo_time = np.zeros((n_samples, self.n_markers))
        children_cell_labels = ["None"] * n_samples

        labels = [item for item, count in zip(cell_type_order, n_per_cell_type) for i in range(count)]

        Psi_b = 0
        if batch in self.overall_batch_effects.keys():
            Psi_b = self.overall_batch_effects[batch]

        start_n = 0
        end_n = 0
        counter = 0
        for c_type in cell_type_order:
            n = n_per_cell_type[counter]
            counter += 1
            if n == 0:
                continue
            end_n += n
            X, expressed_index = self.cell_types[c_type].sample_cell(n, clip)
            Psi_bp = 0
            if batch in self.local_batch_effects.keys():
                Psi_bp = self.local_batch_effects[batch][c_type]
            G, T, children_labels = self.cell_graph.sample_graph(n, c_type)
            expression_matrix[start_n: end_n, :] = X + expressed_index * (G + Psi_b + Psi_bp)
            expressed_index_matrix[start_n: end_n, :] = expressed_index
            pseudo_time[start_n:end_n, :] = T
            children_cell_labels[start_n:end_n] = children_labels
            start_n += n

        # To add temporal effects, we first shuffle the arrays
        indices = np.random.permutation(n_samples)
        expression_matrix = expression_matrix[indices, :]
        expressed_index_matrix = expressed_index_matrix[indices, :]
        labels = [labels[i] for i in indices]
        pseudo_time = pseudo_time[indices, :]
        children_cell_labels = [children_cell_labels[i] for i in indices]

        # Now we add temporal effects]
        time_points = np.linspace(0, 1, n_samples)

        if batch in self.temporal_effects.keys():
            temporal_effects = self.temporal_effects[batch][0](time_points)
            expression_matrix += expressed_index_matrix * temporal_effects[:, np.newaxis]
        if clip:
            expression_matrix = np.clip(expression_matrix, a_min=0, a_max=None)

        E = 0
        if self.background_noise_model is not None:
            E = self.background_noise_model(size=(n_samples, self.n_markers))

        expression_matrix += E

        return expression_matrix, np.array(labels), pseudo_time, np.array(children_cell_labels)

    def sample(self, n_samples,
               cell_abundances=None,
               clip=True):
        if cell_abundances is None:
            self.generate_cell_abundances()
            cell_abundances = self.cell_abundances

        if not np.any([isinstance(i, dict) for i in cell_abundances.values()]):
            cell_abundances_copy = deepcopy(cell_abundances)
            cell_abundances = {}
            for b in range(self.n_batches):
                cell_abundances[b] = cell_abundances_copy

        if isinstance(n_samples, int):
            n_samples = np.repeat(n_samples, self.n_batches)

        expression_matrices = {}
        labels = {}
        pseudo_time = {}
        children_cell_labels = {}
        for b in range(self.n_batches):
            expression_matrices[b], labels[b], pseudo_time[b], children_cell_labels[b] = self.sample_one_batch(
                n_samples[b],
                cell_abundances[b],
                b,
                clip)

        return expression_matrices, labels, pseudo_time, children_cell_labels

