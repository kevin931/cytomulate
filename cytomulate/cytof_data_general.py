# Math computation
import numpy as np

# List/Dict manipulation
from copy import deepcopy

# Trajectory functions
from cytomulate.utilities import trajectories
from cytomulate.cell_graph_general import GeneralCellGraph

# Typing
from typing import Union, Optional, Tuple, Callable


class GeneralCytofData:
    def __init__(self,
                 n_batches: int = 1,
                 background_noise_model: Optional[Callable] = None) -> None:
        """Initialize the GeneralCytofData object

        Parameters
        ----------
        n_batches: int
            Number of batches
        background_noise_model: Optional[Callable]
            The probabilistic model used to generate noise. It should only have one input: size
        """

        self.n_markers = None

        self.background_noise_model = background_noise_model

        # Parameters related to batch effects and temporal effects
        self.n_batches = n_batches
        self.overall_batch_effects = {}
        self.local_batch_effects = {}
        self.temporal_effects = {}
        self.cell_abundances = {}

        # Parameters related to actual cell type models
        self.cell_type_labels_to_ids = {}
        self.cell_type_ids_to_labels = {}
        self.cell_types = {}

        self.cell_graph = GeneralCellGraph()

    def generate_cell_abundances(self,
                                 is_random: bool = True) -> None:
        """Generate cell abundances

        Parameters
        ----------
        is_random: bool
            Whether the cell abundances should be randomly generated
        """
        if is_random:
            # If randomly generate cell abundances,
            # we sample from a dirichlet(1)
            for b in range(self.n_batches):
                self.cell_abundances[b] = {}
                abundance = np.random.dirichlet(np.ones(len(self.cell_types)), size=1).reshape(-1)
                counter = 0
                for c_type in self.cell_types:
                    self.cell_abundances[b][c_type] = abundance[counter]
                    counter += 1
        else:
            # Otherwise, every cell will have an equal probability
            for b in range(self.n_batches):
                self.cell_abundances[b] = {}
                abundance = np.ones(len(self.cell_types))/len(self.cell_types)
                counter = 0
                for c_type in self.cell_types:
                    self.cell_abundances[b][c_type] = abundance[counter]
                    counter += 1


    def generate_overall_batch_effects(self,
                                       variance: float = 0.001) -> None:
        """Generate overall batch effects (main effects)

        Parameters
        ----------
        variance: float
            The variance of the effects
        """
        if self.n_batches == 1:
            # If we only have one batch, the effect is 0
            # (we won't be able to estimate it anyway in real life)
            self.overall_batch_effects[0] = 0
        else:
            # Otherwise, the main effects are generated from a normal distribution
            # We will make sure the effects add up to 0 (ANOVA: Effect coding)
            batch_effects = np.zeros(self.n_batches)
            batch_effects[:(self.n_batches - 1)] = np.random.normal(loc=0, scale=np.sqrt(variance),
                                                                    size=self.n_batches - 1)
            batch_effects[self.n_batches - 1] = -np.sum(batch_effects)
            for b in range(self.n_batches):
                self.overall_batch_effects[b] = batch_effects[b]


    def generate_local_batch_effects(self,
                                     variance: float = 0.001) -> None:
        """Generate local batch effects (interaction effects)
        
        We separate main effects from local effects since some methods are designed only
        to eliminate overall effects
        
        Parameters
        ----------
        variance: float
            The variance of the effects
        """
        if self.n_batches == 1:
            # If we only have one batch, the effects are 0
            # (we won't be able to estimate them anyway in real life)
            self.local_batch_effects[0] = {}
            for c_type in self.cell_type_labels_to_ids:
                self.local_batch_effects[0][c_type] = np.zeros(self.n_markers)
        else:
            # Otherwise, the effects are generated from a normal distribution
            # We will make sure the effects add up to 0 for each cell type as
            # well as for each marker (ANOVA: Effect coding)
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


    def generate_temporal_effects(self,
                                  variance: Optional[float] = None,
                                  coefficients: Optional[Union[list, np.ndarray]] = None,
                                  x: Optional[np.ndarray] = None,
                                  y: Optional[np.ndarray] = None,
                                  **kwargs) -> None:
        """Generate temporal effect

        Parameters
        ----------
        variance: float
            The variance of the end point if using Brownian bridge or polynomial
        coefficients: list or np.ndarray
            The coefficients of the polynomial to be generated
        x: np.ndarray
            The x values used to fit a spline
        y: np.ndarray
            The y values used to fit a spline
        kwargs: Extra parameters for the brownian bridge method or the spline function
        """
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

    def sample_one_batch(self,
                         n_samples: int,
                         cell_abundances: Optional[dict] = None,
                         batch: int = 0,
                         clip: bool = True) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        """Draw random samples for one batch

        Parameters
        ----------
        n_samples: int
            Number of samples
        cell_abundances: dict or None
            A dictionary whose keys are the cell labels. The corresponding values should be
            either the actual number of events for each cell type or the probability of each cell type
        batch: int
            The index of the batch for which we want to draw samples
        clip: bool
            Whether or not the resulting negative expressions should be clipped

        Returns
        -------
        expression_matrix: np.ndarray
            The expression matrix
        labels: np.ndarray
            The array of the corresponding cell type labels
        pseudo_time: np.ndarray
            The array of the positions on the differentiation paths
        children_cell_labels: np.ndarray
            The descendants to which the cells are differentiating towards
        """
        # If cell_abundances is not provided, we use the one stored in the object
        if cell_abundances is None:
            cell_abundances = self.cell_abundances[batch]

        # The idea of sampling is as follows:
        # Instead of drawing samples one by one
        # we sample cell types one by one.
        # And then we permute the expression matrix.

        # We first record the order of the cell types
        cell_type_order = list(cell_abundances.keys())
        # Record the number of events desired for each cell type
        n_per_cell_type = np.zeros(len(cell_abundances), dtype=int)

        if np.all([isinstance(i, int) and (i >= 0) for i in cell_abundances.values()]):
            # The values are actual numbers
            counter = 0
            for c_type in cell_type_order:
                n_per_cell_type[counter] = cell_abundances[c_type]
                counter += 1
        elif np.all([i >= 0 for i in cell_abundances.values()]) and \
                (np.sum(list(cell_abundances.values())) <= 1.5):
            # The values are probability-like
            cell_probabilities = np.zeros(len(cell_abundances))
            counter = 0
            for c_type in cell_type_order:
                cell_probabilities[counter] = cell_abundances[c_type]
                counter += 1
            cell_probabilities /= np.sum(cell_probabilities)
            n_per_cell_type = np.random.multinomial(n_samples, cell_probabilities)
        else:
            raise ValueError('Unknown cell abundances type')

        # The output arrays
        expression_matrix = np.zeros((n_samples, self.n_markers))
        expressed_index_matrix = np.zeros((n_samples, self.n_markers))
        pseudo_time = np.zeros((n_samples, self.n_markers))
        children_cell_labels = ["None"] * n_samples

        labels = [item for item, count in zip(cell_type_order, n_per_cell_type) for i in range(count)]

        # If overall batch effects have not been generated,
        # we set it to 0
        Psi_b = 0
        if batch in self.overall_batch_effects.keys():
            Psi_b = self.overall_batch_effects[batch]

        # We fill up the arrays using the "worm" method (I forget what it is called...)
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
            # If local effects have not been generated, we set them to 0
            Psi_bp = 0
            if batch in self.local_batch_effects.keys():
                Psi_bp = self.local_batch_effects[batch][c_type]
            G, T, children_labels = self.cell_graph.sample_graph(n, c_type)
            # Only expressed markers are subject to batch effects
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

        # If temporal effects have been generated
        if batch in self.temporal_effects.keys():
            temporal_effects = self.temporal_effects[batch][0](time_points)
            # Only expressed markers are subject to temporal effects
            expression_matrix += expressed_index_matrix * temporal_effects[:, np.newaxis]
        if clip:
            expression_matrix = np.clip(expression_matrix, a_min=0, a_max=None)

        E = 0
        if self.background_noise_model is not None:
            E = self.background_noise_model(size=(n_samples, self.n_markers))

        expression_matrix += E

        return expression_matrix, np.array(labels), pseudo_time, np.array(children_cell_labels)

    def sample(self,
               n_samples: Union[int, list, np.ndarray],
               cell_abundances: Optional[dict] = None,
               clip: bool = True) -> Tuple[dict, dict, dict, dict]:
        """Draw random samples for all batches

        Parameters
        ----------
        n_samples: int or list or np.ndarray
            Number of samples for each batch. If an integer is provided, then it will be used for all batches
        cell_abundances: dict or None
            A nested dictionary whose keys are the batches. The corresponding values should be
            a dictionary mapping cell types to cell numbers or probabilities OR
            It can be a plain dictionary whose keys are the cell labels. The corresponding values should be
            either the actual number of events for each cell type or the probability of each cell type
        clip: bool
            Whether or not the resulting negative expressions should be clipped

        Returns
        -------
        expression_matrices: dict
            The dictionary of expression matrices
        labels: dict
            The dictionary of arrays of the corresponding cell type labels
        pseudo_time: dict
            The dictionary of arrays of the positions on the differentiation paths
        children_cell_labels: dict
            The dictionary of descendants to which the cells are differentiating towards
                                
        """
        if cell_abundances is None:
            # If cell_abundances has not been generated, we generate it using the default
            # In emulation, the default will be the observed cell abundances
            # In creation, the default is random
            if len(self.cell_abundances) == 0:
                self.generate_cell_abundances()
            cell_abundances = self.cell_abundances

        # If the dictionary is not nested meaning that the values are actual
        # numbers or probabilities
        # we reuse it for every batch
        if not np.any([isinstance(i, dict) for i in cell_abundances.values()]):
            cell_abundances = {}
            for b in range(self.n_batches):
                cell_abundances[b] = deepcopy(cell_abundances)

        # If n_samples is an integer, we reuse it for every batch
        if isinstance(n_samples, int):
            n_samples = np.repeat(n_samples, self.n_batches)

        # Prepare the output dictionaries
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

