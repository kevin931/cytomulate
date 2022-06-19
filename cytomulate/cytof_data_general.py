# Math computation
import numpy as np

# List/Dict manipulation
from copy import deepcopy

# Trajectory functions
from cytomulate.utilities import trajectories
from cytomulate.cell_graph_general import GeneralCellGraph    

# Typing
from typing import Union, Optional, Tuple, Callable, Dict

OPT_PCK: Dict[str, bool] = {"PyCytoData": True}

try:
    from PyCytoData import PyCytoData
except ImportError:
    OPT_PCK["PyCytoData"] = False


class GeneralCytofData:
    def __init__(self,
                 n_batches: int = 1,
                 background_noise_model: Optional[Union[Callable, dict]] = None) -> None:
        """Initialize the GeneralCytofData object

        Parameters
        ----------
        n_batches: int
            Number of batches
        background_noise_model: Callable or dict
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
                                  coefficients: Optional[Union[dict, list, np.ndarray]] = None,
                                  x: Optional[Union[dict, np.ndarray]] = None,
                                  y: Optional[Union[dict, np.ndarray]] = None,
                                  **kwargs) -> None:
        """Generate temporal effect

        Parameters
        ----------
        variance: float
            The variance of the end point if using Brownian bridge or polynomial
        coefficients: dict, list or np.ndarray
            The coefficients of the polynomial to be generated or a dictionary of coefficients of the polynomials to be generated
        x: dict or np.ndarray
            The x values used to fit a spline or a dictionary of x values used to fit a spline
        y: dict or np.ndarray
            The y values used to fit a spline or dictionary of y values used to fit a spline
        kwargs: Extra parameters for the brownian bridge method or the spline function
        """
        if (coefficients is not None) and (not isinstance(coefficients, dict)):
            # This means that coefficients is a list or an np.ndarray
            coefficients_copy = deepcopy(coefficients)
            coefficients = {}
            for b in range(self.n_batches):
                coefficients[b] = deepcopy(coefficients_copy)

        if (x is not None) and (not isinstance(x, dict)):
            x_copy = deepcopy(x)
            x = {}
            for b in range(self.n_batches):
                x[b] = deepcopy(x_copy)

        if (y is not None) and (not isinstance(y, dict)):
            y_copy = deepcopy(y)
            y = {}
            for b in range(self.n_batches):
                y[b] = deepcopy(y_copy)

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
                         beta_alpha: Union[float, int] = 0.4,
                         beta_beta: Union[float, int] = 1.0) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
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
        beta_alpha: float or int
            The alpha parameter of the beta distribution
        beta_beta: float or int
            The beta parameter of the beta distribution

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
            n_samples = np.sum(list(cell_abundances.values()))
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
            X, expressed_index = self.cell_types[c_type].sample_cell(n)
            # If local effects have not been generated, we set them to 0
            Psi_bp = 0
            if batch in self.local_batch_effects.keys():
                Psi_bp = self.local_batch_effects[batch][c_type]
            G, T, children_labels = self.cell_graph.sample_graph(n, c_type, beta_alpha, beta_beta)
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

        expression_matrix = np.clip(expression_matrix, a_min=0, a_max=None)

        E = 0
        if self.background_noise_model is not None:
            if isinstance(self.background_noise_model, dict):
                E = self.background_noise_model[batch](size=(n_samples, self.n_markers))
            else:
                E = self.background_noise_model(size=(n_samples, self.n_markers))

        expression_matrix += E

        return expression_matrix, np.array(labels), pseudo_time, np.array(children_cell_labels)

    def sample(self,
               n_samples: Union[int, list, np.ndarray],
               cell_abundances: Optional[dict] = None,
               beta_alpha: Union[float, int, dict] = 0.4,
               beta_beta: Union[float, int, dict] = 0.4) -> Tuple[dict, dict, dict, dict]:
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
        beta_alpha: float, int, or dict
            The alpha parameters of the beta distribution
        beta_beta: float, int, or dict
            The beta parameters of the beta distribution

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
            cell_abundances_copy = deepcopy(cell_abundances)
            cell_abundances = {}
            for b in range(self.n_batches):
                cell_abundances[b] = deepcopy(cell_abundances_copy)

        # Then we make sure the keys of the supplied dictionary
        # coincide with the one we have in cell_types
        cell_abundances_copy = deepcopy(cell_abundances)
        cell_abundances = {}
        for b in range(self.n_batches):
            cell_abundances[b] = {}
            for c_type in self.cell_types:
                if c_type in cell_abundances_copy[b]:
                    cell_abundances[b][c_type] = cell_abundances_copy[b][c_type]
                else:
                    cell_abundances[b][c_type] = 0

        # If n_samples is an integer, we reuse it for every batch
        if isinstance(n_samples, int):
            n_samples = np.repeat(n_samples, self.n_batches)

        if isinstance(beta_alpha, float) or isinstance(beta_alpha, int):
            beta_alpha_copy = beta_alpha
            beta_alpha = {}
            for b in range(self.n_batches):
                beta_alpha[b] = beta_alpha_copy

        if isinstance(beta_beta, float) or isinstance(beta_beta, int):
            beta_beta_copy = beta_beta
            beta_beta = {}
            for b in range(self.n_batches):
                beta_beta[b] = beta_beta_copy

        if not isinstance(self.background_noise_model, dict):
            background_noise_model_copy = deepcopy(self.background_noise_model)
            background_noise_model = {}
            for b in range(self.n_batches):
                background_noise_model[b] = deepcopy(background_noise_model_copy)

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
                beta_alpha[b],
                beta_beta[b])

        return expression_matrices, labels, pseudo_time, children_cell_labels
    
    
    def sample_to_pycytodata(self,
                             n_samples: Union[int, list, np.ndarray],
                             cell_abundances: Optional[dict] = None,
                             beta_alpha: Union[float, int, dict] = 0.4,
                             beta_beta: Union[float, int, dict] = 0.4) -> "PyCytoData": #type: ignore
        """Draw random samples for all batches and returns a PyCytoData object.
        
        This method is a wrapper for the ``sample`` method but provides an interface to
        return a ``PyCytoData`` object.

        Parameters
        ----------
        n_samples: int or list or np.ndarray
            Number of samples for each batch. If an integer is provided, then it will be used for all batches
        cell_abundances: dict or None
            A nested dictionary whose keys are the batches. The corresponding values should be
            a dictionary mapping cell types to cell numbers or probabilities OR
            It can be a plain dictionary whose keys are the cell labels. The corresponding values should be
            either the actual number of events for each cell type or the probability of each cell type
        beta_alpha: float, int, or dict
            The alpha parameters of the beta distribution
        beta_beta: float, int, or dict
            The beta parameters of the beta distribution

        Returns
        -------
        pcd: PyCytoData
            A PyCytoData object with the simulated data.
            
        Raises
        -------
        ImportError: No ``PyCytoData`` installation present.
            
        Note
        -----
        The ``PyCytoData`` is not compatible with storing ``pseudo_time`` and ``children_cell_labels``. If you would
        like these information, use the traditional ``sample`` method instead.
        
        Note
        -----
        ``PyCytoData`` is an optional dependency. If an ``ImportError`` is raised, you need to install the
        the package first. Tutorials can be found here: https://pycytodata.readthedocs.io/en/latest/installation.html.
        """

        if not OPT_PCK["PyCytoData"]:
            raise ImportError("Error importing 'PyCytoData'. Install it first if you haven't done so.")

        exprs, labels, _, _ = self.sample(n_samples=n_samples, 
                                          cell_abundances=cell_abundances,
                                          beta_alpha=beta_alpha,
                                          beta_beta=beta_beta)
        
        sample_index: np.ndarray = np.repeat("0", n_samples)
        pcd: PyCytoData = PyCytoData(expression_matrix=exprs[0], cell_types=labels[0], sample_index=sample_index)
        
        if self.n_batches > 1:
            b: int
            for b in range(self.n_batches):
                if b == 0:
                    continue
                sample_index = np.repeat(str(b), n_samples)
                pcd.add_sample(expression_matrix=exprs[b], cell_types=labels[b], sample_index=sample_index)
                
        return pcd