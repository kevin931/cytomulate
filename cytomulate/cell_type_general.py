# Math computation
import numpy as np

# Typing
from typing import Union, Tuple


class GeneralCellType:
    def __init__(self,
                 label: Union[str, int],
                 cell_id: int,
                 n_markers: int) -> None:
        """Initialize the GeneralCellType object

        Parameters
        ----------
        label: str or int
            The label (name) for the cell type
        cell_id: int
            The id number assigned to the cell type
        n_markers: int
            Number of markers used in the experiment
        """
        self.label = label
        self.id = cell_id

        self.n_markers = n_markers
        # We generate the markers as a sequence of numbers
        self.markers = np.arange(self.n_markers)

        # The actual cell type model
        self.model = None

        # cell_mean and cell_covariance are used during cell differentiation
        self.cell_mean = np.zeros(self.n_markers)
        self.cell_covariance = np.zeros((self.n_markers, self.n_markers))
        # zero_probabilities is used for adjustment 
        self.zero_probabilities = np.zeros(n_markers)
        
    def sample_cell(self,
                    n_samples: int) -> Tuple[np.ndarray, np.ndarray]:
        """Draw random samples from the cell type model

        Parameters
        ----------
        n_samples: int
            Number of samples

        Returns
        -------
        X: np.ndarray
            The actual expression matrix
        expressed_index: np.ndarray:  and
            An index array of positive values which
            will be used during the actual sample
            function in the CytofData object
                                
        """
        X = np.zeros((n_samples, self.n_markers))
        X[:, :], _ = self.model.sample(n_samples)
        X = np.clip(X, a_min=0, a_max=None)
        for m in range(self.n_markers):
            n_zero_exp = int((self.zero_probabilities[m]) * n_samples) 
            n_zero_present = np.sum(X[:, m]<0.0001)
            n_zero_needed = np.max([0, n_zero_exp-n_zero_present])
            if n_zero_needed > 0:
                non_zero_ind = np.where(X[:,m]>=0.0001)[0]
                p = 1/(X[non_zero_ind, m])
                p /= np.sum(p)
                # if n_zero_needed is 0, this should yield 
                # [] which when plugged into the next statement
                # shall change nothing 
                ind_to_zero = np.random.choice(non_zero_ind, size=n_zero_needed,
                                                replace=False, p=p)
                X[ind_to_zero, m] = 0
                
        expressed_index = (X > 0)
        return X, expressed_index
