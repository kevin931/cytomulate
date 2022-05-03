# Math computation
import numpy as np

# Typing
from typing import Union, Optional, Any, List, Tuple, Callable


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

    def sample_cell(self,
                    n_samples: int,
                    clip: bool) -> Tuple[np.ndarray, np.ndarray]:
        """Draw random samples from the cell type model

        Parameters
        ----------
        n_samples: int
            Number of samples
        clip: bool
            Whether or not negative values should be clipped

        Returns
        -------
        np.ndarray, np.ndarray: The actual expression matrix and
                                an index array of positive values which
                                will be used during the actual sample
                                function in the CytofData object
        """
        X = np.zeros((n_samples, self.n_markers))
        X[:, :], _ = self.model.sample(n_samples)
        expressed_index = (X > 0)
        if clip:
            X = np.clip(X, a_min=0, a_max=None)
        return X, expressed_index
