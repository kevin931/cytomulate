# Math computation
import numpy as np

# Statistical models
from sklearn.mixture import GaussianMixture

# Typing
from typing import Union

# Superclass
from cytomulate.cell_type_general import GeneralCellType


class EmulationCellType(GeneralCellType):
    def __init__(self,
                 label: Union[str, int],
                 cell_id: int,
                 n_markers: int) -> None:
        """Initialize the EmulationCellType object
        
        Parameters
        ----------
        label: str or int
            The label (name) for the cell type
        cell_id: int
            The id number assigned to the cell type
        n_markers: int
            Number of markers used in the experiment
        """
        super().__init__(label, cell_id, n_markers)


    def fit(self,
            data: np.ndarray,
            max_components: int,
            min_components: int,
            covariance_types: Union[list, tuple]) -> None:
        """Fit cell type models using the data provided
        
        The model selection is done using BIC
        
        Parameters
        ----------
        data: np.ndarray
            The expression matrix corresponding to the target cell type
        max_components: int
            The maximum number of components to be included in the GMM
        min_components: int
            The minimum number of components to be included in the GMM
        covariance_types: list or tuple
            A list of strings specifying the type of covariances to be considered during the model fitting
        """

        # We first get the number of data points
        n = data.shape[0]
        # If the number of component the greater than the number of points
        # we will change them to the number of points
        min_components = np.min([min_components, n])
        max_components = np.min([max_components, n])

        self.cell_mean = np.mean(data, axis=0)
        self.cell_covariance = np.cov(data, rowvar=False)

        # We use BIC (the smaller the better) to perform model selection
        smallest_bic = np.Inf
        current_bic = 0
        best_gm = None
        for n_components in range(min_components, max_components + 1):
            for cv_type in covariance_types:
                gm = GaussianMixture(n_components=n_components,
                                     covariance_type=cv_type).fit(data)
                current_bic = gm.bic(data)
                if current_bic < smallest_bic:
                    smallest_bic = current_bic
                    best_gm = gm

        self.model = best_gm

        # The EM algorithm will not give 1 so we normalize it
        if self.model.n_components == 1:
            self.model.weights_[0] = 1.

