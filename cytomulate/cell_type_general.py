# Math computation
import numpy as np

# Typing
from typing import Union, Optional, Any, List, Callable


class GeneralCellType:
    def __init__(self, label, cell_id, n_markers):
        self.label = label
        self.id = cell_id

        self.n_markers = n_markers
        self.markers = np.arange(self.n_markers)

        self.model = None

    def sample_cell(self, n_samples, clip):
        X = np.zeros((n_samples, self.n_markers))
        X[:, :], _ = self.model.sample(n_samples)
        expressed_index = (X > 0)
        if clip:
            X = np.clip(X, a_min=0, a_max=None)
        return X, expressed_index
