# Math computation
import numpy as np


class CellType:
    def __init__(self, id, name, n_markers):
        self.id = id
        self.name = name
        self.parent = None
        self.children = []
        self.marker_pattern = np.zeros(n_markers)
        self.expression_level = np.zeros(n_markers)
        self.gating_markers = []
