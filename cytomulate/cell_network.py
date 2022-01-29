import numpy as np
import networkx as nx
from networkx.algorithms import tree
from networkx.algorithms.community import greedy_modularity_communities


class CellNetwork:
    def __init__(self):
        self.network = None
        self.trajectories = {}

    def initialize_network(self, cell_types, bead_label):
        self.network = nx.Graph()
        for label in cell_types:
            self.network.add_node(label)

        for label1 in cell_types:
            mean1 = cell_types[label1].observed_mean
            for label2 in cell_types:
                mean2 = cell_types[label2].observed_mean
                if (label1 == label2) or (label1 == bead_label) or (label2 == bead_label):
                    continue
                else:
                    self.network.add_edge(label1, label2, weight=np.linalg.norm(mean1 - mean2))

    def prune_network(self, network_type="complete"):
        if network_type == "complete":
            pass
        else:
            raise ValueError('Unknown network type')


    def generate_trajectories(self):
        pass

    def sample_network(self, n_samples, cell_label, cell_types):
        pass