# Math computation
import numpy as np

# Graph package and functions
import networkx as nx
from utilities import trajectories


class GeneralCellGraph:
    def __init__(self):
        self.graph = None
        self.n_markers = -1
        self.trajectories = {}

    def generate_trajectories(self, cell_types, **kwargs):
        edges = self.graph.edges
        for e in edges:
            from_label = e[0]
            to_label = e[1]
            end_values = cell_types[to_label].observed_mean - cell_types[from_label].observed_mean
            if self.n_markers <= 0:
                self.n_markers = len(end_values)
            self.trajectories[e] = trajectories(end_values=end_values, **kwargs)

    def sample_graph(self, n_samples, cell_label):
        if self.graph is None:
            return 0, 0, ["None"] * n_samples

        children_cell_labels = list(self.graph.successors(cell_label))
        n_children = len(children_cell_labels)
        labels = ["None"] * n_samples

        G = np.zeros((n_samples, self.n_markers))
        pseudo_time = np.zeros((n_samples, self.n_markers))

        if n_children >= 1:
            n_per_child = np.random.multinomial(n_samples, np.ones(n_children)/n_children)
            labels = [item for item, count in zip(children_cell_labels, n_per_child) for i in range(count)]
            start_n = 0
            end_n = 0
            counter = 0
            for c_label in children_cell_labels:
                n = n_per_child[counter]
                counter += 1
                if n == 0:
                    continue
                end_n += n
                for m in range(self.n_markers):
                    p_time = np.random.beta(0.4,1,n)
                    G[start_n: end_n, m] = self.trajectories[(cell_label, c_label)][m](p_time)
                    pseudo_time[start_n: end_n, m] = p_time
                start_n += n

        return G, pseudo_time, labels

    def visualize_graph(self):
        nx.draw(self.graph, with_labels=True)