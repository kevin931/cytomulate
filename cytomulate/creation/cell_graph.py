# Math computation
import numpy as np

# Graph package and functions
import networkx as nx
from networkx.algorithms import tree
from creation.utilities import trajectories

#
import itertools


class CellGraph:
    def __init__(self):
        self.graph = None
        self.n_markers = -1
        self.trajectories = {}
        self.serialized_graph = []

    def initialize_graph(self, cell_types, n_trees):
        if n_trees <= 0:
            n_trees = len(cell_types)

        G = nx.Graph()
        for label in cell_types:
            G.add_node(label)

        for labels in itertools.combinations(cell_types.keys(), 2):
            G.add_edge(labels[0], labels[1], weight=np.random.uniform(size=1)[0])

        nodes_list = np.array(G.nodes)
        np.random.shuffle(nodes_list)
        nodes = np.array_split(nodes_list, n_trees)
        tree_list = []
        for t in nodes:
            if len(t) > 1:
                root = np.random.choice(t)
                mst = tree.minimum_spanning_edges(G.subgraph(t), algorithm="kruskal", weight="weight", data=False)
                mst_G = nx.Graph()
                mst_G.add_edges_from(list(mst))
                tree_list.append(nx.dfs_tree(mst_G, root))
            elif len(t) == 1:
                root = np.random.choice(t)
                mst_G = nx.Graph()
                mst_G.add_node(root)
                tree_list.append(nx.dfs_tree(mst_G, root))
            else:
                continue
        self.graph = nx.compose_all(tree_list)
        self.serialized_graph = list(nx.topological_sort(self.graph))

    def generate_trajectories(self, cell_types, **kwargs):
        edges = self.graph.edges
        for e in edges:
            from_label = e[0]
            to_label = e[1]
            end_values = cell_types[to_label].cell_mean - cell_types[from_label].cell_mean
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
        pass