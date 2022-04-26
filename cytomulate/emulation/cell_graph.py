# Math computation
import numpy as np
import itertools

# Graph package and functions
import networkx as nx
from networkx.algorithms import tree
from networkx.algorithms.community import greedy_modularity_communities

# Superclass
from cell_graph_general import GeneralCellGraph


class EmulationCellGraph(GeneralCellGraph):
    def __init__(self):
        super().__init__()
        self.complete_undirected_graph = None
        self.bead_label = None

    def initialize_graph(self, cell_types, bead_label=None):
        self.complete_undirected_graph = nx.Graph()
        self.bead_label = bead_label
        for label in cell_types:
            self.complete_undirected_graph.add_node(label)

        for labels in itertools.combinations(cell_types.keys(), 2):
            if (labels[0] == bead_label) or (labels[1] == bead_label):
                continue
            else:
                mean1 = cell_types[labels[0]].observed_mean
                mean2 = cell_types[labels[1]].observed_mean
                self.complete_undirected_graph.add_edge(labels[0], labels[1], weight=np.linalg.norm(mean1 - mean2))

    def prune_graph(self, graph_topology="bidirectional complete"):
        if "complete" in graph_topology:
            self.graph = self.complete_undirected_graph.to_directed()
            if graph_topology == "complete":
                nodes = np.array(list(set(self.graph.nodes) - {self.bead_label}))
                np.random.shuffle(nodes)
                for labels in itertools.combinations(nodes, 2):
                    self.graph.remove_edge(labels[0], labels[1])
            else:
                if graph_topology != "bidirectional complete":
                    raise ValueError('Unknown graph type')
        elif "cyclic" in graph_topology:
            nodes = np.array(list(set(self.complete_undirected_graph.nodes) - {self.bead_label}))
            np.random.shuffle(nodes)
            self.graph = nx.DiGraph()
            if self.bead_label is not None:
                self.graph.add_node(self.bead_label)
            for i in range(len(nodes)):
                if i == len(nodes) - 1:
                    self.graph.add_edge(nodes[i], nodes[0])
                    if graph_topology == "bidirectional cyclic":
                        self.graph.add_edge(nodes[0], nodes[i])
                else:
                    self.graph.add_edge(nodes[i], nodes[i+1])
                    if graph_topology == "bidirectional cyclic":
                        self.graph.add_edge(nodes[i+1], nodes[i])
        elif graph_topology in ["tree", "forest"]:
            mst = tree.minimum_spanning_edges(self.complete_undirected_graph, algorithm="kruskal", weight="weight", data=False)
            mst_G = nx.Graph()
            mst_G.add_edges_from(list(mst))
            if graph_topology == "tree":
                root = set(mst_G.nodes).pop()
                self.graph = nx.dfs_tree(mst_G, root)
            elif graph_topology == "forest":
                forest = list(greedy_modularity_communities(mst_G))
                tree_list = []
                for t in range(len(forest)):
                    nodes = list(forest[t])
                    root = np.random.choice(nodes)
                    tree_list.append(nx.dfs_tree(mst_G.subgraph(nodes), root))
                self.graph = nx.compose_all(tree_list)
            else:
                raise ValueError('Unknown graph type')
            if self.bead_label is not None:
                self.graph.add_node(self.bead_label)
        else:
            raise ValueError('Unknown graph type')
