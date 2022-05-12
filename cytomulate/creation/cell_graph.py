# Math computation
import numpy as np
import itertools

# Graph package and functions
import networkx as nx
from networkx.algorithms import tree

# Superclass
from cell_graph_general import GeneralCellGraph

# Typing
from typing import Union, Optional, Any, List, Tuple, Callable


class CreationCellGraph(GeneralCellGraph):
    def __init__(self):
        super().__init__()
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
