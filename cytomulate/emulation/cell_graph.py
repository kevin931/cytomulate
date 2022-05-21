# Math computation
import numpy as np
import itertools

# Graph package and functions
import networkx as nx
from networkx.algorithms import tree
from networkx.algorithms.community import greedy_modularity_communities

# Superclass
from cytomulate.cell_graph_general import GeneralCellGraph

# Typing
from typing import Union, Optional


class EmulationCellGraph(GeneralCellGraph):
    def __init__(self) -> None:
        """Initialize the EmulationCellGraph object

        """
        super().__init__()
        self.complete_undirected_graph = None
        self.bead_label = None


    def initialize_graph(self,
                         cell_types: dict,
                         bead_label: Optional[Union[str, int]] = None) -> None:
        """Build a complete graph among all cell types

        Parameters
        ----------
        cell_types: dict
            A dictionary of CellType objects
        bead_label: str or int
            The label for beads
        """
        self.complete_undirected_graph = nx.Graph()
        self.bead_label = bead_label
        for label in cell_types:
            self.complete_undirected_graph.add_node(label)

        # The edges consist of n choose 2 elements
        for labels in itertools.combinations(cell_types.keys(), 2):
            if (labels[0] == bead_label) or (labels[1] == bead_label):
                # If any of the "cell types" is bead, we skip
                # since beads should not differentiate
                continue
            else:
                # The weights are the l2 distance between two means
                mean1 = cell_types[labels[0]].cell_mean
                mean2 = cell_types[labels[1]].cell_mean
                self.complete_undirected_graph.add_edge(labels[0], labels[1], weight=np.linalg.norm(mean1 - mean2))


    def prune_graph(self,
                    graph_topology: str = "tree") -> None:
        """Construct a tree or a forest from the complete graph

        Parameters
        ----------
        graph_topology: str
            The type of the graph desired. Should be either tree or forest
        """
        if graph_topology in ["tree", "forest"]:
            # We will first construct a minimum spanning tree
            mst = tree.minimum_spanning_edges(self.complete_undirected_graph, algorithm="kruskal", weight="weight", data=False)
            mst_G = nx.Graph()
            mst_G.add_edges_from(list(mst))
            if graph_topology == "tree":
                # If a single tree is desired, we randomly select a node as the root
                # and construct a directed graph
                root = set(mst_G.nodes).pop()
                self.graph = nx.dfs_tree(mst_G, root)
            else:
                # If a forest is desired, we use community detection algorithm
                # to divide the trees
                forest = list(greedy_modularity_communities(mst_G))
                tree_list = []
                for t in range(len(forest)):
                    nodes = list(forest[t])
                    root = np.random.choice(nodes)
                    tree_list.append(nx.dfs_tree(mst_G.subgraph(nodes), root))
                self.graph = nx.compose_all(tree_list)

            # Finally, if bead is present, we add bead to the resulting graph
            if self.bead_label is not None:
                self.graph.add_node(self.bead_label)
        else:
            raise ValueError('Unknown graph type')
