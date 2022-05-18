# Math computation
import numpy as np
import itertools

# Graph package and functions
import networkx as nx
from networkx.algorithms import tree

# Superclass
from cytomulate.cell_graph_general import GeneralCellGraph


class CreationCellGraph(GeneralCellGraph):
    def __init__(self):
        """Initialize the CreationCellGraph object

        """
        super().__init__()
        self.serialized_graph = []

    def initialize_graph(self,
                         cell_types: dict,
                         n_trees: int) -> None:
        """Generate a random graph that consists of one or several trees

        Parameters
        ----------
        cell_types: dict
            A dictionary of CreationCellType objects
        n_trees: int
            Number of trees to be simulated
        """
        # If n_trees is negative, then we assume that the user do not want any tree structure
        # Therefore, every cell type becomes a tree of size 1
        if n_trees <= 0:
            n_trees = len(cell_types)

        G = nx.Graph()
        for label in cell_types:
            G.add_node(label)

        # Since the tree is randomly generated, we will first
        # construct a complete graph with randomly generated weights
        for labels in itertools.combinations(cell_types.keys(), 2):
            G.add_edge(labels[0], labels[1], weight=np.random.uniform(size=1)[0])

        # And then we randomly group cell types together
        # Those that are grouped together will form a tree
        nodes_list = np.array(G.nodes)
        np.random.shuffle(nodes_list)
        nodes = np.array_split(nodes_list, n_trees)
        tree_list = []
        for t in nodes:
            # For each group, we use Kruskal to find a MST
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
        # We use the topological sort on the graph to facilitate future cell type generation
        self.serialized_graph = list(nx.topological_sort(self.graph))
