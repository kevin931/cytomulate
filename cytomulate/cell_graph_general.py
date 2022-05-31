# Math computation
import numpy as np

# Graph package and functions
import networkx as nx
import matplotlib.pyplot as plt
from cytomulate.utilities import trajectories

# Typing
from typing import Union, Tuple


class GeneralCellGraph:
    def __init__(self) -> None:
        """Initialize the GeneralCellGraph object

        """
        self.graph = None
        self.n_markers = -1
        # trajectories would be a dictionary whose keys are edges of a
        # directed graph and whose values would be a list of functions
        # that describes the actual differentiation path
        self.trajectories = {}

    def generate_trajectories(self,
                              cell_types: dict,
                              **kwargs) -> None:
        """Generate the actual differential paths

        Parameters
        ----------
        cell_types: dict
            A dictionary of CellType objects
        kwargs:
            Extra parameters needed for non-default path generation algorithms
        """
        edges = self.graph.edges
        for e in edges:
            from_label = e[0]
            to_label = e[1]
            end_values = cell_types[to_label].cell_mean - cell_types[from_label].cell_mean
            if self.n_markers <= 0:
                self.n_markers = len(end_values)
            self.trajectories[e] = trajectories(end_values=end_values, **kwargs)

    def sample_graph(self,
                     n_samples: int,
                     cell_label: Union[str, int],
                     beta_alpha: Union[float, int] = 0.4,
                     beta_beta: Union[float, int] = 1.0) -> Tuple[np.ndarray, np.ndarray, list]:
        """Draw random samples of a cell type from the cell differentiation graph

        Parameters
        ----------
        n_samples: int
            Number of samples
        cell_label: str or int
            The label of the cell needed
        beta_alpha: float
            The alpha parameter of the beta distribution
        beta_beta: float
            The beta parameter of the beta distribution

        Returns
        -------
        G: np.ndarray
            The additive values of the path
        pseudo_time: np.ndarray
            The pseudo times
        labels: list
           The cell types to which the cell is differentiating
        """
        if len(self.trajectories) < 1:
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
                    p_time = np.random.beta(beta_alpha, beta_beta, n)
                    G[start_n: end_n, m] = self.trajectories[(cell_label, c_label)][m](p_time)
                    pseudo_time[start_n: end_n, m] = p_time
                start_n += n

        return G, pseudo_time, labels


    def visualize_graph(self) -> None:
        """Visualize the cell graph

        """
        connected_components = list(nx.connected_components(self.graph.to_undirected()))
        n_plt = 0
        figs = []
        for nodes in connected_components:
            G = self.graph.subgraph(list(nodes))
            pos = nx.planar_layout(G, scale=20)
            d = dict(G.degree)
            nodelist = list(d.keys())
            nodesize = []
            max_degree = np.max(list(d.values()))

            for k in d:
                predecessors = list(G.predecessors(k))
                if len(predecessors) == 0:
                    nodesize.append(max_degree * 150)
                else:
                    nodesize.append(d[k] * 100)

            colors = []
            for node in G.nodes():
                successors = list(G.successors(node))
                predecessors = list(G.predecessors(node))

                if len(successors) == 0:
                    colors.append("springgreen")
                elif len(predecessors) == 0:
                    colors.append("royalblue")
                else:
                    colors.append("magenta")
            figs.append(plt.figure(n_plt))
            ax = figs[n_plt].add_subplot(1,1,1)
            ax.set_title("Tree " + str(n_plt))
            nx.draw(G, pos,
                    with_labels=True,
                    nodelist=nodelist,
                    node_size=nodesize,
                    node_color=colors)

            n_plt += 1
        plt.show()