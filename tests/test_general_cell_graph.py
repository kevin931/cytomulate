import pytest
import numpy as np
from cytomulate.cell_graph_general import GeneralCellGraph
from cytomulate.creation.cytof_data import CreationCytofData


@pytest.fixture
def Cell_Graph():
    return GeneralCellGraph()


def test_initial(Cell_Graph):
    assert (Cell_Graph.graph is None) and \
           (Cell_Graph.n_markers == -1) and \
           (isinstance(Cell_Graph.trajectories, dict))


def test_generate_trajectories(Cell_Graph):
    y = CreationCytofData()
    y.initialize_cell_types()
    Cell_Graph.graph = y.cell_graph.graph
    Cell_Graph.generate_trajectories(cell_types=y.cell_types)

    assert isinstance(Cell_Graph.trajectories, dict)


def test_sample_graph(Cell_Graph):
    y = CreationCytofData()
    y.initialize_cell_types()
    Cell_Graph.graph = y.cell_graph.graph
    Cell_Graph.generate_trajectories(cell_types=y.cell_types)
    temp = Cell_Graph.sample_graph(n_samples=1, cell_label=0)

    assert isinstance(temp[0], np.ndarray)


def test_visualize_graph(Cell_Graph, mocker):
    plt_mock = mocker.MagicMock()
    mocker.patch("matplotlib.pyplot.show", plt_mock)
    y = CreationCytofData()
    y.initialize_cell_types()
    Cell_Graph.graph = y.cell_graph.graph
    Cell_Graph.visualize_graph()
    plt_mock.assert_called_once()
