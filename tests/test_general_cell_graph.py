import pytest
from cytomulate.cell_graph_general import GeneralCellGraph


@pytest.fixture
def Cell_Graph():
    return GeneralCellGraph()


def test_overall(Cell_Graph):
    assert (Cell_Graph.graph is None) and \
           (Cell_Graph.n_markers == -1) and \
           (isinstance(Cell_Graph.trajectories, dict))
