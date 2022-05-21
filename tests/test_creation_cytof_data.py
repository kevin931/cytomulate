import pytest
import networkx as nx
from cytomulate.creation.cytof_data import CreationCytofData
from cytomulate.utilities import univariate_noise_model


@pytest.fixture
def Cytof_Data():
    """Returns a CreationCytofData instance"""
    model = univariate_noise_model("normal", loc = 0, scale = 1)
    return CreationCytofData(n_batches=2,
                             n_types=20,
                             n_markers=30,
                             n_trees=2,
                             background_noise_model=model)


def test_n_markers(Cytof_Data):
    assert Cytof_Data.n_markers == 30


def test_background_noise_model(Cytof_Data):
    temp = Cytof_Data.background_noise_model((5, 3))
    assert temp.shape == (5, 3)


def test_n_batches(Cytof_Data):
    assert Cytof_Data.n_batches == 2


def test_n_types(Cytof_Data):
    assert len(Cytof_Data.cell_types) == 20


def test_n_trees(Cytof_Data):
    assert nx.number_weakly_connected_components(Cytof_Data.cell_graph.graph) == 2


def test_cell_types(Cytof_Data):
    Cytof_Data.initialize_cell_types()
    temp = Cytof_Data.cell_types[0].sample_cell(2, True)
    assert temp[0].shape == (2, 30)


def test_cell_graph(Cytof_Data):
    Cytof_Data.initialize_cell_types()
    Cytof_Data.generate_cell_graph()
    temp = Cytof_Data.sample(n_samples=1)
    assert isinstance(temp[3], dict)



