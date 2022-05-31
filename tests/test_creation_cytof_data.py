import pytest
import networkx as nx
import numpy as np
from cytomulate.creation.cytof_data import CreationCytofData
from cytomulate.utilities import univariate_noise_model


@pytest.fixture
def Cytof_Data():
    model = univariate_noise_model("normal", loc = 0, scale = 1)
    return CreationCytofData(n_batches=2,
                             n_types=20,
                             n_markers=30,
                             n_trees=2,
                             background_noise_model=model)


@pytest.fixture
def Cytof_Data_No_Diff():
    model = univariate_noise_model("normal", loc = 0, scale = 1)
    return CreationCytofData(n_batches=2,
                             n_types=20,
                             n_markers=30,
                             n_trees=-1,
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
    temp = Cytof_Data.cell_types[0].sample_cell(2)
    assert temp[0].shape == (2, 30)


def test_cell_types_non_tree(Cytof_Data):
    Cytof_Data.cell_graph.graph.add_edge(15, 0)
    Cytof_Data.cell_graph.graph.add_edge(10, 0)
    try:
        Cytof_Data.initialize_cell_types()
    except ValueError:
        assert True


def test_cell_graph(Cytof_Data):
    Cytof_Data.initialize_cell_types()
    Cytof_Data.generate_cell_graph()
    temp = Cytof_Data.cell_graph.sample_graph(n_samples=1, cell_label=0)
    assert isinstance(temp[0], np.ndarray)


def test_cell_graph_no_diff(Cytof_Data_No_Diff):
    Cytof_Data_No_Diff.initialize_cell_types()
    Cytof_Data_No_Diff.generate_cell_graph()
    temp = Cytof_Data_No_Diff.cell_graph.sample_graph(n_samples=1, cell_label=0)
    assert temp[0] == 0


def test_overall_batch_effects(Cytof_Data):
    Cytof_Data.initialize_cell_types()
    Cytof_Data.generate_overall_batch_effects(variance=0.1)
    temp = Cytof_Data.sample(n_samples=1)

    assert (isinstance(temp[0], dict)) and \
           (np.sum(list(Cytof_Data.overall_batch_effects.values())) == 0)


def test_local_batch_effects(Cytof_Data):
    Cytof_Data.initialize_cell_types()
    Cytof_Data.generate_local_batch_effects(variance=0.1)
    temp = Cytof_Data.sample(n_samples=1)

    assert (isinstance(temp[0], dict)) and \
           (np.sum(Cytof_Data.local_batch_effects[0][0]) == 0)


def test_temporal_effects(Cytof_Data):
    Cytof_Data.initialize_cell_types()
    Cytof_Data.generate_temporal_effects(variance=0.1)
    temp = Cytof_Data.sample(n_samples=1)

    assert (isinstance(temp[0], dict)) and \
           (Cytof_Data.temporal_effects[0][0](0) == 0)


def test_cell_abundances(Cytof_Data):
    Cytof_Data.initialize_cell_types()
    Cytof_Data.generate_cell_abundances()

    assert np.sum(list(Cytof_Data.cell_abundances[0].values())) <= 1.5


def test_overall(Cytof_Data):
    Cytof_Data.initialize_cell_types()
    Cytof_Data.generate_cell_graph()
    Cytof_Data.generate_overall_batch_effects(variance=0.1)
    Cytof_Data.generate_local_batch_effects(variance=0.1)
    Cytof_Data.generate_temporal_effects(variance=0.1)
    Cytof_Data.generate_cell_abundances()
    temp = Cytof_Data.sample(n_samples=1)

    assert (isinstance(temp[3], dict))
