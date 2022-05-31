import pytest
import networkx as nx
import numpy as np
from cytomulate.emulation.cytof_data import EmulationCytofData
from cytomulate.utilities import univariate_noise_model


@pytest.fixture
def Cytof_Data():
    model = univariate_noise_model("normal", loc = 0, scale = 1)
    return EmulationCytofData(n_batches=2,
                              background_noise_model=model)


def test_background_noise_model(Cytof_Data):
    temp = Cytof_Data.background_noise_model((5, 3))
    assert temp.shape == (5, 3)


def test_n_batches(Cytof_Data):
    assert Cytof_Data.n_batches == 2


def test_bead_label(Cytof_Data):
    assert Cytof_Data.bead_label is None


def test_cell_types(Cytof_Data):
    expression_matrix = np.random.multivariate_normal(mean=np.zeros(30),
                                                      cov=np.eye(30),
                                                      size = 500)
    labels = np.random.binomial(n=1, p=0.5, size=500)
    Cytof_Data.initialize_cell_types(expression_matrix=expression_matrix,
                                     labels=labels,
                                     max_components=3,
                                     min_components=3,
                                     covariance_types=["diag"])
    temp = Cytof_Data.cell_types[0].sample_cell(2)
    assert temp[0].shape == (2, 30)


def test_cell_graph_tree(Cytof_Data):
    expression_matrix = np.random.multivariate_normal(mean=np.zeros(30),
                                                      cov=np.eye(30),
                                                      size=1000)
    labels = np.random.choice(a=[0,1,2,3,4], size=1000, replace=True)
    Cytof_Data.initialize_cell_types(expression_matrix=expression_matrix,
                                     labels=labels,
                                     max_components=3,
                                     min_components=3,
                                     covariance_types=["diag"])

    Cytof_Data.generate_cell_graph(graph_topology="tree")
    temp = Cytof_Data.cell_graph.sample_graph(n_samples=1, cell_label=0)

    assert (nx.number_weakly_connected_components(Cytof_Data.cell_graph.graph) == 1) and \
           (isinstance(temp[0], np.ndarray))


def test_cell_graph_forest(Cytof_Data):
    expression_matrix = np.random.multivariate_normal(mean=np.zeros(30),
                                                      cov=np.eye(30),
                                                      size=1000)
    labels = np.random.choice(a=[0,1,2,3,4], size=1000, replace=True)
    Cytof_Data.initialize_cell_types(expression_matrix=expression_matrix,
                                     labels=labels,
                                     max_components=3,
                                     min_components=3,
                                     covariance_types=["diag"])

    Cytof_Data.generate_cell_graph(graph_topology="forest")
    temp = Cytof_Data.cell_graph.sample_graph(n_samples=1, cell_label=0)

    assert (nx.number_weakly_connected_components(Cytof_Data.cell_graph.graph) >= 1) and \
           (isinstance(temp[0], np.ndarray))


def test_cell_graph_error(Cytof_Data):
    expression_matrix = np.random.multivariate_normal(mean=np.zeros(30),
                                                      cov=np.eye(30),
                                                      size=1000)
    labels = np.random.choice(a=[0,1,2,3,4], size=1000, replace=True)
    Cytof_Data.initialize_cell_types(expression_matrix=expression_matrix,
                                     labels=labels,
                                     max_components=3,
                                     min_components=3,
                                     covariance_types=["diag"])
    try:
        Cytof_Data.generate_cell_graph(graph_topology="cyclic")
    except ValueError:
        assert True


def test_cell_graph_bead(Cytof_Data):
    expression_matrix = np.random.multivariate_normal(mean=np.zeros(30),
                                                      cov=np.eye(30),
                                                      size=1000)
    labels = np.random.choice(a=[0,1,2,3,4], size=1000, replace=True)
    Cytof_Data.bead_label = 0
    Cytof_Data.initialize_cell_types(expression_matrix=expression_matrix,
                                     labels=labels,
                                     max_components=1,
                                     min_components=1,
                                     covariance_types=["diag"])

    Cytof_Data.generate_cell_graph(graph_topology="forest")
    temp = Cytof_Data.cell_graph.sample_graph(n_samples=1, cell_label=0)

    assert (nx.number_weakly_connected_components(Cytof_Data.cell_graph.graph) >= 1) and \
           (isinstance(temp[0], np.ndarray))


def test_overall_batch_effects(Cytof_Data):
    expression_matrix = np.random.multivariate_normal(mean=np.zeros(30),
                                                      cov=np.eye(30),
                                                      size=1000)
    labels = np.random.choice(a=[0, 1, 2, 3, 4], size=1000, replace=True)
    Cytof_Data.initialize_cell_types(expression_matrix=expression_matrix,
                                     labels=labels,
                                     max_components=3,
                                     min_components=3,
                                     covariance_types=["diag"])
    Cytof_Data.generate_overall_batch_effects(variance=0.1)
    temp = Cytof_Data.sample(n_samples=1)

    assert (isinstance(temp[0], dict)) and \
           (np.sum(list(Cytof_Data.overall_batch_effects.values())) == 0)


def test_local_batch_effects(Cytof_Data):
    expression_matrix = np.random.multivariate_normal(mean=np.zeros(30),
                                                      cov=np.eye(30),
                                                      size=1000)
    labels = np.random.choice(a=[0, 1, 2, 3, 4], size=1000, replace=True)
    Cytof_Data.initialize_cell_types(expression_matrix=expression_matrix,
                                     labels=labels,
                                     max_components=3,
                                     min_components=3,
                                     covariance_types=["diag"])
    Cytof_Data.generate_local_batch_effects(variance=0.1)
    temp = Cytof_Data.sample(n_samples=1)

    assert (isinstance(temp[0], dict)) and \
           (np.sum(Cytof_Data.local_batch_effects[0][0]) == 0)


def test_temporal_effects(Cytof_Data):
    expression_matrix = np.random.multivariate_normal(mean=np.zeros(30),
                                                      cov=np.eye(30),
                                                      size=1000)
    labels = np.random.choice(a=[0, 1, 2, 3, 4], size=1000, replace=True)
    Cytof_Data.initialize_cell_types(expression_matrix=expression_matrix,
                                     labels=labels,
                                     max_components=3,
                                     min_components=3,
                                     covariance_types=["diag"])
    Cytof_Data.generate_temporal_effects(variance=0.1)
    temp = Cytof_Data.sample(n_samples=1)

    assert (isinstance(temp[0], dict)) and \
           (Cytof_Data.temporal_effects[0][0](0) == 0)


def test_cell_abundances(Cytof_Data):
    expression_matrix = np.random.multivariate_normal(mean=np.zeros(30),
                                                      cov=np.eye(30),
                                                      size=1000)
    labels = np.random.choice(a=[0, 1, 2, 3, 4], size=1000, replace=True)
    Cytof_Data.initialize_cell_types(expression_matrix=expression_matrix,
                                     labels=labels,
                                     max_components=3,
                                     min_components=3,
                                     covariance_types=["diag"])
    Cytof_Data.generate_cell_abundances(use_observed=True)

    assert np.sum(list(Cytof_Data.cell_abundances[0].values())) <= 1.5


def test_cell_abundances_random(Cytof_Data):
    expression_matrix = np.random.multivariate_normal(mean=np.zeros(30),
                                                      cov=np.eye(30),
                                                      size=1000)
    labels = np.random.choice(a=[0, 1, 2, 3, 4], size=1000, replace=True)
    Cytof_Data.initialize_cell_types(expression_matrix=expression_matrix,
                                     labels=labels,
                                     max_components=3,
                                     min_components=3,
                                     covariance_types=["diag"])
    Cytof_Data.generate_cell_abundances(use_observed=False,
                                        is_random=True)

    assert np.sum(list(Cytof_Data.cell_abundances[0].values())) <= 1.5


def test_overall_tree(Cytof_Data):
    expression_matrix = np.random.multivariate_normal(mean=np.zeros(30),
                                                      cov=np.eye(30),
                                                      size=1000)
    labels = np.random.choice(a=[0, 1, 2, 3, 4], size=1000, replace=True)
    Cytof_Data.initialize_cell_types(expression_matrix=expression_matrix,
                                     labels=labels,
                                     max_components=3,
                                     min_components=3,
                                     covariance_types=["diag"])
    Cytof_Data.generate_cell_graph(graph_topology="tree")
    Cytof_Data.generate_overall_batch_effects(variance=0.1)
    Cytof_Data.generate_local_batch_effects(variance=0.1)
    Cytof_Data.generate_temporal_effects(variance=0.1)
    Cytof_Data.generate_cell_abundances()
    temp = Cytof_Data.sample(n_samples=1)

    assert (isinstance(temp[3], dict))


def test_overall_forest(Cytof_Data):
    expression_matrix = np.random.multivariate_normal(mean=np.zeros(30),
                                                      cov=np.eye(30),
                                                      size=1000)
    labels = np.random.choice(a=[0, 1, 2, 3, 4], size=1000, replace=True)
    Cytof_Data.initialize_cell_types(expression_matrix=expression_matrix,
                                     labels=labels,
                                     max_components=3,
                                     min_components=3,
                                     covariance_types=["diag"])
    Cytof_Data.generate_cell_graph(graph_topology="forest")
    Cytof_Data.generate_overall_batch_effects(variance=0.1)
    Cytof_Data.generate_local_batch_effects(variance=0.1)
    Cytof_Data.generate_temporal_effects(variance=0.1)
    Cytof_Data.generate_cell_abundances()
    temp = Cytof_Data.sample(n_samples=1)

    assert (isinstance(temp[3], dict))
