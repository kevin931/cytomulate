import pytest
import numpy as np
from cytomulate.creation.cytof_data import CreationCytofData
from cytomulate.cytof_data_general import GeneralCytofData
from cytomulate.utilities import univariate_noise_model


@pytest.fixture
def Cytof_Data():
    """Returns a CreationCytofData instance"""
    model = univariate_noise_model("normal", loc = 0, scale = 1)
    return GeneralCytofData(n_batches=2,
                            background_noise_model=model)


def test_n_markers(Cytof_Data):
    assert Cytof_Data.n_markers is None


def test_background_noise_model(Cytof_Data):
    temp = Cytof_Data.background_noise_model((5, 3))
    assert temp.shape == (5, 3)


def test_n_batches(Cytof_Data):
    assert Cytof_Data.n_batches == 2


def test_generate_cell_abundances(Cytof_Data):
    y = CreationCytofData()
    y.initialize_cell_types()
    Cytof_Data.cell_types = y.cell_types
    Cytof_Data.n_markers = y.n_markers
    Cytof_Data.generate_cell_abundances()
    assert np.sum(list(Cytof_Data.cell_abundances[0].values())) <= 1.5


def test_generate_cell_abundances_equal(Cytof_Data):
    y = CreationCytofData()
    y.initialize_cell_types()
    Cytof_Data.cell_types = y.cell_types
    Cytof_Data.n_markers = y.n_markers
    Cytof_Data.generate_cell_abundances(is_random=False)
    assert np.sum(list(Cytof_Data.cell_abundances[0].values())) <= 1.5


def test_generate_overall_batch_effects(Cytof_Data):
    y = CreationCytofData()
    y.initialize_cell_types()
    Cytof_Data.cell_types = y.cell_types
    Cytof_Data.n_markers = y.n_markers
    Cytof_Data.generate_overall_batch_effects(variance=0.1)

    assert (np.sum(list(Cytof_Data.overall_batch_effects.values())) == 0)


def test_generate_overall_batch_effects_one(Cytof_Data):
    y = CreationCytofData()
    y.initialize_cell_types()
    Cytof_Data.cell_types = y.cell_types
    Cytof_Data.n_markers = y.n_markers
    Cytof_Data.n_batches = 1
    Cytof_Data.generate_overall_batch_effects(variance=0.1)

    assert (np.sum(list(Cytof_Data.overall_batch_effects.values())) == 0)


def test_generate_local_batch_effects(Cytof_Data):
    y = CreationCytofData()
    y.initialize_cell_types()
    Cytof_Data.cell_types = y.cell_types
    Cytof_Data.n_markers = y.n_markers
    Cytof_Data.cell_type_labels_to_ids = y.cell_type_labels_to_ids
    Cytof_Data.cell_type_ids_to_labels = y.cell_type_ids_to_labels
    Cytof_Data.generate_local_batch_effects(variance=0.1)
    assert (np.sum(Cytof_Data.local_batch_effects[0][0]) == 0)


def test_generate_local_batch_effects_one(Cytof_Data):
    y = CreationCytofData()
    y.initialize_cell_types()
    Cytof_Data.cell_types = y.cell_types
    Cytof_Data.n_markers = y.n_markers
    Cytof_Data.cell_type_labels_to_ids = y.cell_type_labels_to_ids
    Cytof_Data.cell_type_ids_to_labels = y.cell_type_ids_to_labels
    Cytof_Data.n_batches = 1
    Cytof_Data.generate_local_batch_effects(variance=0.1)
    assert (np.sum(Cytof_Data.local_batch_effects[0][0]) == 0)


def test_generate_temporal_effects_brownian_bridge(Cytof_Data):
    y = CreationCytofData()
    y.initialize_cell_types()
    Cytof_Data.cell_types = y.cell_types
    Cytof_Data.n_markers = y.n_markers
    Cytof_Data.cell_type_labels_to_ids = y.cell_type_labels_to_ids
    Cytof_Data.cell_type_ids_to_labels = y.cell_type_ids_to_labels
    Cytof_Data.generate_temporal_effects(variance=0.1)
    assert Cytof_Data.temporal_effects[0][0](0) == 0


def test_generate_temporal_effects_polynomial(Cytof_Data):
    y = CreationCytofData()
    y.initialize_cell_types()
    Cytof_Data.cell_types = y.cell_types
    Cytof_Data.n_markers = y.n_markers
    Cytof_Data.cell_type_labels_to_ids = y.cell_type_labels_to_ids
    Cytof_Data.cell_type_ids_to_labels = y.cell_type_ids_to_labels
    Cytof_Data.generate_temporal_effects(variance=0.1, coefficients=[1,2,3])
    assert Cytof_Data.temporal_effects[0][0](0) == 0


def test_generate_temporal_effects_spline(Cytof_Data):
    y = CreationCytofData()
    y.initialize_cell_types()
    Cytof_Data.cell_types = y.cell_types
    Cytof_Data.n_markers = y.n_markers
    Cytof_Data.cell_type_labels_to_ids = y.cell_type_labels_to_ids
    Cytof_Data.cell_type_ids_to_labels = y.cell_type_ids_to_labels
    Cytof_Data.generate_temporal_effects(x={0:np.linspace(0, 1, 10),
                                            1:np.linspace(0, 1, 10)},
                                         y={0:np.zeros(10),
                                            1:np.zeros(10)})
    assert isinstance(Cytof_Data.temporal_effects[0], list)


def test_sample_one_batch_probability(Cytof_Data):
    y = CreationCytofData()
    y.initialize_cell_types()
    Cytof_Data.cell_types = y.cell_types
    Cytof_Data.n_markers = y.n_markers
    Cytof_Data.cell_type_labels_to_ids = y.cell_type_labels_to_ids
    Cytof_Data.cell_type_ids_to_labels = y.cell_type_ids_to_labels
    Cytof_Data.generate_cell_abundances()
    temp = Cytof_Data.sample_one_batch(1)
    assert (isinstance(temp[3], np.ndarray))


def test_sample_one_batch_counts(Cytof_Data):
    y = CreationCytofData(n_types=3)
    y.initialize_cell_types()
    Cytof_Data.cell_types = y.cell_types
    Cytof_Data.n_markers = y.n_markers
    Cytof_Data.cell_type_labels_to_ids = y.cell_type_labels_to_ids
    Cytof_Data.cell_type_ids_to_labels = y.cell_type_ids_to_labels
    temp = Cytof_Data.sample_one_batch(6, cell_abundances={0:1,1:2,2:3})
    assert (isinstance(temp[3], np.ndarray))


def test_sample_one_batch_error(Cytof_Data):
    try:
        y = CreationCytofData(n_types=3)
        y.initialize_cell_types()
        Cytof_Data.cell_types = y.cell_types
        Cytof_Data.n_markers = y.n_markers
        Cytof_Data.cell_type_labels_to_ids = y.cell_type_labels_to_ids
        Cytof_Data.cell_type_ids_to_labels = y.cell_type_ids_to_labels
        temp = Cytof_Data.sample_one_batch(6, cell_abundances={0:-1,1:2,2:3})
        assert (isinstance(temp[3], np.ndarray))
    except ValueError:
        assert True


def test_sample(Cytof_Data):
    y = CreationCytofData()
    y.initialize_cell_types()
    Cytof_Data.cell_types = y.cell_types
    Cytof_Data.n_markers = y.n_markers
    Cytof_Data.cell_type_labels_to_ids = y.cell_type_labels_to_ids
    Cytof_Data.cell_type_ids_to_labels = y.cell_type_ids_to_labels
    temp = Cytof_Data.sample(1)
    assert (isinstance(temp[3], dict))


def test_sample_copy_abundance(Cytof_Data):
    y = CreationCytofData()
    y.initialize_cell_types()
    Cytof_Data.cell_types = y.cell_types
    Cytof_Data.n_markers = y.n_markers
    Cytof_Data.cell_type_labels_to_ids = y.cell_type_labels_to_ids
    Cytof_Data.cell_type_ids_to_labels = y.cell_type_ids_to_labels
    temp = Cytof_Data.sample(6,cell_abundances={0:1,1:2,2:3})
    assert (isinstance(temp[3], dict))

