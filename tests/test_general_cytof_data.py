import pytest
import numpy as np
from cytomulate.creation.cytof_data import CreationCytofData
from cytomulate.cytof_data_general import GeneralCytofData
from cytomulate.utilities import univariate_noise_model

from typing import Dict

OPT_PCK: Dict[str, bool] = {"PyCytoData": True}

try:
    from PyCytoData import PyCytoData
except ImportError:
    OPT_PCK["PyCytoData"] = False

@pytest.fixture
def Cytof_Data():
    """Returns a CreationCytofData instance"""
    model = univariate_noise_model("normal", loc = 0, scale = 1)
    return GeneralCytofData(n_batches=2,
                            background_noise_model=model)

@pytest.fixture
def Cytof_Data1():
    """Returns a CreationCytofData instance"""
    model = {0: univariate_noise_model("normal", loc = 0, scale = 1),
             1: univariate_noise_model("normal", loc = 0, scale = 0.5)}
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


def test_generate_temporal_effects_spline_nondict(Cytof_Data):
    y = CreationCytofData()
    y.initialize_cell_types()
    Cytof_Data.cell_types = y.cell_types
    Cytof_Data.n_markers = y.n_markers
    Cytof_Data.cell_type_labels_to_ids = y.cell_type_labels_to_ids
    Cytof_Data.cell_type_ids_to_labels = y.cell_type_ids_to_labels
    Cytof_Data.generate_temporal_effects(x=np.linspace(0, 1, 10),
                                         y=np.zeros(10))
    assert isinstance(Cytof_Data.temporal_effects[0], list)


def test_sample_one_batch_probability(Cytof_Data1):
    y = CreationCytofData()
    y.initialize_cell_types()
    Cytof_Data1.cell_types = y.cell_types
    Cytof_Data1.n_markers = y.n_markers
    Cytof_Data1.cell_type_labels_to_ids = y.cell_type_labels_to_ids
    Cytof_Data1.cell_type_ids_to_labels = y.cell_type_ids_to_labels
    Cytof_Data1.generate_cell_abundances()
    temp = Cytof_Data1.sample_one_batch(1)
    assert (isinstance(temp[3], np.ndarray))


def test_sample_one_batch_counts(Cytof_Data1):
    y = CreationCytofData(n_types=3)
    y.initialize_cell_types()
    Cytof_Data1.cell_types = y.cell_types
    Cytof_Data1.n_markers = y.n_markers
    Cytof_Data1.cell_type_labels_to_ids = y.cell_type_labels_to_ids
    Cytof_Data1.cell_type_ids_to_labels = y.cell_type_ids_to_labels
    temp = Cytof_Data1.sample_one_batch(6, cell_abundances={0:1,1:2,2:3})
    assert (isinstance(temp[3], np.ndarray))


def test_sample_one_batch_error(Cytof_Data1):
    try:
        y = CreationCytofData(n_types=3)
        y.initialize_cell_types()
        Cytof_Data1.cell_types = y.cell_types
        Cytof_Data1.n_markers = y.n_markers
        Cytof_Data1.cell_type_labels_to_ids = y.cell_type_labels_to_ids
        Cytof_Data1.cell_type_ids_to_labels = y.cell_type_ids_to_labels
        temp = Cytof_Data1.sample_one_batch(6, cell_abundances={0:-1,1:2,2:3})
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
    y = CreationCytofData(n_types=3)
    y.initialize_cell_types()
    Cytof_Data.cell_types = y.cell_types
    Cytof_Data.n_markers = y.n_markers
    Cytof_Data.cell_type_labels_to_ids = y.cell_type_labels_to_ids
    Cytof_Data.cell_type_ids_to_labels = y.cell_type_ids_to_labels
    temp = Cytof_Data.sample(6, cell_abundances={0:1,1:2,2:3})
    assert (isinstance(temp[3], dict))


def test_sample_missing_abundance(Cytof_Data):
    y = CreationCytofData(n_types=3)
    y.initialize_cell_types()
    Cytof_Data.cell_types = y.cell_types
    Cytof_Data.n_markers = y.n_markers
    Cytof_Data.cell_type_labels_to_ids = y.cell_type_labels_to_ids
    Cytof_Data.cell_type_ids_to_labels = y.cell_type_ids_to_labels
    temp = Cytof_Data.sample(4, cell_abundances={0:1,2:3})
    assert temp[0][0].shape[0] == 4


@pytest.mark.parametrize("n_batches",[1,2])
def test_sample_to_pycytodata(n_batches: int):
    if OPT_PCK["PyCytoData"]:
        y: CreationCytofData = CreationCytofData(n_batches=n_batches)
        y.initialize_cell_types()
        pcd: PyCytoData = y.sample_to_pycytodata(100)
        assert isinstance(pcd, PyCytoData)
        assert pcd.n_cell_types <= 10
        assert pcd.n_samples == n_batches
        assert pcd.n_channels == 20
    