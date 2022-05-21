import pytest
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
    assert hasattr(Cytof_Data, "generate_cell_abundances")


def test_generate_overall_batch_effects(Cytof_Data):
    assert hasattr(Cytof_Data, "generate_overall_batch_effects")


def test_generate_local_batch_effects(Cytof_Data):
    assert hasattr(Cytof_Data, "generate_local_batch_effects")


def test_generate_temporal_effects(Cytof_Data):
    assert hasattr(Cytof_Data, "generate_temporal_effects")


def test_sample_one_batch(Cytof_Data):
    assert hasattr(Cytof_Data, "sample_one_batch")


def test_sample(Cytof_Data):
    assert hasattr(Cytof_Data, "sample")

