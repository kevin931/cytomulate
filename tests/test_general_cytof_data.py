import pytest
import networkx as nx
import numpy as np
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
