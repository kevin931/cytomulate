import cytomulate
import numpy as np

import pytest

from typing import List, Any, Union, Callable


class TestFileIO():
    
    @classmethod
    def setup_class(cls):
        existing_data: "np.ndarray" = np.random.rand(3000, 35)
        cell_type_indicator: List["np.ndarray"] = []
        cell_type_indicator.append(np.repeat(list(range(15)), 100))
        cell_type_indicator.append(np.repeat(list(range(15)), 200))
        cls.cytof_data_model: "cytomulate.cytof_data.CytofData"  = cytomulate.cytof_data.CytofData(n_batches=2,
                                                                                                   n_trees=5, 
                                                                                                   n_cells_per_batch=[1500, 3000],
                                                                                                   n_cell_types=15, 
                                                                                                   n_markers=35)
        cls.cytof_data_data: "cytomulate.cytof_data.CytofData"  = cytomulate.cytof_data.CytofData(n_batches=2,
                                                                                                  n_trees=5,
                                                                                                  n_cells_per_batch=[1500, 3000],
                                                                                                  n_cell_types=15,
                                                                                                  n_markers=35,
                                                                                                  expression_matrix=existing_data)
        
    
    @pytest.mark.parametrize("attr,value",
            [("n_markers", 35),
             ("n_batches", 2),
             ("n_trees", 5),
             ("n_cell_types", 15),
             ("n_cells_per_batch", [1500, 3000])]
            )
    def test_instance_attributes(self, attr: str, value: Any):
        assert getattr(self.cytof_data_model, attr) == value
        assert getattr(self.cytof_data_data, attr) == value
        
    
    @pytest.mark.parametrize("attr,var_type",
        [("forest", None),
         ("background_noise_level", None),
         ("cell_type_indicator", None)]
        )  
    def test_instance_attributes_none(self, attr: str, var_type: Any):
        assert getattr(self.cytof_data_model, attr) is var_type
        assert getattr(self.cytof_data_data, attr) is var_type
        
    
    @pytest.mark.parametrize("attr,var_type",
        [("cytof_data", dict),
         ("cell_type_proportions", np.ndarray),
         ("temporal_effect", list),
         ("batch_effect_by_cell_types_markers", dict),
         ("batch_effect_overall", np.ndarray),
         ("cell_types", dict)]
        )  
    def test_instance_attributes_type(self, attr: str, var_type: Any):
        assert isinstance(getattr(self.cytof_data_model, attr), var_type)
        assert isinstance(getattr(self.cytof_data_data, attr), var_type)
        
        
    def test_cytof_data_model(self):
        assert self.cytof_data_data.simulation_mode == "Data"
        assert isinstance(self.cytof_data_data.expression_matrix, np.ndarray)
    
    
    def test_cytof_data_data(self):
        assert self.cytof_data_model.simulation_mode == "Model"
        assert self.cytof_data_model.expression_matrix is None
        
        
    @pytest.mark.parametrize("batch", ["0", "1"])  
    def test_batch_effect_by_cell_types_markers(self, batch):
        batch_keys_model: List[str] = list(self.cytof_data_model.batch_effect_by_cell_types_markers.keys())
        batch_keys_data: List[str] = list(self.cytof_data_data.batch_effect_by_cell_types_markers.keys())
        batch_name: str = "batch" + batch
        assert batch_name in batch_keys_model
        assert batch_name in batch_keys_data
        assert isinstance(self.cytof_data_model.batch_effect_by_cell_types_markers[batch_name], np.ndarray)
        assert isinstance(self.cytof_data_data.batch_effect_by_cell_types_markers[batch_name], np.ndarray)
        
        
    @pytest.mark.parametrize("batch", ["0", "1"])  
    def test_cytof_data_batches(self, batch):
        batch_keys_model: List[str] = list(self.cytof_data_model.cytof_data.keys())
        batch_keys_data: List[str] = list(self.cytof_data_data.cytof_data.keys())
        batch_name: str = "batch" + batch
        assert batch_name in batch_keys_model
        assert batch_name in batch_keys_data
        assert isinstance(self.cytof_data_model.cytof_data[batch_name], dict)
        assert isinstance(self.cytof_data_data.cytof_data[batch_name], dict)
        assert isinstance(self.cytof_data_model.cytof_data[batch_name]["cell_type_indices"], np.ndarray)
        assert isinstance(self.cytof_data_data.cytof_data[batch_name]["cell_type_indices"], np.ndarray)
    
    
    def test_generate_tempral_effect(self):
        self.cytof_data_model.generate_temporal_effect()
        self.cytof_data_data.generate_temporal_effect()
        assert len(self.cytof_data_model.temporal_effect) > 0
        assert len(self.cytof_data_data.temporal_effect) > 0
        for func in self.cytof_data_model.temporal_effect:
            assert callable(func)
        for func in self.cytof_data_data.temporal_effect:
            assert callable(func)
    
    
    def test_generate_cell_type_proportions(self):
        self.cytof_data_model.generate_cell_type_proportions()
        self.cytof_data_data.generate_cell_type_proportions()
        assert isinstance(self.cytof_data_model.cell_type_proportions, np.ndarray)
        assert self.cytof_data_model.cell_type_proportions.shape == (2, 15)
        assert isinstance(self.cytof_data_data.cell_type_proportions, np.ndarray)
        assert self.cytof_data_data.cell_type_proportions.shape == (2, 15)
        
        
    def test_initialize_cell_types_model(self):
        self.cytof_data_model.initialize_cell_types()
        assert len(self.cytof_data_model.cell_types) == 15
        for i in range(15):
            assert isinstance(self.cytof_data_model.cell_types[i], cytomulate.cell_type.CellType)
            assert self.cytof_data_model.cell_types[i].id == i
            assert self.cytof_data_model.cell_types[i].name == i
            assert self.cytof_data_model.cell_types[i].n_markers == 35  
            
            
    def test_initialize_cell_types_data_clustering(self):
        self.cytof_data_data.initialize_cell_types()
        assert self.cytof_data_data.n_cell_types > 0
        for i in range(self.cytof_data_data.n_cell_types):
            assert isinstance(self.cytof_data_model.cell_types[i], cytomulate.cell_type.CellType)
            assert self.cytof_data_model.cell_types[i].id == i
            assert self.cytof_data_model.cell_types[i].name == i
            assert self.cytof_data_model.cell_types[i].n_markers == 35
            
            
    def test_initialize_cell_types_data_types_given(self):
        existing_data: "np.ndarray" = np.random.rand(500, 35)
        cell_types: List[str] = ["a", "b", "c", "d", "e"]
        cell_type_indicator: "np.ndarray" = np.repeat(cell_types, 100)
        cytof_data_types_given: "cytomulate.cytof_data.CytofData"  = cytomulate.cytof_data.CytofData(n_batches=1,
                                                                                                     n_trees=2,
                                                                                                     n_cells_per_batch=[500],
                                                                                                     n_cell_types=5,
                                                                                                     n_markers=35,
                                                                                                     expression_matrix=existing_data,
                                                                                                     cell_type_indicator=cell_type_indicator)
        cytof_data_types_given.initialize_cell_types()
        assert cytof_data_types_given.n_cell_types == 5
        for i in range(cytof_data_types_given.n_cell_types):
            assert isinstance(cytof_data_types_given.cell_types[i], cytomulate.cell_type.CellType)
            assert cytof_data_types_given.cell_types[i].id == i
            assert cytof_data_types_given.cell_types[i].name == cell_types[i]
            assert cytof_data_types_given.cell_types[i].n_markers == 35
            
            
    ## TODO: Initialize cell types fit_model method 
    ## TODO: Fix the Forest class. Now, this is broken.    
    # def test_grow_forest(self):
    #     self.cytof_data_model.grow_forest()
    #     self.cytof_data_data.grow_forest()
    #     assert isinstance(self.cytof_data_model.forest, cytomulate.forest.Forest)
    #     assert isinstance(self.cytof_data_data.forest, cytomulate.forest.Forest)
    
    
    
    @classmethod
    def teardown_class(cls):
        pass
