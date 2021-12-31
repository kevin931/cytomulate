import cytomulate
import numpy as np
from sklearn.mixture import GaussianMixture

import pytest

from typing import List, Any, Callable


class TestCellType():
    
    @classmethod
    def setup_class(cls):
        cls.cell: cytomulate.cell_type.CellType = cytomulate.cell_type.CellType(2, "a_2", 35)
        cls.parent: cytomulate.cell_type.CellType = cytomulate.cell_type.CellType(1, "a_1", 35)
        expression_matix: "np.ndarray" = np.random.rand(100, 35)
        cls.parent.fit_model(dat = expression_matix)
        cls.parent.markers_pattern = np.random.choice([0,1], size = 35).reshape(1,35)
        cls.parent.gating_markers = {2,3}
        cls.child: cytomulate.cell_type.CellType = cytomulate.cell_type.CellType(3, "a_3", 35)
        cls.child.fit_model(dat = expression_matix)
        
        
    @pytest.mark.parametrize("attr,value",
            [("id", 2),
             ("name", "a_2"),
             ("n_markers", 35)]
            )
    def test_instance_attributes_value(self, attr: str, value: Any):
        assert getattr(self.cell, attr) == value
        
        
    @pytest.mark.parametrize("attr",
        ["parent",
         "model_for_expressed_markers",
         "background_noise_level"]) 
    def test_instance_attribute_none(self, attr: str):
        assert getattr(self.cell, attr) is None
        
        
    @pytest.mark.parametrize("attr,type",
        [("children", list),
         ("marker_pattern", np.ndarray),
         ("expression_level", np.ndarray),
         ("variance_level", np.ndarray),
         ("parent_cell_type", dict),
         ("children_cell_types", dict),
         ("differential_paths_to_children", dict),
         ("markers_pattern", np.ndarray),
         ("gating_markers", set),
         ("overall_mean", np.ndarray),
         ("overall_var", np.ndarray)]) 
    def test_instance_attribute_type(self, attr: str, type: Any):
        assert isinstance(getattr(self.cell, attr), type)
        
    
    @pytest.mark.parametrize("method,expected",
        [("inherit_markers_pattern", "This is the root."),
         ("inherit_model", "This is the root.")]) 
    def test_parent_cell_type_exception(self, method: str, expected: str):
        try:
            getattr(self.cell, method)()
        except Exception as e:
            assert expected in str(e)
        else:
            raise
            
            
    def test_fit_model(self):
        expression_matix: "np.ndarray" = np.random.rand(500, 35)
        self.cell.fit_model(dat = expression_matix)
        assert isinstance(self.cell.model_for_expressed_markers, GaussianMixture)
        assert isinstance(self.cell.overall_mean, np.ndarray)
        assert self.cell.overall_var is None
            
            
    def test_generate_initial_expression(self):
        expression = self.cell.generate_initial_expressions()
        assert isinstance(expression, np.ndarray)
        assert expression.shape == (35,)
        assert expression.dtype == np.dtype("float64")
        
        
    def test_inherit_model(self):
        self.cell.parent_cell_type[1] = self.parent
        self.cell.inherit_model()
        assert isinstance(self.cell.model_for_expressed_markers, GaussianMixture)
        assert isinstance(self.parent.model_for_expressed_markers, GaussianMixture)
        assert np.array_equal(self.cell.model_for_expressed_markers.means_, self.parent.model_for_expressed_markers.means_)
    
    
    @pytest.mark.parametrize("mutation,additional_markers",
        [(0.2, 2),
         (0, 2),
         (1, 2),
         (0.2, 0),
         (0.2, 33)]) 
    def test_inherit_markers_pattern(self, mutation: float, additional_markers: int):
        self.cell.inherit_markers_pattern(mutation_probability=mutation, n_additional_gating_markers=additional_markers)
        for marker in {2, 3}:  
            assert marker in self.cell.gating_markers
        assert len(self.cell.gating_markers) == 2 + additional_markers
        assert np.sum(np.not_equal(self.cell.markers_pattern, self.parent.markers_pattern)) == int(33*mutation)
        
        
    @pytest.mark.parametrize("mutation,additional_markers,expected",
        [(1.1, 2, "Mutation probability must be between 0 and 1."),
         (-0.1, 2, "Mutation probability must be between 0 and 1."),
         (0.2, -1, "Marker number has to be positive."),
         (0.2, 100, "Maker number exceeds the number of remaining markers.")]) 
    def test_inherit_markers_patten_error(self, mutation: float, additional_markers: int, expected: str):
        try:
            self.cell.inherit_markers_pattern(mutation_probability=mutation, n_additional_gating_markers=additional_markers)
        except ValueError as e:
            assert str(e) == expected
        else:
            raise
            
    
    @pytest.mark.parametrize("alpha,beta",
        [(0.4, 1),
         (1, 1.2)])         
    def test_generate_from_path(self, alpha: float, beta: float):
        self.cell.children.append(self.child)
        self.cell.children_cell_types[3] = self.child
        
        differentiating_path: List[Callable] = []
        for i in range(35):
            original_value: float = self.cell.markers_pattern[0][i]
            end_value: float = 0
            if np.isclose(original_value, 0):
                end_value = 1
            differentiating_path.append(cytomulate.utilities.smooth_brownian_bridge(original_value, end_value))
            
        self.cell.differential_paths_to_children[3] = differentiating_path
        differentiating_point: np.ndarray = self.cell.generate_from_paths(3, alpha=alpha, beta=beta)
        assert isinstance(differentiating_point, np.ndarray)
        assert differentiating_point.shape == (1, 35)
        
        
    @pytest.mark.parametrize("differentiate,alpha,beta",
        [(True, 0.4, 1),
         (True, 1, 1.2),
         (False, 0.4, 1)])     
    def test_generate_final_expression(self, differentiate: bool, alpha: float, beta: float):
        final_expressions: np.ndarray = self.cell.generate_final_expressions(differentiate=differentiate, alpha=alpha, beta=beta)
        assert isinstance(final_expressions, np.ndarray)
        assert final_expressions.shape == (35, )
        assert final_expressions.dtype == np.dtype("float64")
        
        
    @pytest.mark.parametrize("alpha,beta",
        [(-0.4, 1),
         (1, -0.4)])     
    def test_beta_values_error(self, alpha: float, beta: float):
        try:
            self.cell.generate_final_expressions(alpha=alpha, beta=beta)
        except ValueError:
            assert True
        else:
            raise
    
    
    @classmethod
    def teardown_class(cls):
        pass
