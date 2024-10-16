import pytest
from ase.io import read
import matplotlib.pyplot as plt
from crystallite_size_calculator.crystallite_size import ComputeCrystalSizes
import numpy as np



@pytest.fixture
def mof_5():
    return 'tests/test_data/EDUSIF_fair.cif'

def test_compute_crystallite_sizes(mof_5):
    '''
    Test the compute_crystallite_sizes function with a sample ase atom object.
    '''
    size_engine = ComputeCrystalSizes(mof_5)
    r, g_r  = size_engine.compute_rdf_from_diffraction_pattern()
    assert len(r) == 185
    assert len(g_r) == 185
    assert size_engine.compute_crystallite_size_from_envelope_function() == 7.133
    assert size_engine.size_strain_from_williamson_hall_method() == (np.float64(-0.013), np.float64(17.876))
    assert size_engine.compute_strain_from_williamson_hall_method() == (np.float64(-0.013), np.float64(17.876))
    assert size_engine.size_and_strain_from_warren_averbach_method() == (np.float64(0.034), np.float64(6.382))
    assert size_engine.size_from_scherrer_eq().all() == np.array([0.36958657, 0.45323925]).all()
    assert size_engine.size_from_modified_scherrer_eq() == 6.871

