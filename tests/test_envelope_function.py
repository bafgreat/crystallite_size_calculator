import pytest
from ase.io import read
import matplotlib.pyplot as plt
from crystallite_size_calculator.crystallite_size import ComputeCrystalSizes



@pytest.fixture
def mof_5():
    return 'tests/test_data/EDUSIF_fair.cif'

# def test_compute_rdf(mof_5):
#     '''
#     Test the compute_rdf function with a sample ase atom object.
#     '''
#     r, g_r = envelope_function_approach.compute_rdf(mof_5, d=1)
#     rmax, gmax = envelope_function_approach.find_maxima(r, g_r)
#     d_crys = envelope_function_approach.fit_envelope_function(g_r, r, num_intervals=10, initial_guess=(10))
#     print (d_crys)

#     # plt.scatter(rmax, gmax)


def test_compute_crystallite_sizes(mof_5):
    '''
    Test the compute_crystallite_sizes function with a sample ase atom object.
    '''
    size_engine = ComputeCrystalSizes(mof_5)
    r, g_r  = size_engine.compute_rdf_from_diffraction_pattern()
    assert len(r) ==  185
    assert len(g_r) == 185
    print (size_engine.compute_crystallite_size_from_envelope_function())
    print (size_engine.size_strain_from_williamson_hall_method())
    print (size_engine.compute_strain_from_williamson_hall_method())
    print (size_engine.size_and_strain_from_warren_averbach_method())

