import pytest
from ase.io import read
import matplotlib.pyplot as plt
from crystallite_size_calculator import envelope_function_approach



@pytest.fixture
def mof_5():
    return read('tests/test_data/EDUSIF_fair.cif')

def test_compute_rdf(mof_5):
    '''
    Test the compute_rdf function with a sample ase atom object.
    '''
    r, g_r = envelope_function_approach.compute_rdf(mof_5, d=3)
    rmax, gmax = envelope_function_approach.find_maxima(r, g_r)
    d_crys = envelope_function_approach.fit_envelope_function(g_r, r, num_intervals=10, initial_guess=(10))
    print (d_crys)

    # plt.scatter(rmax, gmax)