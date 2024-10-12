import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from scipy.fft import fft
from ase.io import read
from pymatgen.analysis.diffraction.xrd import XRDCalculator
from pymatgen.analysis.diffraction.neutron import NDCalculator
from pymatgen.io.ase import AseAtomsAdaptor
from rdfpy import rdf



class PXRDProcessor:
    """
    A class to process Powder X-ray Diffraction (PXRD)
    data and compute the pair distribution function (PDF).

    **Attributes:**
        - wavelength: float
            Wavelength of the X-ray source (in Angstroms).
        - angles_2theta : np.ndarray
            Array of 2θ diffraction angles (in degrees).
        - intensity :np.ndarray
            Array of intensity values corresponding to each 2θ angle.
    """

    def __init__(self, angles_2theta, intensity, wavelength):
        """
        Initialize the PXRDProcessor with PXRD data
        and the X-ray source wavelength.

        **Parameters:**
            - angles_2theta : np.ndarray
                Array of 2θ diffraction angles (in degrees).
            - intensity :np.ndarray
                Array of intensity values corresponding to each 2θ angle.
            - wavelength :float
                Wavelength of the X-ray source (in Angstroms).
        """
        self.angles_2theta = angles_2theta
        self.intensity = intensity
        self.wavelength = wavelength

    def compute_q(self):
        """
        Convert 2θ angles (in degrees) to the scattering vector magnitude Q.

        **Returns:**
            - q: np.ndarray
                Array of scattering vector magnitudes Q.
        """
        theta = np.radians(self.angles_2theta / 2)
        q = (4 * np.pi / self.wavelength) * np.sin(theta)
        return q

    def interpolate_data(self, q):
        """
        Interpolate the intensity data to ensure evenly spaced
        Q values for Fourier transform.

        **Parameters:**
            - q : np.ndarray
                Array of scattering vector magnitudes Q.

        **Returns:**
            - q_uniform : np.ndarray
                Uniformly spaced Q values.
            - intensity_interpolated : np.ndarray
                Interpolated intensity values.
        """
        q_min, q_max = q.min(), q.max()
        q_uniform = np.linspace(q_min, q_max, len(q))
        intensity_interpolated = interp1d(q, self.intensity,
                                          kind='linear',
                                          fill_value="extrapolate"
                                          )
        return q_uniform, intensity_interpolated(q_uniform)

    def compute_structure_factor(self, intensity_uniform):
        """
        Compute the structure factor S(Q) by subtracting
        background and normalizing intensity.

        **Parameters:**
            - intensity_uniform :np.ndarray
                Array of interpolated intensity values.

        **Returns:**
            - s_q : np.ndarray
                The structure factor S(Q).
        """
        background = np.min(intensity_uniform)
        s_q = (intensity_uniform - background) / np.max(intensity_uniform - background)
        return s_q

    def fourier_transform(self, q_uniform, s_q):
        """
        Perform Fourier transform on the structure factor
        S(Q) to compute the radial distribution function g(r).

        **parameters:**
            - q_uniform :np.ndarray
                Uniformly spaced Q values.
            - s_q :np.ndarray)
                Structure factor S(Q).

        **Returns:**
            - r : np.ndarray
                Array of r values (atomic distances).
            - g_r :
                Array of g(r) values (radial distribution function).
        """
        delta_q = q_uniform[1] - q_uniform[0]
        g_r = fft(s_q * q_uniform) * delta_q
        r = np.fft.fftfreq(len(q_uniform), delta_q)
        return r, g_r.real

    def compute_pdf(self):
        """
        Compute the pair distribution function (PDF) from PXRD data.

        **Returns:**
            - r : np.ndarray
                Array of r values (atomic distances).
            - g_r :
                Array of g(r) values (radial distribution function).
        """
        q = self.compute_q()
        q_uniform, intensity_uniform = self.interpolate_data(q)
        s_q = self.compute_structure_factor(intensity_uniform)
        r, g_r = self.fourier_transform(q_uniform, s_q)
        return r, g_r


    def plot_pdf(self, r, g_r):
        """
        Plot the pair distribution function g(r) as a
        function of atomic distances r.

        **parameters:**
            - r : np.ndarray
                Array of r values (atomic distances).
            g_r : np.ndarray
                Array of g(r) values (radial distribution function).
        """
        plt.figure(figsize=(8, 5))
        plt.plot(r, g_r, label='Radial Distribution Function g(r)')
        plt.xlabel('r (Å)')
        plt.ylabel('g(r)')
        plt.title('Radial Distribution Function (PDF) from PXRD Data')
        plt.legend()
        plt.xlim(0, 10)
        plt.show()


def load_cif(cif_filename):
    """
    Load a crystal structure from a CIF file using ase and convert
    to pymatgen structure
    **Rarameters:**
        - cif_filename : str
        The filename of the CIF file.
    **Returns:**
        - pymatgen.Structure : The pymatgen structure.
    """
    ase_atoms = read(cif_filename)
    structure = AseAtomsAdaptor.get_structure(ase_atoms)
    return structure


def compute_pxrd(cif_filename, wavelengths='CuKa'):
    """
    Compute the powder X-ray diffraction

    **Parameters:**
        - cif_filename : str
            The filename of the CIF file.

    **Returns:**
        - r : np.ndarray
            Array of r values (atomic distances).
        - g_r : np.ndarray
            Array of g(r) values (radial distribution function).
    """
    structure = load_cif(cif_filename)
    xrd_calculator = XRDCalculator(wavelengths=wavelengths)
    pattern = xrd_calculator.get_pattern(structure)
    two_theta = np.array(pattern.x)
    intensities = np.array(pattern.y)
    hkl = pattern.hkls
    return two_theta, intensities, hkl


def compute_neutron_diffraction(cif_filename, wavelengths='CuKa'):
    """
    Compute the powder X-ray diffraction

    **Parameters:**
        - cif_filename : str
            The filename of the CIF file.

    **Returns:**
        - r : np.ndarray
            Array of r values (atomic distances).
        - g_r : np.ndarray
            Array of g(r) values (radial distribution function).
    """
    structure = load_cif(cif_filename)
    nd_calculator = NDCalculator(wavelengths=wavelengths)
    pattern = nd_calculator.get_pattern(structure)
    two_theta = np.array(pattern.x)
    intensities = np.array(pattern.y)
    hkl = pattern.hkls
    return two_theta, intensities, hkl


def compute_rdf_from_structure(ase_atoms, dr=0.05, d=5):
    """
    This function computes the radial distribution function, which
    describes the variation in  atomic density as a function of
    radial distance, r, within a material.

    The function takes as input an ase atom object and create a (d*d*d)
    supercell so that there are enough atoms in the crystal and then computes
    the radial distribution function for this supercell. Moreover, noise
    can added to the coordinates to make the resulting function smoother.

    **Parameters:**
        - ase_atoms : ase atom object
        - dr : float
            Determines the spacing between successive radii over
            which g(r) is computed.
        - d : int
            Determines the number of unit cells in each
            dimension in the supercell.
        - noise : bool
            If True, add noise to the coordinates.

    **Returns:**
        - g_r : (n_radii) np.array
            radial distribution function values g(r).
        -radii : (n_radii) np.array
            radii over which g(r) is computed
    **Example:**
        >>> g_r, radii = crystallite_size_calculator.\
            envelope_function_approach.\
                compute_rdf(coords, dr=0.05)
        >>> plt.plot(r, g_r)
    """
    supercell = ase_atoms * (d, d, d)
    coords = supercell.positions
    # noise = np.random.normal(loc=0.0, scale=0.05, size=(coords.shape))
    g_r, radii = rdf(coords, dr)
    return radii, g_r
