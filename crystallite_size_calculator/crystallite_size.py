"""
A python module to estimate the average crystallite sizes
and crystal szie distritubtion computed from the
"""

from crystallite_size_calculator import envelope_function_approach, strain_calculator, xrd_calculator

from crystallite_size_calculator.xrd_calculator import PXRDProcessor
from crystallite_size_calculator import filetyper


class ComputeCrystalSizes:
    """
    A class to compute the crystal sizes and microstrains
    of crystal structures and powder diffraction patterns.
    """
    def __init__(self, cif_file=None, two_theta=None, intensity=None, wavelength=1.5406):
        self.cif_file = cif_file
        self.two_theta = two_theta
        self.intensity = intensity
        self.wavelength = wavelength
        if self.cif_file is not None:
            structure = filetyper.load_pystructure_from_cif(cif_file)
            two_theta, intensity, _ = xrd_calculator.get_pxrd_from_structure(structure)
            self.two_theta = two_theta
            self.intensity = intensity

    def compute_rdf_from_diffraction_pattern(self):
        """
        Compute the radial distribution function (RDF)
        from the diffraction pattern.

        **Parameters:**
            - two_theta : np.ndarray
                Array of 2θ angles.
            - intensities : np.ndarray
                Array of intensity values corresponding to each 2θ angle.
            - wavelength : float
                Wavelength of the X-ray source (in Angstroms).

        **Returns:**
            - r : np.ndarray
                Array of r values (atomic distances).
            - g_r : np.ndarray
                Array of g(r) values (radial distribution function).
        """
        pxrd_processor = PXRDProcessor(self.two_theta, self.intensity, self.wavelength)
        r, g_r = pxrd_processor.compute_pdf()
        return r, g_r

    def compute_crystallite_size_from_envelope_function(self):
        """
        Compute the average crystallite size (D)
        and crystallite size distribution (G(r))
        using the envelope function approach.

        **Returns:**
            - d_crys : float
                The estimated average crystallite size.
            - g_r : np.ndarray
                Array of G(r) values (radial distribution function).
        """
        r, g_r = self.compute_rdf_from_diffraction_pattern()
        d_crys = envelope_function_approach.fit_envelope_function(g_r, r, num_intervals=10, initial_guess=(1))
        return round(d_crys, 3)

    def size_strain_from_williamson_hall_method(self):
        """
        Compute the crystallite size and microstrain
        using the Williamson-Hall method.

        **Returns:**
            - strain : float
                The calculated microstrain.
            - d_crys : float
                The estimated average crystallite size.
        """
        strain, d_crys = strain_calculator.williamson_hall_method(self.two_theta, self.intensity, self.wavelength)
        return round(strain, 3), round(d_crys, 3)

    def compute_strain_from_williamson_hall_method(self):
        """
        Compute the crystallite size and microstrain
        using the Williamson-Hall method.

        **Returns:**
            - strain : float
                The calculated microstrain.
            - d_crys : float
                The estimated average crystallite size.
        """
        crystallite_size = self.compute_crystallite_size_from_envelope_function()
        fwhm_data, peak_positions = strain_calculator.estimate_fwhm_from_pxrd_no_profiling(self.two_theta, self.intensity)
        crsy = strain_calculator.strain_from_williamson_hall_method(fwhm_data, peak_positions, self.wavelength, crystallite_size)
        strain, d_crys = strain_calculator.williamson_hall_method(self.two_theta, self.intensity, self.wavelength)
        return  round(strain, 3), round(d_crys, 3)

    def size_and_strain_from_warren_averbach_method(self):
        """
        Compute the crystallite size and microstrain
        using the Warren-Averbach method.

        **Returns:**
            - strain : float
                The calculated microstrain.
            - d_crys : float
                The estimated average crystallite size.
        """
        crystallite_size = self.compute_crystallite_size_from_envelope_function()
        strain, d_crys = strain_calculator.warren_averbach_method(self.two_theta,
                                                                  self.intensity,
                                                                  self.wavelength,
                                                                  crystallite_size
                                                                  )
        return round(strain, 3), round(d_crys, 3)






