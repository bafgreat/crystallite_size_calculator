"""
    A module to perform strain analysis using both the Williamson-Hall and Warren-Averbach methods on
    X-ray diffraction (XRD) data. It also provides functionality to estimate peak FWHM from PXRD data

    Methods
    -------
    williamson_hall_method(fwhm_data, theta_data, wavelength, crystallite_size):
        Calculates strain using the Williamson-Hall method and plots the Williamson-Hall plot.

    warren_averbach_method(fwhm_data, peak_order, crystallite_size):
        Calculates strain using the Warren-Averbach method and plots FWHM vs peak order.

    get_diffraction_pattern(structure, wavelength):
        Generates a diffraction pattern from a given crystal structure using Pymatgen's XRD calculator.

    estimate_fwhm_from_pxrd(two_theta, intensities, prominence=0.1, height_threshold=0.05):
        Estimates the full width at half maximum (FWHM) from powder XRD data.

    run_analysis(structure_file, crystallite_size):
        Runs the strain analysis using both methods (Williamson-Hall and Warren-Averbach) on a given structure file.

"""
import numpy as np
import matplotlib.pyplot as plt
from pymatgen.analysis.diffraction.xrd import XRDCalculator
from scipy.signal import find_peaks
from pymatgen.core.structure import Structure
from scipy.special import erf, wofz
from scipy.optimize import curve_fit



def williamson_hall_method(two_theta, intensities, wavelength):
    """
    A function to calculate crystallite size and strain using the Williamson-Hall method.

    The Williamson-Hall method decomposes peak broadening into contributions from
    crystallite size and lattice strain. It uses the Full Width at Half Maximum (FWHM)
    and Bragg angles of the diffraction peaks to estimate the crystallite size (D) and
    the lattice strain (ε).

    The basic form of the Williamson-Hall equation is:

        β_total * cos(θ) = (k * λ / D) + 4 * ε * sin(θ)

    Where:
        - β_total: The total peak broadening (in radians) after correcting for instrumental broadening.
        - θ: The Bragg angle (half the diffraction angle, in radians).
        - k: A shape factor, generally approximated to 0.9.
        - λ: The wavelength of the X-rays (in Ångströms).
        - D: The crystallite size (in nanometers).
        - ε: The microstrain in the lattice.

    The function generates a Williamson-Hall plot (β_total * cos(θ) vs. 4 * sin(θ)),
    performs a linear regression on the data points, and extracts the slope (proportional to strain)
    and the intercept (related to crystallite size).

    Parameters:
    ----------
    two_theta : np.array
        Array of 2-theta values (in degrees) for the PXRD pattern.
    intensities : np.array
        Array of intensity values corresponding to the 2-theta values.
    wavelength : float
        X-ray wavelength in Angstroms (e.g., 1.5406 Å for Cu K-alpha).

    Returns:
    -------
    strain : float
        The calculated strain value (ε).
    crystallite_size : float
        The calculated crystallite size in nanometers (D).
    """

    # fwhm_data, peak_positions = estimate_fwhm_from_pxrd(two_theta, intensities)
    fwhm_data, peak_positions = estimate_fwhm_from_pxrd_no_profiling(two_theta, intensities)


    fwhm_radians = np.radians(fwhm_data)

    theta_data = peak_positions / 2
    theta_radians = np.radians(theta_data)

    cos_theta = np.cos(theta_radians)
    sin_theta = np.sin(theta_radians)

    y = fwhm_radians * cos_theta

    x = 4 * sin_theta

    coeffs = np.polyfit(x, y, 1)
    slope = coeffs[0]
    intercept = coeffs[1]

    strain = slope / 4
    k = 0.9
    crystallite_size = (k * wavelength) / intercept


    # plt.figure(figsize=(8, 6))
    # plt.plot(x, y, 'bo', label='Data')
    # plt.plot(x, np.polyval(coeffs, x), 'r-', label=f'Fit: Slope (Strain) = {strain:.4e}, Intercept (Crystallite Size) = {crystallite_size:.4f} nm')
    # plt.xlabel(r'$4 \sin(\theta)$')
    # plt.ylabel(r'$\beta \cos(\theta)$')
    # plt.title('Williamson-Hall Plot')
    # plt.legend()
    # plt.show()
    return strain, crystallite_size


def strain_from_williamson_hall_method(fwhm_data, theta_data, wavelength, crystallite_size):
    """
    A function to calculate strain using the Williamson-Hall method.
    This method decomposes the broadening of diffraction peaks into
    contributions from crystallite size and lattice strain. The broadening
    is characterized using the Full Width at Half Maximum (FWHM) of
    the diffraction peaks.

    The basic form of the Williamson-Hall equation is:

        β_total * cos(θ) = (k * λ / D) + 4 * ε * sin(θ)

    Where:
        - β_total: The total peak broadening (in radians) after correcting for
            instrumental broadening.
        - θ: The Bragg angle (half the diffraction angle, in radians).
        - k: A shape factor, genarally approxiamted to 0.9 for all systems.
        - λ: The wavelength of the X-rays (in nm or Ångstroms).
        - D: The crystallite size (in nm or Ångstroms)
        - ε: The microstrain in the lattice, also estimated by this method.

    The method involves plotting β_total * cos(θ) on the y-axis and 4 * sin(θ)
    on the x-axis, and performing a linear regression on the data points.
    The slope of the line is proportional to the lattice strain (ε),
    while the y-intercept is related to the crystallite size (D).

    **Parameters:**
        - fwhm_data : np.array
            Full-width at half maximum (in radians) for each peak.
        - theta_data : np.array
            Bragg angles (in degrees) for each peak.
        - wavelength : float
            X-ray wavelength in Angstroms.
        - crystallite_size : float
            Crystallite size in nanometers.

    **Returns:**
        - strain : float
            The calculated strain value.
    """
    theta_radians = np.radians(theta_data)
    cos_theta = np.cos(theta_radians)
    sin_theta = np.sin(theta_radians)

    size_term = (0.9 * wavelength) / (crystallite_size * cos_theta)
    strain_term = fwhm_data * cos_theta - size_term

    x = 4 * sin_theta
    y = fwhm_data * cos_theta

    coeffs = np.polyfit(x, y, 1)
    strain = coeffs[0]

    # plt.figure(figsize=(8, 6))
    # plt.plot(x, y, 'bo', label='Data')
    # plt.plot(x, np.polyval(coeffs, x), 'r-', label=f'Fit: Slope (Strain) = {strain:.4e}')
    # plt.xlabel(r'$4 \sin(\theta)$')
    # plt.ylabel(r'$\beta \cos(\theta)$')
    # plt.title('Williamson-Hall Plot')
    # plt.legend()
    # plt.show()

    return strain


def voigt_profile(x, amplitude, center, sigma, gamma):
    """
    Defines a Voigt profile for fitting peaks in X-ray diffraction data.
    The Voigt profile is a convolution of a Lorentzian and a Gaussian function,
    often used to model diffraction peaks that exhibit both Lorentzian
    (lifetime broadening) and Gaussian (instrumental broadening) components.

    **Parameters:**
        - x : np.array
            Array of x-values (typically 2-theta values in degrees)
            over which the Voigt profile is calculated.
        - amplitude : float
            The amplitude (height) of the peak.
        - center : float
            The center of the peak (in degrees 2-theta).
        - sigma : float
            The width of the Gaussian component of the Voigt
            profile (related to instrumental broadening).
        - gamma : float
            The width of the Lorentzian component of the
            Voigt profile (related to sample effects such as strain or defects).

    **Returns:**
        - np.array
            Array of y-values representing the Voigt profile
            at each x-value. This can be used to model or fit diffraction peaks.
    """
    z = ((x - center) + 1j * gamma) / (sigma * np.sqrt(2))
    return amplitude * np.real(wofz(z)) / (sigma * np.sqrt(2 * np.pi))


def estimate_fwhm_from_pxrd(two_theta, intensities, prominence=0.1, height_threshold=0.05):
    """
    A function to estimate the Full Width at Half Maximum (FWHM)
    and peak positions from PXRD data. This function identifies
    peaks using a peak detection algorithm and fits the peaks with
    a Voigt profile to accurately calculate the FWHM for each peak.
    The function estimates the initial guesses for the Gaussian
    width (sigma) and Lorentzian width (gamma) based on the data.
    To achive this the initial estimates for `sigma` (Gaussian width)
    are derived from the peak's Full Width at Half Maximum (FWHM).
    Moreover the Gamma (Lorentzian width) is initially
    estimated as a fraction of the FWHM.

    **Parameters:**
        - two_theta : np.array
            Array of 2-theta values (in degrees), representing the angles at which
            diffraction intensities were measured.
        - intensities : np.array
            Array of intensity values corresponding to the 2-theta values, representing the PXRD pattern.
        - prominence : float, optional
            The minimum prominence of peaks to be considered. This parameter helps differentiate significant
            peaks from background noise. Default is 0.1.
        - height_threshold : float, optional
            The minimum height of peaks to be considered, as a fraction of the maximum intensity in the pattern.
            Default is 0.05 (i.e., 5% of the maximum intensity).

    **Returns:**
        - fwhm_data : np.array
            Array of FWHM values (in degrees 2-theta) for each detected peak.
        - peak_positions : np.array
            Array of 2-theta positions (in degrees) for each detected peak.
    """
    peaks, properties = find_peaks(intensities, prominence=prominence, height=height_threshold * np.max(intensities))
    fwhm_data = []
    peak_positions = []

    for i, peak in enumerate(peaks):
        left_base = properties["left_bases"][i]
        right_base = properties["right_bases"][i]
        peak_region_x = two_theta[left_base:right_base]
        peak_region_y = intensities[left_base:right_base]

        amplitude = peak_region_y.max()
        center = two_theta[peak]

        # Estimate the FWHM by finding points near half max intensity
        half_max = amplitude / 2
        closest_to_half_max = np.where(np.isclose(peak_region_y, half_max, atol=0.1 * half_max))[0]

        if len(closest_to_half_max) >= 2:
            estimated_fwhm = peak_region_x[closest_to_half_max[-1]] - peak_region_x[closest_to_half_max[0]]
            sigma = estimated_fwhm / (2 * np.sqrt(2 * np.log(2)))  # FWHM to sigma conversion
        else:
            estimated_fwhm = (peak_region_x[-1] - peak_region_x[0]) / 2  # Use full width as an estimate
            sigma = estimated_fwhm / (2 * np.sqrt(2 * np.log(2)))

        # Gamma as a fraction of FWHM
        gamma = estimated_fwhm / 2 if 'estimated_fwhm' in locals() else 0.1

        try:
            # Fit the Voigt profile using 'trf' method
            popt, _ = curve_fit(
                voigt_profile,
                peak_region_x,
                peak_region_y,
                p0=[amplitude, center, sigma, gamma],
                method='trf',
                max_nfev=5000,  # Increase max_nfev
                bounds=([0, center - 2, 0, 0], [np.inf, center + 2, np.inf, np.inf])
            )
            fitted_sigma = popt[2]
            fitted_gamma = popt[3]

            # Calculate FWHM using Voigt parameters
            fwhm = 0.5346 * (2 * fitted_gamma) + np.sqrt(0.2166 * (2 * fitted_gamma)**2 + (2 * fitted_sigma)**2)
        except RuntimeError:
            # Fallback to estimated FWHM if fitting fails
            fwhm = estimated_fwhm

        fwhm_data.append(fwhm)
        peak_positions.append(center)

    return np.array(fwhm_data), np.array(peak_positions)

def warren_averbach_method(two_theta,
                           intensities,
                           crystallite_size,
                           wavelength,
                           instrumental_fwhm=None,
                           prominence=0.1,
                           height_threshold=0.05):
    """
    A function to calculate strain using the Warren-Averbach method.
    This method is used to separate the contributions of crystallite
    size and microstrain to the broadening of X-ray diffraction peaks.
    The analysis is performed by applying Fourier transforms to the
    diffraction peak profiles, allowing for a more detailed description
    of strain distribution.

    The method assumes that the broadening of diffraction peaks is
    influenced by two factors:
        1. Crystallite size (broadening due to finite crystallite size).
        2. Microstrain (broadening due to strain distribution within
        the crystallites).

    The Fourier coefficients A(L), which describe the peak shape as a
    function of the Fourier length L, are given by:
        A(L) = exp(-L/D) * exp(-2 * pi^2 * <epsilon^2(L)> * L^2)

    Where:
        - A(L): The Fourier coefficient at distance L.
        - D: The average crystallite size (in the direction of analysis).
        - <epsilon^2(L)>: The mean square strain as a function of distance L.
        - L: The Fourier length (distance over which the strain is analyzed).

    The crystallite size contribution to the Fourier coefficients is:
        A_s(L) = exp(-L/D)

    The strain contribution to the Fourier coefficients is:
        A_e(L) = exp(-2 * pi^2 * <epsilon^2(L)> * L^2)

    The total broadening of the diffraction peaks can be described as:
        beta_total^2 = beta_size^2 + beta_strain^2
    Where:
        - beta_total: The total peak broadening (in radians).
        - beta_size: The broadening due to crystallite size.
        - beta_strain: The broadening due to microstrain distribution.

    This function first fits the diffraction peaks using a Voigt profile
    (a convolution of Lorentzian and Gaussian profiles), calculates the
    Full Width at Half Maximum (FWHM), and subtracts the instrumental
    broadening. The Fourier coefficients are calculated from the corrected
    FWHM values, and a logarithmic fit is used to separate the contributions
    of crystallite size and strain.
    **Steps:**
        1. **Peak Detection and Fitting:**
           The function uses a Voigt profile to fit the diffraction peaks in
           the PXRD data. This fitting helps to accurately extract the
           Full Width at Half Maximum (FWHM) for each peak.

        2. **Instrumental Broadening Correction:**
           The instrumental broadening is subtracted from the measured FWHM to
           isolate the broadening due to crystallite size and strain.
           The corrected FWHM values are converted from degrees to radians.
           N.B
           From a theoretical approach, one can estimate the FWHM values the
           the pxrd of LaB6 or Si as standard. Simply use pymatgen to
           determine the PXRD of either of these systems and the compute their
           FWHM and use as instrumental_fwhm.

        3. **Fourier Analysis:**
           The interplanar spacing (d) is calculated using Bragg’s Law:
               d = λ / (2 * sin(θ))
           where θ is half of the 2-theta value. The Fourier coefficients
           A(L) are then calculated based on the corrected FWHM and
           crystallite size.

        4. **Logarithmic Fit:**
           The function performs a logarithmic fit on the Fourier coefficients:
               ln(A(L)) = -L/D - 2 * pi^2 * <epsilon^2> * L^2
           From this fit, the crystallite size (D) and the strain
           (<epsilon^2>) are extracted.

        5. **Strain and Crystallite Size Calculation:**
           The strain is calculated as the square root of the second-order
           coefficient from the logarithmic fit:
               strain = sqrt(-coefficient[0] / (2 * pi^2))
           The crystallite size is calculated from the first-order
           coefficient of the fit.

    **Parameters:**
        - two_theta : np.array
            Array of 2-theta values (in degrees), representing the angles
            at which diffraction intensities are measured.
        - intensities : np.array
            Array of intensity values corresponding to the 2-theta values,
            representing the PXRD pattern.
        - crystallite_size : float
            An initial guess for the crystallite size (in nanometers).
            The final crystallite size will be calculated from the fit.
        - wavelength : float
            The wavelength of the X-rays used in the experiment (in Ångströms).
            Common values are around 1.5406 Å for Cu K-alpha radiation.
        - instrumental_fwhm : numpy array
            The instrumental broadening (in degrees), which must be subtracted
            from the measured FWHM to isolate the sample's broadening.
            This can be obtained from a standard sample (e.g., Si or LaB6).
            if instrumental_fwhm is not None, default value of
            [26.22637089, 36.93188358,  9.7547036 ] is used computed from LaB6.
            with a wavelength of 1.5406
        - prominence : float, optional
            Minimum prominence of peaks to be considered in the analysis
            (relative to surrounding noise). Default is 0.1.
        - height_threshold : float, optional
            Minimum height of peaks to be considered, as a fraction of
            the maximum intensity. Default is 0.05.

    **Returns:**
        - strain : float
            The calculated strain value, representing the microstrain
            distribution in the crystallites.
        - d_crys : float
            The calculated crystallite size, refined from the Fourier
            analysis of the diffraction peaks.
    """
    # if instrumental_fwhm is None:
    #     instrumental_fwhm = np.array([26.22637089, 36.93188358,  9.7547036 ])

    fwhm_data, peak_positions = estimate_fwhm_from_pxrd_no_profiling(two_theta,
                                                        intensities,
                                                        prominence,
                                                        height_threshold)

    # fwhm_corrected = np.sqrt(fwhm_data**2 - instrumental_fwhm**2)
    # we will not use the corrected instrumental_fwhm

    # fwhm_data_radians = np.radians(fwhm_corrected)
    fwhm_value = np.sqrt(fwhm_data)
    fwhm_data_radians = np.radians(fwhm_value)
    theta = np.radians(peak_positions / 2)
    d_spacing = wavelength / (2 * np.sin(theta))
    l_values = d_spacing
    beta_total = fwhm_data_radians
    a_l = np.exp(-l_values / crystallite_size) *\
        np.exp(-2 * np.pi**2 * (beta_total**2) * l_values**2)

    # Logarithmic fit to separate size and strain contributions
    log_a_l = np.log(a_l)
    # Linear fit of the form: ln(A(L)) = -L/D - 2 * pi^2 * <epsilon^2> * L^2
    coefficients = np.polyfit(l_values, log_a_l, 2)

    d_crys = -1 / coefficients[1]
    strain = np.sqrt(-coefficients[0] / (2 * np.pi**2))

    # Plot the Fourier coefficients and the fit
    # plt.figure(figsize=(8, 6))
    # plt.plot(l_values, a_l, 'bo', label='Fourier Coefficients A(L)')
    # plt.plot(l_values, np.exp(np.polyval(coefficients, l_values)), 'r--', label='Fit')
    # plt.xlabel('Fourier Length (L)')
    # plt.ylabel('Fourier Coefficient A(L)')
    # plt.title('Warren-Averbach Analysis')
    # plt.legend()
    # plt.show()
    return strain, d_crys


def estimate_fwhm_from_pxrd_no_profiling(two_theta, intensities, prominence=0.1, height_threshold=0.05):
    """
    Estimates the full width at half maximum (FWHM) for peaks in PXRD data.

    **Parameters:**
        -  two_theta : array-like
            Array of 2-theta values (in degrees).
        - intensities : array-like
            Array of intensity values corresponding to the 2-theta values.
        - prominence : float, optional
            Minimum prominence of peaks to be considered. Default is 0.1.
        - height_threshold : float, optional
            Minimum height of peaks to be considered (relative to the maximum intensity). Default is 0.05.

    **Returns:**
        - fwhms : list of floats
            List of estimated FWHMs for each detected peak.
        - peak_positions : list of floats
            List of 2-theta positions of the detected peaks.
    """
    normalized_intensities = intensities / np.max(intensities)
    peak_indices, peak_properties = find_peaks(
        normalized_intensities, height=height_threshold, prominence=prominence)

    fwhms = []
    peak_positions = []

    for peak_index in peak_indices:
        half_max = intensities[peak_index] / 2.0

        left_indices = np.where(intensities[:peak_index] < half_max)[0]
        right_indices = np.where(intensities[peak_index:] < half_max)[0]

        if len(left_indices) == 0 or len(right_indices) == 0:
            continue

        left_idx = left_indices[-1]
        right_idx = right_indices[0] + peak_index

        fwhm = two_theta[right_idx] - two_theta[left_idx]
        fwhms.append(fwhm)
        peak_positions.append(two_theta[peak_index])


    return np.array(fwhms), np.array(peak_positions)

