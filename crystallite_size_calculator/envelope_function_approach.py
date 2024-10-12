'''
The envelope function approach for accurately determining the crystallite sizes
and crystallite size distributions. This approach involves transforming
diffraction data into radial distribution functions (G(r)),
identifying maxima and fitting an envelope function to estimate
the average crystallite size (ACS). The module is implemented by following the
description found the manuscript publisehed by  Thorsten M. Gesinga and Lars Robben.
J. Appl. Cryst. (2024). 57, 1466–1476
'''

import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt


def find_maxima(r, g_r):
    """
    A function to identify all peaks in the pair distribution function G(r).
    Each peak in the G(r) function represents the most probable atomic
    distances within the structure. This function detects the
    positions of these peaks based on changes in the sign of the
    derivative of G(r).

    **parameter:**
        - r np.array: Array of r values representing radial distances.
        - g_r np.array: Array of G(r) values representing the pair
            distribution function.

    **Returns:**
        tuple: A tuple containing:
            - r_max np.array: Array of r values corresponding
                to the maxima in G(r).
            - g_max np.array: Array of G(r) values at the maxima.

    **Example:**
        >>> r_max, g_max = find_maxima(r, g_r)
        >>> plt.scatter(r_max, g_max)
    """
    g_r, r = remove_nan(g_r, r)
    maxima = (np.diff(np.sign(np.diff(g_r))) < 0).nonzero()[0] + 1
    return r[maxima], g_r[maxima]


def envelope_function(r, d_crys):
    """
    This function implements a general computation of
    the envelope function that is used to describe the
    decay of the G(r) function in order to deternine
    the average crystallite size.


    **Parameters:**
        - r (float): Radial distance.
        -d_crys (float): Diameter of the crystallite.
    **Returns:**
        float: Value of the envelope function at the given r.

    **Example:**
        >>> r_values = np.linspace(0, 100, 1000)
        >>> envelope = envelope_function(r_values, D=50)
        >>> plt.plot(r_values, envelope)
    """
    return (1 - 3/2 * r/d_crys + 1/2 * (r/d_crys)**3) * (r < d_crys)


def fit_envelope(r_max, g_max):
    """
    Funtion to fit the envelope function to the peaks of the G(r) data to
    estimate the average crystallite size. This function takes the peaks of
    G(r) and fits an envelope function to extract the average crystallite
    size, D. The envelope function represents the decay of atomic
    correlation beyond a certain distance, which is related to
    the crystallite size.

    **Parameters:**
        - r_max np.array: Array of r values
            corresponding to the maxima in G(r).
        - g_max np.array: Array of G(r) values at the maxima.

    **Returns:**
        - float: The estimated average crystallite size D in nanometers.

    **Example:**
        >>> D = fit_envelope(r_max, g_max)
        >>> print(f"Estimated crystallite size: {D:.2f} nm")
    """
    # Curve fitting the envelope function to the maxima
    lower_bounds = [0, 0]   # Lower bounds for [D, n]
    upper_bounds = [np.inf, 3]  # Upper bounds for [D, n]

    # Perform the curve fitting
    popt, _ = curve_fit(envelope_function, r_max, g_max,
                        bounds=(lower_bounds, upper_bounds)
                        )

    d_crys = popt[0]
    return d_crys


def fit_envelope_approximation(r, g_robs, interval_length=10):
    """
    Approximate the envelope function by extracting maxima from the G(r) data
    and fitting the envelope function.

    **Parameters:**
       - r np.array:
        Array of radial distances.
       - g_robs np.array:
        Observed G(r) values.
       - interval_length: (int, optional):
        The length of the intervals for
            maxima extraction. Defaults to 10.

    **Returns:**
        - d_crys :float
            The estimated average crystallite size (D).
    """
    # Normalize Gobs(r)
    g_rnorm = g_robs / max(g_robs)

    # Divide Gnorm into intervals and find the maxima in each interval
    num_intervals = len(g_rnorm) // interval_length
    r_max, g_max = [], []

    for i in range(num_intervals):
        start = i * interval_length
        end = (i + 1) * interval_length
        max_idx = np.argmax(g_rnorm[start:end])
        r_max.append(r[start + max_idx])
        g_max.append(g_rnorm[start + max_idx])
    # Fit the envelope function to the extracted maxima
    popt, _ = curve_fit(envelope_function, r_max, g_max,
                        bounds=([0, 0], [np.inf, 3])
                        )
    d_crys = popt[0]  # Extract the fitted crystallite size (D)

    return d_crys


def spherical_envelope_function(r, d_crys):
    """
    This function computes the spherical envelope function for a single
    crystallite. The spherical envelope function models the
    decay in correlations due to the finite size of the
    crystallites. The equation used here is derived from
    Howell et al. (2006) for a single spherical crystallite.

    **Parameters:**
        - r : np.array
            The radial distance (r) values, typically from G(r) function.
        - d_crys : float
            The diameter of the crystallites
            (in nm or other appropriate units).

    **Returns:**
        - np.array
            The value of the envelope function for each r value, describing
            the decay of correlations due to finite crystallite size.
    """
    return np.heaviside(d_crys - r, 0.5) * (1 - (3/2)*(r/d_crys) + (1/2)*(r/d_crys)**3)


def extract_maxima(g_r, radii, num_intervals=5):
    """
    Thi function extracts the maxima of the radial distribution
    function G(r) over specified intervals. This function divides
    the observed G(r) data into several intervals and finds the maximum
    value within each interval. These maxima approximate the envelope
    function, which is used to fit and determine crystallite size.

    **Parameters:**
        - g_r : np.array
            The observed values of the radial distribution function G(r),
            which is calculated from the diffraction data.
        - radii : np.array
            The corresponding radial distance values for G(r).
        - num_intervals : int
            The number of intervals to split the G(r) data into.
            A reasonable number depends on the resolution of
            G(r) and the sample data.

    **Returns:**
        - tuple:
        - max_r_values : np.array
            Radial distance values where maxima occur in each interval.
        - max_G_values : np.array
            The corresponding G(r) maxima for the intervals.
    """
    max_r_values = []
    max_g_values = []
    g_r, radii = remove_nan(g_r, radii)
    interval_length = len(radii) // num_intervals

    for i in range(num_intervals):
        interval_start = i * interval_length
        interval_end = (i + 1) * interval_length
        r_interval = radii[interval_start:interval_end]
        g_interval = g_r[interval_start:interval_end]

        max_idx = np.argmax(g_interval)
        max_r_values.append(r_interval[max_idx])
        max_g_values.append(g_interval[max_idx])
    return np.array(max_r_values), np.array(max_g_values)


def remove_nan(g_r, radii):
    """
    Removes entries from g_r that contain np.nan and ensures the corresponding
    entries in radii are removed to keep both arrays synchronized.

    **Parameters:**
        - g_r :np.array
            Array of G(r) values, may contain np.nan.
        - radii:  np.array
            Array of radial distances corresponding to g_r.

    **Returns:**
        - tuple: A tuple containing:
            - g_r_clean :np.array
                The cleaned array of G(r) values without np.nan.
            - radii_clean : np.array
                The cleaned array of radii values, matching the cleaned g_r.

    Raises:
        ValueError: If g_r and radii are not of the same length.

    Example:
        >>> g_r = np.array([1.2, np.nan, 3.4, 4.5, np.nan])
        >>> radii = np.array([0.5, 1.0, 1.5, 2.0, 2.5])
        >>> g_r_clean, radii_clean = remove_nan(g_r, radii)
        >>> print(g_r_clean)
        >>> print(radii_clean)
    """
    if len(g_r) != len(radii):
        raise ValueError("g_r and radii must have the same length.")

    valid_mask = ~np.isnan(g_r)

    g_r_clean = g_r[valid_mask]
    radii_clean = radii[valid_mask]
    return g_r_clean, radii_clean


def fit_envelope_function(g_r, radii, num_intervals, initial_guess=(5)):
    """
    This function is used to determine the apparent average crystallite size (ACS)
    based on the fitting of the envelope function to the observed G(r) data.
    It fits the spherical envelope function to the extracted
    maxima of G(r) to determine crystallite size. The function extracts the
    maxima of the radial distribution function G(r) by dividing G(r) into
    intervals, then uses those maxima to fit the spherical envelope function.
    The crystallite size (D) is determined by finding the optimal parameter
    that minimizes the difference between the observed
    maxima and the theoretical envelope function.

    **Parameters:**
        - g_r : np.array
            The observed radial distribution function G(r) values.
        - r_values : np.array
            The corresponding radial distance values.
        - num_intervals :int
            Number of intervals to divide G(r) for extracting maxima.
        - initial_guess : tuple
            Initial guess for the crystallite size D for curve fitting.

    **Returns:**
        - d_crys :float
            The optimized crystallite diameter D (average crystallite size).


    """
    # Extract maxima from G(r) function
    r_max, g_max = extract_maxima(g_r, radii, num_intervals)

    # Perform curve fitting to determine crystallite size D
    popt, pcov = curve_fit(spherical_envelope_function,
                           r_max,
                           g_max,
                           p0=initial_guess
                           )
    d_crys = popt[0]
    return d_crys
