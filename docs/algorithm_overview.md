# Crystallite Size Determination Using the Envelope Function Approach

This document describes a step-by-step process for determining crystallite sizes and their distribution based on powder diffraction data using an envelope function approach. The algorithm is based on methods outlined by Howell et al. (2006) and further discussed in the context of pair distribution functions (PDFs). The approach avoids limitations inherent in the Scherrer equation and provides a more detailed crystallite size distribution.

## Algorithm Steps

### 1. **Collect Powder Diffraction Data**
The process begins with the collection of powder diffraction data using standard X-ray diffraction (XRD) techniques. This diffraction data provides insights into the atomic structure and crystallite sizes in the sample.

### 2. **Generate the Pair Distribution Function \( G(r) \)**
Once the diffraction data is collected, the next step is to compute the radial distribution function, or pair distribution function \( G(r) \), which describes the probability of finding an atom at a distance \( r \) from a reference atom. This is usually done using specialized software such as `pdfgetX3` or `PDFgui`.

### 3. **Approximate the Envelope Function from \( G(r) \)**
The crystallite size information is encoded in the envelope function \( f_{\text{env}}(r) \), which describes the decay of correlations due to the finite size of crystallites. The relationship between the observed \( G(r) \), the ideal infinite-crystal PDF \( G_{\infty}(r) \), and the envelope function is given by:
\[
G_{\text{obs}}(r) = G_{\infty}(r) \cdot f_{\text{env}}(r)
\]
The goal is to approximate the envelope function by extracting maxima from \( G(r) \).

### 4. **Extract Maxima from \( G(r) \)**
The observed \( G(r) \) function is divided into several intervals, and in each interval, the maximum value of \( G(r) \) is extracted. These maxima represent the envelope function at specific radial distances. By reducing the complexity of the data to its maxima, the envelope function can be more easily approximated.

### 5. **Fit the Envelope Function**
The extracted maxima are then fitted to a theoretical model of the envelope function for spherical crystallites. The envelope function is typically modeled as:
\[
f_{\text{env}}(r; D) = 1 - \frac{3}{2} \left( \frac{r}{D} \right) + \frac{1}{2} \left( \frac{r}{D} \right)^3
\]
where \( D \) is the average crystallite size (diameter). The fitting process involves minimizing the difference between the extracted maxima and this theoretical function to find the optimal crystallite size \( D \).

### 6. **Determine the Apparent Average Crystallite Size (ACS)**
The crystallite size \( D \), obtained through the fitting procedure, represents the apparent average crystallite size (ACS). This method provides a more robust estimation of crystallite size than traditional approaches like the Scherrer equation, especially for cases where the crystallite size distribution is broad.

### 7. **Verify the Results**
The crystallite size obtained through this method can be compared with results from other techniques, such as transmission electron microscopy (TEM) or Rietveld refinement, to verify the accuracy. Additionally, the results can be cross-checked with the Maximum Observable Apparent Crystallite Size (MOACS), which depends on the X-ray wavelength and instrument properties.

---

### Advantages of the Envelope Function Approach:
- **Accurate Crystallite Size Estimation**: The method avoids assumptions made by the Scherrer equation, providing more accurate estimates of crystallite size, especially when there's a broad distribution.
- **Size Distribution**: In addition to determining the average size, this approach can be extended to estimate the distribution of crystallite sizes.
- **Robust to Noise**: By focusing on the maxima of \( G(r) \), the method is less sensitive to noise in the diffraction data.

This approach allows for detailed crystallite size analysis and is particularly useful in cases where traditional diffraction methods fall short.
