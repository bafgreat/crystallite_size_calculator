# Crystallite Size Calculator

**crystallite-size-calculator** is a Python module designed to estimate
the average crystallite sizes and microstrain from powder X-ray diffraction (PXRD) data.
It provides methods based on the envelope function, Williamson-Hall method, Warren-Averbach method,
 and Scherrer equation to accurately compute crystallite sizes and microstrains.

## Features

- Compute the radial distribution function (RDF) from PXRD data.
- Estimate average crystallite size using the envelope function approach.
- Calculate crystallite size and microstrain using the Williamson-Hall method.
- Compute crystallite size using the Scherrer and modified Scherrer equations.
- Support for advanced strain and crystallite size analysis using the Warren-Averbach method.

## Installation

To install the **crystallite-size-calculator** module, use `pip`:

```bash
pip install crystallite_size_calculator
```

## Eaxamples
# Crystallite Size Calculator - CIF File Examples

Below are examples demonstrating how to use the **crystallite-size-calculator** module to load CIF files, compute crystallite sizes, and analyze microstrain.

### 1. Loading a CIF File and Computing Crystallite Size

You can load a crystal structure from a CIF file and compute the crystallite size using the envelope function approach.

```python
from crystallite_size_calculator.crystallite_size import ComputeCrystalSizes

# Provide the path to your CIF file
cif_file = "example_structure.cif"

# Initialize the ComputeCrystalSizes class with the CIF file
calculator = ComputeCrystalSizes(cif_file=cif_file)

# Compute the crystallite size using the envelope function approach
d_crys = calculator.compute_crystallite_size_from_envelope_function()
print(f"Crystallite size from CIF file: {d_crys:.2f} nm")

# Compute crystallite size and strain using the Williamson-Hall method
strain, d_crys = calculator.size_strain_from_williamson_hall_method()
print(f"Crystallite size: {d_crys:.2f} nm, Microstrain: {strain:.4f}")

# Compute crystallite size using the Scherrer equation
d_crys = calculator.size_from_scherrer_eq()
print(f"Crystallite size using Scherrer equation: {d_crys:.2f} nm")

# Compute crystallite size and strain using the Warren-Averbach method
strain, d_crys = calculator.size_and_strain_from_warren_averbach_method()
print(f"Crystallite size: {d_crys:.2f} nm, Microstrain: {strain:.4f}")


# Compute the radial distribution function (RDF) from the CIF file
r, g_r = calculator.compute_rdf_from_diffraction_pattern()
# Print the first few r and g(r) values
print("r values:", r[:5])
print("g(r) values:", g_r[:5])
```