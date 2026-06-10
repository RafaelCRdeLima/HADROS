# UHE Spectral Models

The UHE source is now factored into two independent pieces:

```text
emissivity_UHE(r, theta, E) = source_morphology(r, theta) * spectral_weight(E)
```

The transport pipeline interacts only with the generic spectral interface:

```text
uhe_spectral_weight(E_GeV, source_params, spectral_params)
```

This keeps the spectral model separate from the Kerr geodesics, DIS opacity,
density backgrounds, and source morphology definitions.

## Implemented Models

### Monochromatic

```text
SPECTRAL_MODEL=monochromatic
```

This is the backward-compatible single-energy mode. It evaluates the existing
ray-tracing pipeline at the selected `ENU` without integrating over a broad
energy distribution.

For exact continuity with previous image normalization, the legacy source
energy factor controlled by `SOURCE_POWERLAW` and `SOURCE_EMAX_GEV` is retained
inside this single-energy mode.

Metadata:

```text
spectral_model = monochromatic
```

### Power Law

```text
SPECTRAL_MODEL=powerlaw
```

Definition:

```text
dN/dE proportional to E^(-gamma)
```

with:

```text
gamma = SPECTRAL_GAMMA
```

### Power Law With Exponential Cutoff

```text
SPECTRAL_MODEL=powerlaw_cutoff
```

Definition:

```text
dN/dE proportional to E^(-gamma) exp(-E/Ecut)
```

with:

```text
gamma = SPECTRAL_GAMMA
Ecut  = SPECTRAL_ECUT_GEV
```

## Spectral Integration

For broad spectra, the code uses logarithmic energy bins between:

```text
SPECTRAL_E_MIN_GEV
SPECTRAL_E_MAX_GEV
```

with:

```text
SPECTRAL_N_BINS
```

For each energy bin, the existing single-energy Kerr/DIS transport calculation
is evaluated. The code then reports:

```text
<tau>    = integral Phi(E) tau(E) dE / integral Phi(E) dE
<Psurv>  = integral Phi(E) Psurv(E) dE / integral Phi(E) dE
I        = integral I(E) dE
```

where `Phi(E)` is the selected spectral weight.

## Makefile Parameters

```text
SPECTRAL_MODEL
SPECTRAL_GAMMA
SPECTRAL_ECUT_GEV
SPECTRAL_E_MIN_GEV
SPECTRAL_E_MAX_GEV
SPECTRAL_N_BINS
```

Safe defaults:

```text
SPECTRAL_MODEL=monochromatic
SPECTRAL_GAMMA=$(SOURCE_POWERLAW)
SPECTRAL_ECUT_GEV=$(SOURCE_EMAX_GEV)
SPECTRAL_E_MIN_GEV=1.0e5
SPECTRAL_E_MAX_GEV=1.0e12
SPECTRAL_N_BINS=8
```

## Example Commands

Backward-compatible single-energy image:

```bash
make image-from-small-cache SPECTRAL_MODEL=monochromatic NTHREADS=2
```

Power-law spectral image:

```bash
make image-from-small-cache \
  SPECTRAL_MODEL=powerlaw \
  SPECTRAL_GAMMA=2.0 \
  SPECTRAL_E_MIN_GEV=1e5 \
  SPECTRAL_E_MAX_GEV=1e12 \
  SPECTRAL_N_BINS=8 \
  NTHREADS=2
```

Power law with cutoff:

```bash
make image-from-small-cache \
  SPECTRAL_MODEL=powerlaw_cutoff \
  SPECTRAL_GAMMA=2.0 \
  SPECTRAL_ECUT_GEV=1e12 \
  SPECTRAL_N_BINS=8 \
  NTHREADS=2
```

Validation:

```bash
make validate_spectra NTHREADS=2 PYTHON=python3
```

Outputs:

```text
output/spectra/spectral_validation.csv
output/spectra/spectral_validation.md
plots/spectra/spectral_models.png
plots/spectra/survival_probability_vs_energy.png
plots/spectra/spectral_weighted_vs_monochromatic.png
```

## Energy-Band Composite Image

The multiband image product is an energy-band composite image. It combines
several monochromatic UHE image slices into three false-color energy intervals:

```text
blue  = low-energy UHE band
green = intermediate-energy UHE band
red   = high-energy UHE band
```

For the default command, the exact bands are:

```text
blue  = 1e5-1e7 GeV
green = 1e7-1e10 GeV
red   = 1e10-1e12 GeV
```

These are false colors encoding neutrino energy intervals. They are not
physical photon colors.

Command:

```bash
make plot_multiband_image \
  SPECTRAL_MODEL=powerlaw_cutoff \
  SPECTRAL_E_MIN_GEV=1e5 \
  SPECTRAL_E_MAX_GEV=1e12 \
  SPECTRAL_N_BINS=8 \
  NTHREADS=2
```

Outputs:

```text
plots/spectra/observed_multiband_image.png
output/spectra/observed_multiband_flux.csv
```

The CSV includes metadata lines specifying the product name, the false-color
interpretation, and the exact energy bands.

## Adding a New Spectrum

To add a new spectrum:

1. Add a new value to `UHESpectralModel`.
2. Extend `parse_uhe_spectral_model`.
3. Implement the new branch in `uhe_spectral_weight`.

No changes should be needed in:

```text
DIS cross sections
sigma tables
optical-depth integration
Kerr geodesics
density backgrounds
source morphologies
```

Future models can include broken power laws, log-parabolas,
Fermi-Dirac-like spectra, tabulated spectra, or user-defined spectra.
