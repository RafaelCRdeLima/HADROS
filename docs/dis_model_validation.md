# DIS Model Validation

This validation freezes the astrophysical model and varies only the
neutrino-nucleon DIS cross-section table used by the UHE opacity calculation.

The goal is to quantify how strongly UHE observables depend on the assumed DIS
model:

- `GBW`
- `IIM`
- `PDF_reference`

No Kerr geodesics, camera model, density background, UHE source morphology,
UHE spectral model, MeV module, dashboard architecture, or ParaView export is
modified by this workflow.

## Sigma Tables

The local sigma tables are read by `src/sigma_table.cpp`.

Expected columns:

```text
Enu_GeV sigma_GeV_minus2 sigma_cm2
```

The loader uses log-log interpolation and rejects requests outside the tabulated
energy range. It does not extrapolate.

The current local files are scalar charged-current neutrino-nucleon tables:

```text
data/sigma/sigma_nuN_CC_GBW.dat
data/sigma/sigma_nuN_CC_IIM.dat
data/sigma/sigma_nuN_CC_PDF_reference.dat
```

The tables do not contain separate CC/NC, neutrino/antineutrino, or
structure-function columns. The `use_F3` flag is recorded in image metadata, but
xF3 inclusion cannot be inferred from the table alone.

## PDF Reference Curve

`PDF_reference` is a documented literature-scale charged-current reference
curve, tabulated locally on the same energy grid as GBW/IIM:

```text
sigma_CC(E) = 5.53e-36 * (E_GeV)^0.363 cm^2
```

This power-law reference is intended as a controlled PDF-like benchmark for the
validation suite. It is not a full CTW/CSMS PDF uncertainty band and should not
be described as such.

Literature context:

- Gandhi, Quigg, Reno & Sarcevic, Phys. Rev. D 58, 093009; UHE neutrino
  interactions over roughly 1e9-1e21 eV.
- Connolly, Thorne & Waters, Phys. Rev. D 83, 113009; MSTW PDF-based CC/NC
  neutrino cross sections and parametrizations for 1e4-1e12 GeV.
- Cooper-Sarkar, Mertsch & Sarkar, JHEP 08, 042; Standard Model high-energy
  neutrino cross sections and PDF uncertainties.

## Workflow

Run:

```bash
make validate_dis_models NTHREADS=2 PYTHON=python3
```

The workflow:

1. audits available DIS tables;
2. writes `data/sigma/sigma_nuN_CC_PDF_reference.dat`;
3. compares cross sections and ratios;
4. generates one 32x32 Kerr geodesic cache for validation;
5. runs the same astrophysical setup for GBW, IIM, and PDF_reference;
6. compares `tau`, `P_surv`, image intensity, and observed-spectrum proxies;
7. writes a human-readable interpretation.

## Outputs

Data products:

```text
output/dis_validation/dis_table_audit.md
output/dis_validation/sigma_model_comparison.csv
output/dis_validation/dis_observable_comparison.csv
output/dis_validation/dis_model_validation_summary.md
```

Figures:

```text
plots/dis_validation/sigma_nuN_model_comparison.png
plots/dis_validation/sigma_nuN_model_ratios.png
plots/dis_validation/observed_spectrum_dis_models.png
plots/dis_validation/attenuation_ratio_dis_models.png
plots/dis_validation/uhe_image_dis_model_comparison.png
plots/dis_validation/uhe_image_dis_model_ratio.png
```

## Interpretation

If two DIS tables differ by a factor `f` in sigma at fixed energy and the medium
is optically thin, dimensionless optical depth changes approximately by the
same factor while image intensity changes weakly. In optically thick regimes,
small sigma differences can produce large survival-probability differences.

Safe claims:

- The suite isolates DIS-table dependence because the ray cache, density
  background, source prescription, observer, spin, and energy are fixed inside
  each regime.
- The results can support statements about sensitivity to UHE DIS cross-section
  assumptions.

Claims to avoid:

- Do not claim that `PDF_reference` is a full modern PDF uncertainty band.
- Do not claim CC+NC coverage unless corresponding tables are added.
- Do not claim xF3 status from the sigma table alone.
