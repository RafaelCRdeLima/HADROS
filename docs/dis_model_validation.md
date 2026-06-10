# DIS Model Validation

This validation freezes the astrophysical model and varies only the
neutrino-nucleon DIS cross-section table used by the UHE opacity calculation.

Compared models:

- `GBW`
- `IIM`
- `literature_powerlaw_scale`
- `CTW_reference`

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

Current local files:

```text
data/sigma/sigma_nuN_CC_GBW.dat
data/sigma/sigma_nuN_CC_IIM.dat
data/sigma/sigma_nuN_CC_literature_powerlaw_scale.dat
data/sigma/sigma_nuN_CC_CTW_reference.dat
```

The propagation tables are scalar charged-current neutrino-nucleon tables. They
do not contain separate NC, antineutrino, or structure-function columns. The
`use_F3` flag is recorded in image metadata, but xF3 inclusion cannot be
inferred from the table alone.

## CTW Reference

`CTW_reference` uses the published central `nu N` charged-current table from:

```text
Connolly, Thorne & Waters,
Phys. Rev. D 83, 113009,
arXiv:1102.0691,
Table I.
```

Status:

- Source: published table values, not plot digitization.
- PDF set: MSTW 2008.
- Channel: charged current only.
- Particle: neutrino, not antineutrino.
- Target: isoscalar nucleon.
- Energy range: `1e4 <= E_nu/GeV <= 1e12`.
- Units: `cm^2`, with the corresponding `GeV^-2` column computed using
  `1 GeV^-2 = 0.389379338e-27 cm^2`.

This is a defensible literature/PDF-based reference for the current paper-level
DIS comparison. It is still not a full CTW uncertainty-band implementation.

## Literature Power-Law Scale

The old approximate curve has been renamed:

```text
literature_powerlaw_scale
```

with:

```text
sigma_CC(E) = 5.53e-36 * E_GeV^0.363 cm^2
```

It is retained only as a historical scale check. It must not be called
`PDF_reference` in paper text.

## Workflow

Run:

```bash
make validate_dis_models NTHREADS=2 PYTHON=python3
```

The workflow:

1. audits available DIS tables;
2. writes the CTW and power-law reference tables;
3. compares cross sections and ratios;
4. generates one 32x32 Kerr geodesic cache for validation;
5. runs the same astrophysical setup for all DIS models;
6. compares `tau`, `P_surv`, image intensity, and observed-spectrum proxies;
7. writes a human-readable interpretation.

## Outputs

Data products:

```text
output/dis_validation/dis_table_audit.md
output/dis_validation/reference_model_audit.md
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
- The CTW comparison supports statements about sensitivity to the assumed UHE
  charged-current neutrino-nucleon cross section relative to a published
  MSTW-2008 Standard Model DIS calculation.

Claims to avoid:

- Do not call `literature_powerlaw_scale` a PDF reference.
- Do not claim CTW PDF uncertainty-band coverage yet.
- Do not claim CC+NC or neutrino+antineutrino coverage unless corresponding
  tables are added.
- Do not claim xF3 status from the sigma table alone.
