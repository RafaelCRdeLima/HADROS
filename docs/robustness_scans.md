# Robustness and sensitivity scans

Point 3 tests whether the qualitative opacity conclusions survive controlled
parameter variations. The scans are validation-scale studies, not production
convergence tests.

## Scientific motivation

The paper argues that UHE neutrino opacity depends on:

- baryonic column density;
- geometry;
- UHE source morphology;
- DIS model choice.

The scan suite checks whether trends in optical depth, survival probability,
and image intensity are stable under reasonable changes in source and density
parameters.

## Command

```bash
make robustness_scans NTHREADS=2 \
  PYTHON=/home/rafael/micromamba/envs/dis/bin/python
```

Outputs:

```text
output/scans/scan_summary.csv
output/scans/scan_summary.md
plots/scans/
```

## Measured observables

Every scan case records:

- `mean_tau`;
- `max_tau`;
- `mean_P_surv`;
- `total_intensity`;
- `valid_rays`.

## Scan definitions

Scan A, source model dependence:

Fixed Gaussian background, spin, observer, and energy. Compare:

- `inner_ring`;
- `jet_base`;
- `funnel_wall`;
- `density_weighted`;
- `shock_layer`.

Question: do conclusions survive changing source morphology?

Scan B, density background dependence:

Fixed `funnel_wall` source. Compare:

- `gaussian`;
- `powerlaw`;
- `powerlaw_funnel`;
- `powerlaw_funnel_envelope`.

Question: are conclusions robust against background morphology?

Scan C, energy dependence:

For each source model, scan UHE energies and record `tau(E)` and `P_surv(E)`.

Question: does opacity increase with energy over the sampled UHE range?

Scan D, observer inclination:

Scan near-polar to near-equatorial observer inclinations.

Question: how strongly does viewing angle affect optical depth and image
intensity?

Scan E, black-hole spin:

Run `a = 0.0, 0.5, 0.9, 0.99` with fixed source and background.

Question: does spin significantly modify opacity maps or mainly morphology and
ray sampling?

Scan F, parameter stability:

For `density_weighted`, vary `q` and `s`.

For `funnel_wall`, vary wall angle and angular width.

Question: do observables change smoothly without numerical jumps?

## Interpretation guidelines

Safe statements are supported when trends persist across scan families at
validation resolution. Avoid claiming production convergence, astrophysical
uniqueness, or self-consistent source physics from these scans alone.

The `shock_layer` source is the implemented name for the density-gradient
source:

```text
shock_layer == density_gradient source
```
