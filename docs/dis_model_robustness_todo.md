# DIS-model robustness TODO

DIS-model robustness is not part of the completed Point-3 robustness scan suite.

Point 3 currently tests sensitivity to:

- UHE source morphology;
- density-background morphology;
- neutrino energy;
- observer inclination;
- black-hole spin;
- selected source parameters.

It does not yet establish robustness against the high-energy DIS prescription.

A dedicated DIS-model comparison should be run before claiming DIS robustness in
the paper. At minimum, this should compare:

- GBW;
- IIM;
- one or more modern PDF-based neutrino-nucleon cross-section tables, if
  available.

The comparison should use the same geodesic cache, density background, source
model, observer, spin, and energy grid, varying only the DIS input table/model.

Required observables:

- `mean_tau`;
- `max_tau`;
- `mean_P_surv`;
- `total_intensity`;
- energy dependence of tau and survival probability.

Until this is done, safe wording is:

```text
The current robustness scans test geometry, density, source morphology,
inclination, spin, and energy dependence. Robustness to the DIS model remains a
separate validation task.
```
