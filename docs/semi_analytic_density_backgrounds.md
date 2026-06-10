# Semi-analytic density backgrounds

This module implements controlled semi-analytic density morphologies for UHE
neutrino opacity studies. The profiles are parameterized density fields used
for post-processing and sensitivity tests.

These backgrounds must not be described as hydrodynamical simulations,
realistic collapsar solutions, or equilibrium accretion flows. They are
controlled environments designed to test how optical depth and survival
probability respond to density morphology.

## Coordinates and notation

All radii are in units of `r_g = GM/c^2`. The polar angle is `theta`, with the
equatorial plane at `theta = pi/2`.

The implementation computes:

```text
rho(r, theta) = max(rho_raw(r, theta), rho_floor)
```

where `rho_floor >= 0` prevents exactly zero-density regions in funnels.

## Common factors

The Gaussian disk/torus shape is:

```text
G(r, theta) =
  exp[-((r - r0)/sigma_r)^2] *
  exp[-((theta - pi/2)/H_over_R)^2]
```

with the legacy hard cuts:

```text
4 < r < 18
abs(theta - pi/2) < 0.45
```

The power-law disk/torus shape is:

```text
P(r, theta) =
  (r/r0)^(-radial_power) *
  exp[-(cos(theta)/H_over_R)^2] *
  [1 - exp(-((r - R_min)/sigma_r)^2)] *
  exp[-(r/R_max)^4]
```

inside `R_min <= r <= R_max`; outside that interval it is zero.

The bipolar funnel depletion factor is:

```text
F(theta) =
  clamp(1 - funnel_depletion *
    [exp(-(theta/funnel_theta)^2) +
     exp(-((pi - theta)/funnel_theta)^2)],
    0, 1)
```

The radial envelope is:

```text
E(r) = envelope_rho0 * (r/r0)^(-envelope_alpha)
```

for `r0 <= r <= R_max`; otherwise it is zero.

## Implemented profiles

`gaussian`:

```text
rho_raw = rho0 * G(r, theta)
```

This is the legacy thick-torus control model.

`powerlaw`:

```text
rho_raw = rho0 * P(r, theta)
```

This gives an equatorially concentrated, radially declining controlled
background.

`gaussian_funnel`:

```text
rho_raw = rho0 * G(r, theta) * F(theta)
```

This keeps the legacy Gaussian torus while adding polar depletion.

`powerlaw_funnel`:

```text
rho_raw = rho0 * P(r, theta) * F(theta)
```

This combines the power-law disk/torus with a depleted polar funnel.

`gaussian_envelope`:

```text
rho_raw = rho0 * G(r, theta) + E(r)
```

This adds a controlled external envelope to the Gaussian torus.

`powerlaw_envelope`:

```text
rho_raw = rho0 * P(r, theta) + E(r)
```

This adds a controlled external envelope to the power-law disk/torus.

`powerlaw_funnel_envelope`:

```text
rho_raw = rho0 * P(r, theta) * F(theta) + E(r)
```

This is the most structured Point-1 background: radial power law, polar funnel,
and external envelope.

`collapsar_ndaf_like`:

```text
rho_raw = rho0 * C(r, theta)
```

where `C(r, theta)` is a smooth, equatorially concentrated power-law-like
morphology with inner taper, outer taper, and an inner dense enhancement. It is
intended as a controlled semi-analytic collapsar/NDAF-like preset with larger
emitting mass and higher density over a broader volume than the fiducial
Gaussian torus.

This is not a hydrodynamical collapsar simulation, an equilibrium accretion
solution, or a calibrated NDAF model. It is a literature-guided
semi-analytic density morphology for post-processing tests.

## Makefile parameters

Density parameters:

- `DENSITY_PROFILE`: profile name.
- `TORUS_RHO0`: `rho0`, density normalization in g/cm^3.
- `TORUS_R0_RG`: `r0`.
- `TORUS_SIGMA_RG`: `sigma_r`.
- `TORUS_H_OVER_R`: vertical angular thickness parameter.
- `TORUS_RADIAL_POWER`: radial power-law index.
- `FUNNEL_DEPLETION`: polar depletion amplitude between 0 and 1.
- `FUNNEL_THETA_DEG`: funnel angular width in degrees.
- `ENVELOPE_RHO0`: envelope normalization in g/cm^3.
- `ENVELOPE_ALPHA`: envelope radial power-law index.
- `TORUS_R_MIN_RG`: minimum radial domain for power-law/envelope profiles.
- `TORUS_R_MAX_RG`: maximum radial domain for power-law/envelope profiles.
- `RHO_FLOOR`: explicit density floor in g/cm^3.

Camera and black-hole parameters:

- `ASPIN`: Kerr spin parameter.
- `MBH_MSUN`: black-hole mass in solar masses.
- `CAM_R_OBS_RG`: observer distance.
- `CAM_THETA_DEG`: observer inclination in degrees.
- `CAM_FOV_DEG`: camera field of view in degrees.
- `CAM_NX`, `CAM_NY`: production camera grid.
- `CAM_R_MAX_RG`: maximum integration radius.
- `CAM_STEP`: ray-integration step parameter.

Current UHE source parameters:

- `SOURCE_MODEL`: phenomenological UHE source prescription.
- `SOURCE_R_RG`: current inner ring radius.
- `SOURCE_SIGMA_RG`: current inner ring radial width.
- `SOURCE_THETA_DEG`: current inner ring angular width.
- `SOURCE_POWERLAW`: spectral power-law index.
- `SOURCE_EMAX_GEV`: source cutoff energy in GeV.
- `SOURCE_NORM`: source normalization.
- `SOURCE_FUNNEL_THETA_DEG`: funnel-wall half-opening angle.
- `SOURCE_DENSITY_Q`: density exponent for `density_weighted`.
- `SOURCE_RADIAL_S`: radial exponent for `density_weighted`.
- `SOURCE_RHO_REF`: density reference for `density_weighted`; negative values
  use the automatic equatorial reference.
- `SOURCE_CUTOFF_MIN`, `SOURCE_CUTOFF_MAX`: numerical safety cutoffs for the
  `density_weighted` spatial factor.
- `SOURCE_GRADIENT_DR_RG`: radial finite-difference step for `shock_layer`.
- `SOURCE_GRADIENT_DTHETA_DEG`: angular finite-difference step for `shock_layer`.

Validation and runtime parameters:

- `NTHREADS`: requested OpenMP thread count.
- `OMP_NUM_THREADS`: exported OpenMP thread count.
- `SMALL_NX`, `SMALL_NY`: small validation grid.
- `SMALL_R_MAX_RG`, `SMALL_CAM_STEP`: small validation ray settings.
- `SMALL_TORUS_NR`, `SMALL_TORUS_NTH`: density diagnostic grid.
- `SMALL_CACHE_PATH`: small CPU cache path.
- `MEDIUM_NX`, `MEDIUM_NY`: medium CUDA validation grid.
- `MEDIUM_R_MAX_RG`, `MEDIUM_CAM_STEP`: medium validation ray settings.
- `MEDIUM_CACHE_PATH`: medium CPU cache path used for CPU/CUDA comparison.
- `GPU_CACHE_RAYS`: CUDA cache rays per batch.
- `GPU_CACHE_MAX_POINTS`: maximum cached points per ray.
- `GPU_TOL`: CUDA geodesic RKF45 tolerance.
- `NVCC`: CUDA compiler path.
- `PYTHON`: Python interpreter used by plotting/validation scripts.
- `USE_F3`: records whether the DIS table uses the `xF3` contribution.

## Output metadata

Image outputs record the fields required for reproducibility:

```text
profile_type, rho0, r0, sigma_r, H_over_R, radial_power,
funnel_depletion, funnel_theta_deg, rho_floor, envelope_rho0,
envelope_alpha, R_min, R_max, spin, observer_distance,
observer_inclination, DIS_model, use_F3
```

The code also writes unit-suffixed aliases such as `rho0_gcm3`, `r0_rg`,
`sigma_r_rg`, `rho_floor_gcm3`, `R_min_rg`, and `R_max_rg`.

## Validation summary

Gaussian backward compatibility:

```text
mean_optical_depth          relative difference = 0
maximum_optical_depth       relative difference = 2.09e-16
mean_survival_probability   relative difference = 0
total_image_intensity       relative difference = 0
rays_reaching_torus         relative difference = 0
```

The Gaussian mode reproduces the legacy Gaussian formula on the controlled
32x32 cache within numerical precision.

Density floor:

```text
rho = max(rho_raw, rho_floor)
```

was verified for `gaussian_funnel`, `powerlaw_funnel`, and
`powerlaw_funnel_envelope` by forcing `rho_floor = 1e-12` and confirming that
the minimum sampled density is exactly `1e-12`.

CPU validation:

- Build target passes.
- Gaussian backward compatibility passes.
- Diagnostic plots are generated in `plots/validation_density_backgrounds/`.

CUDA validation:

- CUDA compilation works with `/home/rafael/micromamba/envs/dis/bin/nvcc`.
- 8x8 CUDA smoke test passes for opacity/survival metrics.
- 32x32 medium CUDA image validation passes for `gaussian` and
  `powerlaw_funnel_envelope` with nonzero total intensity.

Medium CUDA report:

```text
output/validation/cuda_cpu_comparison_medium.txt
```

Medium 32x32 results:

```text
gaussian:
  mean_tau relative difference        = -1.27e-16
  max_tau relative difference         = 0
  mean_P_surv relative difference     = 0
  total_intensity relative difference = 4.50e-6
  valid_rays CPU/CUDA                 = 904 / 904

powerlaw_funnel_envelope:
  mean_tau relative difference        = 0
  max_tau relative difference         = 0
  mean_P_surv relative difference     = 0
  total_intensity relative difference = 4.50e-6
  valid_rays CPU/CUDA                 = 904 / 904
```

The medium CUDA validation is not production resolution. Before using CUDA
figures in the paper, rerun the same comparison at the adopted science
resolution.

## Point-1 status

Point 1 is complete for the controlled semi-analytic density-background
framework and its CPU/CUDA validation cycle at validation resolution.
