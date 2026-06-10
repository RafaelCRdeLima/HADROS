# UHE source models

Point 2 implements a selectable family of phenomenological UHE source
prescriptions. These models are controlled emissivity morphologies used to test
whether opacity maps and survival probabilities are robust against source
location assumptions.

They must not be described as realistic particle acceleration, first-principles
cosmic-ray transport, or self-consistent jet physics.

## Common notation

Let `j_UHE(r, theta, phi, E)` be the comoving UHE emissivity proxy. The current
radiative-transfer code uses an energy dependence of the form:

```text
j_E proportional to E^(-source_powerlaw) * exp(-E/source_emax)
```

The implementation keeps that spectral factor and varies only the spatial
factor:

```text
j_UHE = source_norm * S(r, theta, phi) *
        E^(-source_powerlaw) * exp(-E/source_emax)
```

where `S` is one of the source morphologies below.

## A. Inner equatorial ring

Status: implemented as `SOURCE_MODEL=inner_ring`. This is the default and
preserves the original source behavior.

Mathematical definition:

```text
S_ring(r, theta) =
  exp[-((r - source_r)/source_sigma_r)^2] *
  exp[-((theta - pi/2)/source_theta_width)^2]
```

Physical interpretation:

This is a compact equatorial injection region inside the torus. It is a
controlled proxy for UHE production near the inner engine, not a resolved
emission model.

Expected impact on opacity maps:

The rays that cross the equatorial torus before reaching the observer should
carry the strongest attenuation. Inclined observers should see a strong
dependence on line-of-sight passage through dense equatorial material.

Tunable parameters:

- `source_r`
- `source_sigma_r`
- `source_theta_width`
- `source_norm`
- `source_powerlaw`
- `source_emax`

Implementation status:

This model is the baseline branch in CPU and CUDA.

## B. Compact jet-base source

Status: implemented as `SOURCE_MODEL=jet_base`.

Mathematical definition:

For a bipolar source near the spin axis:

```text
S_jet_base(r, theta) =
  exp[-((r - source_r)/source_sigma_r)^2] *
  { exp[-(theta/source_theta_width)^2]
    + exp[-((pi - theta)/source_theta_width)^2] }
```

Physical interpretation:

This places UHE production near the base of a polar outflow. It is a
parameterized jet-base proxy, not a jet simulation.

Expected impact on opacity maps:

Polar lines of sight may become brighter or less attenuated, especially in
models with evacuated funnels. The contrast between equatorial and polar
viewing angles should increase.

Tunable parameters:

- `source_r`
- `source_sigma_r`
- `source_theta_width`
- optional north/south asymmetry
- `source_norm`

Implementation status:

Implemented in CPU and CUDA as a bipolar compact source around the rotation
axis.

## C. Funnel-wall source

Status: implemented as `SOURCE_MODEL=funnel_wall`.

Mathematical definition:

For emission concentrated around a funnel half-opening angle `theta_wall`:

```text
S_wall(r, theta) =
  exp[-((r - source_r)/source_sigma_r)^2] *
  { exp[-((theta - theta_wall)/source_theta_width)^2]
    + exp[-((theta - (pi - theta_wall))/source_theta_width)^2] }
```

Physical interpretation:

This represents UHE production along a shear/interface region between a
low-density polar funnel and the denser torus/envelope. It is a controlled
interface proxy, not a shock calculation.

Expected impact on opacity maps:

The source lies close to density gradients, so optical-depth maps should show
strong angular structure. Some rays begin near low-density regions while others
immediately cross the torus wall.

Tunable parameters:

- `theta_wall`
- `source_theta_width`
- `source_r`
- `source_sigma_r`
- optional radial taper or maximum radius

Implementation status:

Implemented in CPU and CUDA. Metadata records `source_funnel_theta_deg`.

## D. Density-weighted source

Status: implemented as `SOURCE_MODEL=density_weighted`.

Mathematical definition:

```text
S_density(r, theta) =
  [rho(r, theta)/rho0]^q *
  (r/r0)^(-s)
```

with optional radial limits:

```text
S_density = 0 outside source_r_min <= r <= source_r_max
```

Physical interpretation:

This tests UHE production tied to local matter density and compactness. It is a
phenomenological emissivity proxy.

Expected impact on opacity maps:

For `q > 0`, emission is weighted toward denser regions, so the observed image
may become more attenuated and more sensitive to DIS opacity. Larger `s`
concentrates the source inward.

Tunable parameters:

- `source_density_power = q`
- `source_radial_power = s`
- `source_r_min`
- `source_r_max`
- `source_norm`

Implementation status:

Implemented in CPU and CUDA using the active density background and the
parameters `SOURCE_DENSITY_Q` and `SOURCE_RADIAL_S`.

## E. Density-gradient source

Status: implemented as `SOURCE_MODEL=shock_layer`.

Naming convention: `shock_layer` is the implemented model name for the
density-gradient source. In text, use "shock-layer/density-gradient source" on
first mention to make the equivalence unambiguous.

Mathematical definition:

```text
S_gradient(r, theta) proportional to |grad rho|
```

A practical dimensionless form is:

```text
S_gradient =
  sqrt[(d ln rho / d ln r)^2 +
       (d ln rho / d theta)^2] *
  radial_taper(r)
```

or, for direct density contrast:

```text
S_gradient = |grad rho| / rho0
```

Physical interpretation:

This targets interfaces where the density changes rapidly, such as funnel
walls or torus/envelope transitions. It is a controlled interface-emission
proxy.

Expected impact on opacity maps:

Emission should concentrate along sharp structural boundaries. The resulting
maps may highlight silhouette-like edges and strong viewing-angle dependence.

Tunable parameters:

- finite-difference step in `r`
- finite-difference step in `theta`
- optional smoothing scale
- radial limits
- gradient normalization

Implementation status:

Implemented in CPU and CUDA using finite differences on the active density
background. The calculation uses the floored density consistently.

## Makefile parameters

- `SOURCE_MODEL`: `inner_ring`, `funnel_wall`, `jet_base`, `shock_layer`, or
  `density_weighted`.
- `SOURCE_R_RG`: characteristic source radius.
- `SOURCE_SIGMA_RG`: radial width for localized source models.
- `SOURCE_THETA_DEG`: angular width.
- `SOURCE_POWERLAW`: UHE spectral power-law index.
- `SOURCE_EMAX_GEV`: UHE exponential cutoff.
- `SOURCE_NORM`: source normalization.
- `SOURCE_FUNNEL_THETA_DEG`: funnel-wall half-opening angle.
- `SOURCE_DENSITY_Q`: density exponent for `density_weighted`.
- `SOURCE_RADIAL_S`: radial exponent for `density_weighted`.
- `SOURCE_GRADIENT_DR_RG`: radial finite-difference step for `shock_layer`.
- `SOURCE_GRADIENT_DTHETA_DEG`: angular finite-difference step for
  `shock_layer`.
- `SOURCE_RHO_REF`: density reference for `density_weighted`; values below zero
  use `rho(SOURCE_R_RG, pi/2)` automatically. The safe Makefile default uses
  `TORUS_RHO0` to avoid accidentally normalizing by the density floor.
- `SOURCE_CUTOFF_MIN`, `SOURCE_CUTOFF_MAX`: numerical cutoffs applied to the
  `density_weighted` and `shock_layer` spatial factors. The default upper
  cutoff is `1e2`.

## Metadata

Image outputs record:

```text
source_model
source_r_rg
source_sigma_rg
source_theta_deg
source_powerlaw
source_emax_GeV
source_norm
source_funnel_theta_deg
source_rho_ref
rho_ref
source_q
source_s
source_density_power_q
source_radial_power_s
source_cutoff_min
source_cutoff_max
source_gradient_dr_rg
source_gradient_dtheta_deg
```

## Validation status

Implemented checks:

- CPU build passes.
- CUDA build passes.
- All five source models run on the CPU using the medium validation cache.
- CUDA route was smoke-tested for `funnel_wall` using the medium CUDA cache.
- Morphology diagnostic plots are generated by:

```bash
make validate_source_plots
```

Output directory:

```text
plots/validation_uhe_sources/
```

Recommended additional validation:

Each source model should include:

- backward compatibility for `SOURCE_MODEL=inner_ring`;
- source-emissivity diagnostic maps;
- metadata recording all source parameters;
- CPU/CUDA consistency checks if CUDA paths compute emissivity independently;
- at least one comparison across `gaussian` and `powerlaw_funnel_envelope`
  density backgrounds.
