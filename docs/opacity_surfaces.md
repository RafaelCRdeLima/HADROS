# Opacity Surfaces

Point 4 characterizes the optical-depth structure of the medium. The extracted
surfaces are diagnostic UHE neutrino opacity surfaces, analogous to
photospheres, defined from DIS attenuation through the parameterized density
background.

Opacity surfaces are properties of the medium. For fixed density field,
geometry, energy, and DIS cross section, they do not depend on the UHE source
morphology.

## Mathematical Definition

For the current axisymmetric implementation:

```text
rho = rho(r, theta)
```

the radial outward DIS optical depth is:

```text
tau_E(r, theta) = int_r^Rmax [rho(r', theta) / m_u] sigma_DIS(E) dr'
```

where `dr'` is converted from `r_g` to cm using the selected black-hole mass.

The opacity surface is:

```text
r_tau(theta, E)
```

such that:

```text
tau_E(r_tau, theta) = tau_target
```

The standard targets are:

```text
tau_target = 2/3, 1, 3
```

The code is organized around `r_tau(theta)` files so that a later extension can
support:

```text
r_tau(theta, phi)
```

without changing the interpretation of the axisymmetric products.

## Extraction Method

For each angular sample `theta`, the code samples a uniform radial grid between
`R_min` and `R_max`. It computes cumulative outward optical depth from `R_max`
inward:

```text
tau(R_max, theta) = 0
```

and uses trapezoidal radial integration. The crossing radius is obtained by
linear interpolation in `tau` between the two neighboring radial samples that
bracket the target value.

If the total optical depth never reaches the requested target, the file records:

```text
r_tau_rg = -1
crossing_found = 0
```

No artificial crossing is created.

## Output Columns

Each surface file contains metadata plus:

```text
theta_rad theta_deg r_tau_rg tau_total classification crossing_found
```

The `classification` code is:

```text
0  transparent    tau < 0.1
1  transition     0.1 <= tau < 1
2  semi-opaque    1 <= tau < 3
3  opaque         tau >= 3
```

## Command

```bash
make opacity_surfaces TAU_SURFACE_VALUE=1.0 \
  PYTHON=/home/rafael/micromamba/envs/dis/bin/python
```

Useful controls:

```text
ENU
DENSITY_PROFILE
TORUS_RHO0
TORUS_R0_RG
TORUS_SIGMA_RG
TORUS_H_OVER_R
TORUS_RADIAL_POWER
FUNNEL_DEPLETION
FUNNEL_THETA_DEG
ENVELOPE_RHO0
ENVELOPE_ALPHA
TORUS_R_MIN_RG
TORUS_R_MAX_RG
RHO_FLOOR
TAU_SURFACE_VALUE
OPACITY_SURFACE_NTHETA
OPACITY_SURFACE_NR
```

## Generated Products

Surface files:

```text
output/opacity_surfaces/tau_surface_tau067_E1e5.dat
output/opacity_surfaces/tau_surface_tau1_E1e5.dat
output/opacity_surfaces/tau_surface_tau3_E1e5.dat
output/opacity_surfaces/tau_surface_tau067_E1e7.dat
output/opacity_surfaces/tau_surface_tau1_E1e7.dat
output/opacity_surfaces/tau_surface_tau3_E1e7.dat
output/opacity_surfaces/tau_surface_tau067_E1e9.dat
output/opacity_surfaces/tau_surface_tau1_E1e9.dat
output/opacity_surfaces/tau_surface_tau3_E1e9.dat
output/opacity_surfaces/tau_surface_tau067_E1e11.dat
output/opacity_surfaces/tau_surface_tau1_E1e11.dat
output/opacity_surfaces/tau_surface_tau3_E1e11.dat
output/opacity_surfaces/tau_surface_tau067_E1e12.dat
output/opacity_surfaces/tau_surface_tau1_E1e12.dat
output/opacity_surfaces/tau_surface_tau3_E1e12.dat
```

Compatibility copies are also written:

```text
output/opacity_surfaces/tau_surface_tau067.dat
output/opacity_surfaces/tau_surface_tau1.dat
output/opacity_surfaces/tau_surface_tau3.dat
```

Summary and validation:

```text
output/opacity_surfaces/energy_summary.csv
output/opacity_surfaces/background_summary.csv
output/opacity_surfaces/source_independence_check.txt
output/opacity_surfaces/validation_report.txt
output/opacity_surfaces/opacity_surface_summary.md
```

Figures:

```text
plots/opacity_surfaces/tau_surface_comparison.png
plots/opacity_surfaces/energy_dependence_tau_surface.png
plots/opacity_surfaces/opacity_phase_diagram.png
plots/opacity_surfaces/p_surv_energy_theta.png
plots/opacity_surfaces/source_independence_test.png
plots/opacity_surfaces/background_surface_comparison.png
plots/opacity_surfaces/mean_tau_energy.png
plots/opacity_surfaces/opacity_classification_energy_theta.png
```

## Source-Independence Test

The validation labels the same fixed-medium extraction with:

```text
inner_ring
jet_base
funnel_wall
density_weighted
shock_layer
```

and compares `r_tau(theta)` across all labels. The expected and required result
is:

```text
max_abs_delta_r_tau = 0
```

This proves that the opacity-surface product is tied to density, geometry,
energy, and the DIS table, not to the emissivity prescription.

## ParaView Workflow

The opacity module exports extracted surfaces as legacy VTK PolyData files:

```text
output/opacity_surfaces/paraview/tau_surface_tau1_E1e9.vtk
```

These files are surfaces of revolution generated from the extracted
axisymmetric `r_tau(theta)` curves. They are suitable for ParaView surface
inspection.

The classification product is also exported:

```text
output/opacity_surfaces/paraview/opacity_classification_regions.vtk
```

Do not interpret this as a volumetric tau field. Exporting tau and `P_surv` as
3D local scalar fields requires a well-defined line-of-sight or radial
optical-depth convention and is not faked here.

## Limitations

- Current extraction is axisymmetric: `r_tau(theta)`.
- The path is radial and outward, not a full null-geodesic camera integral.
- The result depends on the selected DIS table and density morphology.
- Spin and observer angle affect camera/geodesic optical depths, but this
  radial diagnostic surface does not yet include a geodesic-dependent surface
  extraction.
- A meaningful UHE neutrinosphere exists only where the sampled medium reaches
  the requested optical depth. Transparent configurations correctly return no
  crossing.
