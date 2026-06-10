# ParaView Export

## Motivation

This export is intended for visual inspection of the current BH_Torus_RTX
semi-analytic density and UHE source prescriptions in ParaView. It converts
the existing axisymmetric fields into a regular 3D Cartesian sampling so that
volume rendering, slices, and isosurfaces can be inspected interactively.

This is a 3D Cartesian sampling of an axisymmetric semi-analytic model, not a
fully 3D hydrodynamical simulation.

## Coordinate Transformation

The model fields are evaluated as axisymmetric functions of spherical radius
and polar angle:

```text
rho = rho(r, theta)
j_UHE = j_UHE(r, theta)
```

For each Cartesian grid point:

```text
r     = sqrt(x^2 + y^2 + z^2)
theta = acos(z / r)
phi   = atan2(y, x)
```

At `r = 0`, the exporter sets:

```text
theta = 0
phi   = 0
```

The current model ignores `phi` internally, but `phi` is exported for reference
and for future compatibility with truly 3D backgrounds.

## VTK Format

The first implementation writes a legacy VTK ASCII structured-points file:

```text
output/paraview/bh_torus_fields.vtk
```

The Cartesian domain is:

```text
x, y, z in [-PARAVIEW_BOX_RG, +PARAVIEW_BOX_RG]
```

with resolution controlled by:

```text
PARAVIEW_NX
PARAVIEW_NY
PARAVIEW_NZ
PARAVIEW_BOX_RG
```

## Exported Fields

The exporter writes:

```text
density_gcm3
log10_density
uhe_emissivity
log10_uhe_emissivity
r_rg
theta
phi
normalized_source
```

`normalized_source` is the dimensionless source morphology normalized by the
maximum value on the exported grid.

The logarithmic fields use numerical floors so that the output contains no
NaN or Inf values.

## Command

Example:

```bash
make -C BH_Torus_RTX paraview_fields \
  DENSITY_PROFILE=powerlaw_funnel_envelope \
  SOURCE_MODEL=funnel_wall \
  PARAVIEW_NX=64 PARAVIEW_NY=64 PARAVIEW_NZ=64 \
  PARAVIEW_BOX_RG=80 \
  NTHREADS=4
```

The target respects:

```text
NTHREADS
OMP_NUM_THREADS
```

## Opening in ParaView

Open:

```text
BH_Torus_RTX/output/paraview/bh_torus_fields.vtk
```

Suggested visualizations:

- Volume rendering of `log10_density`.
- Slices in the `xy`, `xz`, and `yz` planes.
- Isosurfaces of `density_gcm3`.
- Volume rendering of `uhe_emissivity`.
- Overlaid density isosurfaces and source-emissivity volume rendering.

## Limitations

This export does not introduce new physics and does not modify the physical
model. It only samples the existing axisymmetric density and UHE emissivity
prescriptions on a Cartesian grid.

It should not be described as a hydrodynamical simulation, a realistic
collapsar simulation, or a fully 3D accretion flow.

Exporting tau and P_surv fields requires a well-defined line-of-sight or
radial optical-depth convention and will be handled separately in the
opacity-surface module.

Future work can add truly `phi`-dependent 3D backgrounds and vector/tensor
fields once those quantities exist in the physical model.
