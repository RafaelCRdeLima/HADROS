# Phenomenological UHE source prescriptions

The code implements selectable phenomenological UHE source prescriptions through
`SOURCE_MODEL`. They are controlled emissivity morphologies for opacity studies,
not realistic particle acceleration, first-principles cosmic-ray transport, or
self-consistent jet physics.

Available models:

- `inner_ring`: original compact equatorial source.
- `funnel_wall`: emission near the polar funnel boundary.
- `jet_base`: compact bipolar source near the rotation axis and small radii.
- `density_weighted`: distributed source with `j_UHE` proportional to
  `rho^q r^{-s}`.
- `shock_layer`: density-gradient source. This is equivalent to the
  "density_gradient" source discussed in the paper notes.

For clarity, write:

```text
shock_layer == density_gradient source
```

The density-weighted source records and uses:

- `rho_ref`
- `source_q`
- `source_s`
- `source_cutoff_min`
- `source_cutoff_max`

The cutoffs prevent extreme emissivity spikes during parameter scans for
`density_weighted` and `shock_layer`.
The default scan configuration uses `rho_ref = TORUS_RHO0` and
`source_cutoff_max = 1e2`.

Diagnostics:

```bash
make validate_source_plots
```

Detailed mathematical definitions are in:

```text
docs/uhe_source_models_design.md
```
