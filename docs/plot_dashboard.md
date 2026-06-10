# Plot Dashboard

The plot dashboard is a static HTML index for HADROS visual products. It
is an output-management tool only. It does not run simulations and does not
modify the DIS cross sections, optical-depth physics, Kerr geodesics, density
models, UHE source models, UHE spectral models, or MeV physics.

## Command

Generate the dashboard with:

```bash
make dashboard PYTHON=python3
```

Outputs:

```text
dashboard/index.html
dashboard/plot_manifest.json
```

Open `dashboard/index.html` directly in a browser. No server or internet access
is required.

## Purpose

The project now produces plots from several modules. The dashboard makes it
easier to answer:

- which plots currently exist,
- where they are stored,
- which module they belong to,
- which command probably generated them,
- which data products are associated with them,
- which ParaView files are available.

## Categories

The dashboard groups products into:

- Kerr images,
- density backgrounds,
- UHE source morphologies,
- robustness scans,
- opacity surfaces,
- spectral models,
- MeV neutrino physics,
- ParaView exports,
- other unclassified products.

## Discovery Rules

The generator scans:

```text
plots/
output/
```

Image-like products:

```text
.png
.jpg
.jpeg
.svg
.pdf
```

Data-like products:

```text
.csv
.dat
.txt
.md
.vtk
.vti
.vtu
```

Plots are categorized from their relative path and filename. For example,
`plots/scans/source_model_comparison.png` is assigned to robustness scans, while
`plots/mev_physics/mev_opacity_vs_energy.png` is assigned to MeV neutrino
physics.

The command field is filled only for known filename/path patterns. If no known
pattern exists, the dashboard writes:

```text
unknown
```

## Manifest

The machine-readable manifest is:

```text
dashboard/plot_manifest.json
```

Each entry contains:

```json
{
  "title": "...",
  "category": "...",
  "plot_path": "...",
  "data_path": "...",
  "description": "...",
  "related_module": "...",
  "command": "...",
  "created_at": "..."
}
```

ParaView files are included as data products even when they do not have an image
preview.

## Adding New Plot Patterns

To teach the dashboard about a new product family, edit:

```text
scripts/build_plot_dashboard.py
```

The relevant functions are:

```text
categorize()
describe_kerr()
describe_density()
describe_source()
describe_scan()
describe_opacity()
describe_spectral()
describe_mev()
describe_paraview()
```

Add a filename/path pattern to `categorize()` and a short caption in the
matching `describe_*()` helper. Keep descriptions concise and avoid inventing a
command unless the generating target is known.

## Validation

The generator validates that all indexed relative paths exist before writing the
dashboard. Run:

```bash
make dashboard PYTHON=python3
```

Then check:

```text
dashboard/index.html
dashboard/plot_manifest.json
```

The HTML uses relative links and no external JavaScript or CSS dependencies.

## Limitations

The dashboard infers categories from file paths and names, so unfamiliar
filenames may be placed under `Other`.

Associated data files are matched heuristically by directory, stem, and module.
If the exact data provenance is not encoded in the filename, the dashboard may
leave `data_path` empty or list nearby summary products.

The dashboard is not a provenance database. For publication figures, still
record the exact Makefile command and parameter set in the paper notes or output
metadata.
