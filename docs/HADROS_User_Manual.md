# HADROS User Manual

High-energy Astrophysical DIS Radiation and Observation Simulator

This manual is written for students and new users. It explains first the
physical ideas behind HADROS, then the main commands, and only after that the
internal program flow. HADROS is a research code for controlled opacity studies,
so the manual is careful about what the code can and cannot claim.

![HADROS workflow](figures/manual_pipeline_scheme.png)

## 1. Overview

HADROS combines Kerr ray tracing, semi-analytic astrophysical backgrounds, UHE
neutrino DIS opacity, phenomenological source and spectral prescriptions,
diagnostic MeV neutrino post-processing, ParaView export, and a static dashboard
for organizing plots.

At a high level the pipeline is:

```text
source prescription
  -> Kerr ray tracing
  -> density / temperature / Ye background
  -> opacity and radiative transfer
  -> images, plots, CSV tables, dashboard
```

HADROS does not perform hydrodynamic evolution. It is not full neutrino
radiation hydrodynamics, not a calibrated collapsar simulation, and not a
complete PDF-uncertainty package. The current backgrounds are semi-analytic
density morphologies, the UHE sources are phenomenological prescriptions, and
the MeV luminosities are diagnostic proxies.

## 2. Physical Theories

### 2.1 Kerr ray tracing

The spacetime around a rotating black hole is modeled with the Kerr metric. The
main physical parameters are the black-hole mass `MBH_MSUN` and the
dimensionless spin `ASPIN`. The observer is placed at a large radius
`CAM_R_OBS_RG` and inclination `CAM_THETA_DEG`. The camera launches rays
backward through the Kerr geometry and stores geodesics in a cache.

For neutrinos in the UHE range, the rest mass is negligible compared with the
energy, so the propagation is approximated by null geodesics:

```text
ds^2 = 0
```

The ray path is therefore controlled by the Kerr geometry. The energy
dependence in UHE images mostly enters through opacity,
`sigma_nuN(E)`, not through chromatic lensing. In other words, a low-energy and
high-energy neutrino see essentially the same null geodesic, but not the same
survival probability.

### 2.2 Semi-analytic density backgrounds

The matter field is a controlled semi-analytic background, not a hydrodynamic
solution. Available density profiles include:

- `gaussian`
- `powerlaw`
- `gaussian_funnel`
- `powerlaw_funnel`
- `gaussian_envelope`
- `powerlaw_envelope`
- `powerlaw_funnel_envelope`
- `collapsar_ndaf_like`

Every profile applies an explicit density floor:

```text
rho = max(rho_raw, rho_floor)
```

This prevents numerical zero-density funnels from becoming artificially
transparent.

| Regime | Purpose | Typical density | Interpretation | Limitation |
|---|---|---:|---|---|
| `fiducial_uhe_default` | UHE opacity baseline | low, often near `1e-2 g cm^-3` in validation examples | transparent controlled test background | not a collapsar disk |
| `fiducial_mev_density` | MeV diagnostic baseline | around `1e9-1e10 g cm^-3` in audits | weakly neutrino-cooled compact torus | low luminosity is expected |
| `collapsar_ndaf_like` | literature-guided high-density preset | significant regions near `1e10-1e12 g cm^-3` | semi-analytic NDAF/collapsar-like thermodynamic regime | not hydrodynamics and not tuned to a target luminosity |

### 2.3 UHE neutrino DIS opacity

The UHE optical depth is computed from the baryon column along a ray:

```text
tau_DIS(E) = integral n_b(s) sigma_nuN(E) ds
P_surv(E) = exp[-tau_DIS(E)]
```

Here `n_b` is the baryon number density and `sigma_nuN(E)` is the neutrino
nucleon DIS cross section. HADROS currently compares:

- `GBW`, a saturation-inspired DIS model table;
- `IIM`, another small-x saturation-inspired model table;
- `CTW_reference`, the central `nu N` charged-current table from Connolly,
  Thorne & Waters, Phys. Rev. D 83, 113009, Table I;
- `literature_powerlaw_scale`, an older approximate scale curve retained only
  for comparison.

`CTW_reference` is charged-current `nu N` only. It does not include neutral
current interactions, antineutrino tables, regeneration, or PDF uncertainty
bands.

![Optical depth concept](figures/manual_tau_concept.png)

### 2.4 UHE source morphologies

The UHE emissivity is separated into a source morphology and a spectral model.
The source morphologies are phenomenological:

- `inner_ring`: compact equatorial inner-ring baseline;
- `funnel_wall`: emission near the boundary of the low-density funnel;
- `jet_base`: compact source near the axis and small radii;
- `shock_layer` or density-gradient source: emission where density gradients
  are large;
- `density_weighted`: distributed source with
  `j_UHE proportional to rho^q r^-s`.

These models are not first-principles particle acceleration simulations. They
are controlled prescriptions for testing how source placement changes the
observed opacity maps.

### 2.5 UHE spectral models

The UHE spectral interface currently supports:

```text
monochromatic
powerlaw
powerlaw_cutoff
```

The power-law spectrum is

```text
dN/dE proportional to E^-gamma
```

and the cutoff power law is

```text
dN/dE proportional to E^-gamma exp(-E/Ecut)
```

The default remains `SPECTRAL_MODEL=monochromatic` for backward compatibility.
Energy-band composite images use false colors: blue for the low-energy UHE
band, green for the intermediate band, and red for the high-energy band. These
colors encode energy intervals; they are not physical photon colors.

### 2.6 Optical-depth surfaces

HADROS can extract diagnostic opacity surfaces such as:

```text
tau = 2/3
tau = 1
tau = 3
```

These surfaces are properties of the medium, geometry, neutrino energy, and DIS
cross section. They do not depend on the UHE source morphology. The current
implementation is axisymmetric and extracts `r_tau(theta)`. It is designed so a
future `r_tau(theta, phi)` extension can be added for true 3D backgrounds.

### 2.7 Thermal MeV neutrino module

The MeV module is separate from UHE/DIS. It is a diagnostic local emissivity and
radiative-transfer model, not a calibrated absolute luminosity prediction.

Implemented channels include URCA-like processes:

```text
e- + p -> n + nu_e
e+ + n -> p + anti_nu_e
```

pair annihilation:

```text
e- + e+ -> nu + anti-nu
```

and nucleon-nucleon bremsstrahlung:

```text
N + N -> N + N + nu + anti-nu
```

The MeV opacity includes absorption and scattering proxies with approximate
energy scaling like `sigma ~ E_MeV^2`. The local transfer equation is:

```text
dI/ds = j - alpha I
```

where `j` is emissivity and `alpha` is opacity. CPU MeV is the canonical
physical reference implementation. CUDA MeV remains legacy/toy until it is
ported and validated.

![UHE and MeV modules](figures/manual_uhe_vs_mev_physics.png)

### 2.8 MeV diagnostic upgrades

The MeV diagnostic workflow includes temperature profiles, electron-fraction
profiles, Fermi-Dirac-like spectral integration, MeV neutrinosphere extraction
at `tau_MeV = 2/3`, electron-degeneracy diagnostics, opacity-component
decomposition, and luminosity-proxy audits.

Luminosities from this module should be described as diagnostic proxy values
unless a separate calibration against literature or simulations is performed.

### 2.9 ParaView and 3D visualization

The ParaView export samples the current axisymmetric semi-analytic model onto a
3D Cartesian grid. It can export density, log-density, UHE emissivity,
log-emissivity, radius, polar angle, azimuth, and normalized source morphology.

This is a 3D Cartesian sampling of an axisymmetric model. It is not a fully 3D
hydrodynamical simulation. HADROS does not export a fake local tau field;
opacity requires a line-of-sight or radial integration convention.

## 3. Main Makefile Commands

### 3.1 Safe build and validation

Show available targets:

```bash
make
```

Compile CPU binaries:

```bash
make build NTHREADS=2
```

Run a cheap CPU validation:

```bash
make validate_small NTHREADS=2
```

Run a tiny CPU/CUDA smoke test:

```bash
make validate_small_cuda NTHREADS=2 \
  NVCC=/home/rafael/micromamba/envs/dis/bin/nvcc \
  PYTHON=/home/rafael/micromamba/envs/dis/bin/python
```

Run the medium CPU/CUDA validation:

```bash
make validate_medium_cuda NTHREADS=2 \
  NVCC=/home/rafael/micromamba/envs/dis/bin/nvcc \
  PYTHON=/home/rafael/micromamba/envs/dis/bin/python
```

### 3.2 Running images

Monochromatic UHE image using a cache:

```bash
make image-from-cache ENU=1e11 DENSITY_PROFILE=gaussian NTHREADS=2
```

Spectral UHE image:

```bash
make image-from-cache SPECTRAL_MODEL=powerlaw \
  SPECTRAL_GAMMA=2.0 SPECTRAL_E_MIN_GEV=1e5 \
  SPECTRAL_E_MAX_GEV=1e12 SPECTRAL_N_BINS=8 NTHREADS=2
```

Energy-band composite UHE image:

```bash
make plot_multiband_image PYTHON=python3
```

MeV physical image:

```bash
make image-from-cache MEV_MODEL=physical MEV_FLAVOR=anti_nu_e NTHREADS=2
```

MeV multiband image:

```bash
make mev_multiband_image NTHREADS=2 PYTHON=python3
```

### 3.3 Running robustness scans

```bash
make robustness_scans NTHREADS=2 PYTHON=python3
```

Outputs go to `output/scans/` and `plots/scans/`.

### 3.4 DIS validation

```bash
make validate_dis_models NTHREADS=2 PYTHON=python3
```

This compares GBW, IIM, `literature_powerlaw_scale`, and `CTW_reference`.
Outputs go to `output/dis_validation/` and `plots/dis_validation/`.

### 3.5 Opacity surfaces

```bash
make opacity_surfaces TAU_SURFACE_VALUE=1.0 PYTHON=python3
```

Outputs go to `output/opacity_surfaces/` and `plots/opacity_surfaces/`.

### 3.6 MeV diagnostics

```bash
make validate_mev_physics PYTHON=python3
make mev_neutrinosphere PYTHON=python3
make mev_luminosity PYTHON=python3
make audit_mev_luminosity PYTHON=python3
make audit_torus_regime PYTHON=python3
make audit_collapsar_ndaf_like PYTHON=python3
```

### 3.7 ParaView export

```bash
make paraview_fields \
  DENSITY_PROFILE=powerlaw_funnel_envelope \
  SOURCE_MODEL=funnel_wall \
  PARAVIEW_NX=64 PARAVIEW_NY=64 PARAVIEW_NZ=64 \
  PARAVIEW_BOX_RG=80 NTHREADS=4
```

Output:

```text
output/paraview/bh_torus_fields.vtk
```

### 3.8 Dashboard

```bash
make dashboard PYTHON=python3
```

Open:

```text
dashboard/index.html
```

The dashboard indexes existing plots and output products. It does not run
simulations.

### 3.9 Thread control

HADROS defaults to safe thread use:

```text
NTHREADS ?= 4
OMP_NUM_THREADS ?= $(NTHREADS)
```

Use explicit limits:

```bash
make -j2 build
make validate_small NTHREADS=2
make run_production NTHREADS=4
```

Production runs are explicit. The default `make` target prints help and does
not launch a full simulation.

## 4. Outputs and Directory Map

![Directory map](figures/manual_directory_map.png)

- `data/`: input data tables, especially DIS sigma tables.
- `include/`: C++ headers.
- `src/`: core C++ implementations.
- `apps/`: command-line applications built by the Makefile.
- `scripts/`: plotting, validation, and audit scripts.
- `docs/`: manuals and scientific notes.
- `plots/`: generated figures.
- `output/`: generated tables, images, ray caches, and validation reports.
- `dashboard/`: static HTML dashboard.

Common places to look:

- Images: `output/images/`
- Plot figures: `plots/`
- DIS reports: `output/dis_validation/`
- MeV reports: `output/validation/`, `output/mev_luminosity/`,
  `output/mev_neutrinosphere/`
- ParaView files: `output/paraview/`
- Dashboard: `dashboard/index.html`

## 5. Program Workflow

The internal workflow is:

1. Choose a physical setup: density profile, source model, spectrum, observer,
   spin, and energy.
2. Generate or load a Kerr ray cache.
3. Evaluate the local background along each ray.
4. Compute local UHE and/or MeV emissivity.
5. Compute opacity.
6. Integrate radiative transfer along rays.
7. Write data tables.
8. Generate plots.
9. Index products in the dashboard.

```text
Makefile variables
  -> executable arguments
  -> C++ model objects
  -> ray-by-ray integration
  -> output/*.dat and output/*.csv
  -> scripts/*.py plots
  -> dashboard/index.html
```

## 6. Recommended Workflows for Students

### Beginner: run dashboard only

```bash
make dashboard PYTHON=python3
```

Open `dashboard/index.html` and explore existing products.

### Beginner: run small validation

```bash
make validate_small NTHREADS=2 PYTHON=python3
```

This checks that the code compiles, creates a small ray cache, and runs a cheap
image calculation.

### Intermediate: compare DIS models

```bash
make validate_dis_models NTHREADS=2 PYTHON=python3
```

Read `output/dis_validation/dis_model_validation_summary.md`.

### Intermediate: generate MeV diagnostic plots

```bash
make validate_mev_physics PYTHON=python3
```

Then inspect `plots/mev_physics/`.

### Advanced: create collapsar_ndaf_like comparison

```bash
make audit_collapsar_ndaf_like PYTHON=python3
```

This compares the baseline and literature-guided semi-analytic preset without
changing the physics modules.

### Advanced: export to ParaView

```bash
make paraview_fields DENSITY_PROFILE=collapsar_ndaf_like NTHREADS=4
```

Open `output/paraview/bh_torus_fields.vtk` in ParaView.

## 7. Scientific Limitations

- Backgrounds are semi-analytic, not hydrodynamic evolutions.
- HADROS does not solve full neutrino radiation hydrodynamics.
- The MeV luminosities are proxy/diagnostic quantities, not calibrated
  absolute predictions.
- CUDA MeV remains legacy/toy and is not equivalent to CPU physical MeV.
- `CTW_reference` is central charged-current `nu N` only.
- UHE neutral-current interactions and regeneration are not yet included.
- There is no true 3D hydrodynamic background yet.
- ParaView export samples axisymmetric fields into 3D; it does not create a
  fully 3D simulation.

## 8. Troubleshooting

If `make` seems to run too long, stop the process and use a validation target
with smaller resolution, such as `validate_small`.

If CPU usage is too high, pass `NTHREADS=2` or `NTHREADS=4`, and use limited
make parallelism such as `make -j2 build`.

If CUDA is not found, check:

```bash
which nvcc
nvcc --version
nvidia-smi
```

Then pass the explicit `NVCC` path if needed.

If plots are missing, regenerate the relevant plotting target, then rebuild the
dashboard with `make dashboard PYTHON=python3`.

If stale caches are suspected after a crash, do not trust old ray or image
caches automatically. See `README_RECOVERY_20260609.md`.

If `NaN` or `Inf` appears in outputs, check density floors, energy ranges, and
whether the requested energy is inside the sigma-table range.

If output files become huge, avoid production targets and use validation
resolutions first.

## 9. Glossary

DIS: deep inelastic scattering, the high-energy neutrino-nucleon interaction
used for UHE opacity.

UHE: ultra-high energy, here referring to GeV-to-EeV scale neutrinos.

MeV: mega-electronvolt; thermal neutrino energies in compact-object disks are
often in the MeV range.

Optical depth: dimensionless attenuation measure, `tau`.

Survival probability: `P_surv = exp(-tau)`.

Neutrinosphere: diagnostic surface where optical depth reaches a chosen value,
often `tau = 2/3`.

Kerr metric: spacetime geometry around a rotating black hole.

Ray tracing: integration of camera rays or geodesics through spacetime.

NDAF: neutrino-dominated accretion flow.

Collapsar: massive-star collapse scenario that can form a black hole and dense
accretion environment.

Emissivity: local production rate of radiation or neutrinos.

Opacity: local attenuation coefficient or interaction strength.

CTW: Connolly, Thorne & Waters; HADROS uses their central `nu N` CC table as a
published reference.

GBW: Golec-Biernat-Wusthoff saturation-inspired DIS model.

IIM: Iancu-Itakura-Munier saturation-inspired DIS model.
