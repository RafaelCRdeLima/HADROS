# MeV Diagnostic Workflow

This workflow uses the canonical CPU MeV physics module for diagnostic and
morphological post-processing. It does not change the UHE/DIS opacity pipeline.

The products described here are not calibrated absolute neutrino luminosity
predictions. They are controlled diagnostics for comparing MeV thermal
structure, MeV opacity, and UHE opacity morphology.

## 1. Build

```bash
make build NTHREADS=2
```

## 2. Validate the CPU MeV Physics Module

```bash
make validate_mev_physics NTHREADS=2 PYTHON=python3
```

Main outputs:

```text
output/validation/mev_physics_validation.txt
plots/mev_physics/mev_emissivity_channels.png
plots/mev_physics/mev_opacity_vs_energy.png
plots/mev_physics/mev_channel_dominance_map.png
plots/mev_physics/mev_flavor_comparison.png
plots/mev_physics/mev_transfer_limits.png
plots/mev_physics/mev_image_physical_vs_toy.png
```

This validation checks positivity, finite outputs, temperature scaling, density
scaling, flavor behavior, and radiative-transfer limits.

## 3. Generate a MeV Energy-Band Composite Image

```bash
make mev_multiband_image NTHREADS=2 PYTHON=python3
```

Main outputs:

```text
plots/mev_physics/mev_multiband_false_color_image.png
output/mev_neutrinosphere/mev_multiband_flux.csv
```

The colors are false colors:

```text
blue  = low-energy MeV band, 3-8 MeV
green = intermediate-energy MeV band, 8-20 MeV
red   = high-energy MeV band, 20-50 MeV
```

They encode MeV neutrino energy intervals, not physical photon colors.

## 4. Generate Diagnostic MeV Neutrinospheres

```bash
make mev_neutrinosphere PYTHON=python3
```

Main outputs:

```text
output/validation/mev_realism_upgrade_validation.txt
output/mev_neutrinosphere/mev_tau_surface_tau067_E10MeV.dat
output/mev_neutrinosphere/mev_tau_surface_tau1_E10MeV.dat
output/mev_neutrinosphere/mev_tau_surface_tau3_E10MeV.dat
plots/mev_physics/mev_temperature_profiles.png
plots/mev_physics/mev_ye_profiles.png
plots/mev_physics/mev_fermi_dirac_spectral_weights.png
plots/mev_physics/mev_neutrinosphere_tau067.png
plots/mev_physics/mev_neutrinosphere_energy_dependence.png
plots/mev_physics/mev_tau_phase_diagram.png
plots/mev_physics/mev_vs_uhe_opacity_surfaces.png
```

The MeV tau surfaces use MeV weak absorption plus scattering opacity. They do
not use DIS cross sections.

## 5. Regenerate the Plot Dashboard

## 5. Run Next-Generation MeV Diagnostics

Electron degeneracy maps:

```bash
make validate_mev_degeneracy PYTHON=python3
```

Outputs:

```text
output/validation/mev_degeneracy_validation.txt
plots/mev_physics/mev_electron_degeneracy_map.png
plots/mev_physics/mev_eta_vs_density_temperature.png
```

Opacity component decomposition:

```bash
make validate_mev_opacity_components PYTHON=python3
```

Outputs:

```text
output/validation/mev_opacity_components_validation.txt
plots/mev_physics/mev_opacity_components_vs_energy.png
plots/mev_physics/mev_opacity_components_rhoT_map.png
```

Diagnostic luminosity proxy:

```bash
make mev_luminosity NTHREADS=2 PYTHON=python3
```

Outputs:

```text
output/validation/mev_luminosity_validation.txt
output/mev_luminosity/mev_luminosity_summary.csv
output/mev_luminosity/mev_luminosity_summary.md
plots/mev_physics/mev_luminosity_vs_temperature.png
plots/mev_physics/mev_luminosity_vs_density.png
plots/mev_physics/mev_luminosity_channel_breakdown.png
```

Run all three:

```bash
make validate_mev_upgrades NTHREADS=2 PYTHON=python3
```

These diagnostics improve physical interpretability but remain post-processing
proxies. In particular, the luminosity values should not be described as
calibrated absolute luminosities.

Audit the current luminosity proxy before interpreting its absolute scale:

```bash
make audit_mev_luminosity NTHREADS=2 PYTHON=python3
```

Main outputs:

```text
output/mev_luminosity/luminosity_budget.csv
output/mev_luminosity/luminosity_budget.md
output/mev_luminosity/unit_audit.md
output/mev_luminosity/realism_scan.csv
output/mev_luminosity/realism_gap_analysis.md
output/mev_luminosity/audit_summary.md
plots/mev_physics/emitting_volume_contribution.png
plots/mev_physics/emissivity_weighted_temperature_distribution.png
plots/mev_physics/emissivity_weighted_density_distribution.png
plots/mev_physics/luminosity_vs_rho0.png
plots/mev_physics/luminosity_vs_temperature.png
plots/mev_physics/luminosity_vs_Ye.png
```

The audit determines whether a low luminosity proxy comes from unit mistakes,
small emitting volume, low temperature/density, or the intended diagnostic
normalization.

Compare the named weakly cooled and collapsar/NDAF-like regimes:

```bash
make audit_collapsar_ndaf_like PYTHON=python3
```

The comparison uses three explicit names:

```text
fiducial_uhe_default
```

Historical low-density UHE opacity fiducial, `rho0 = 1e-2 g/cm3`.

```text
fiducial_mev_density
```

Weakly neutrino-cooled MeV luminosity-audit fiducial, `rho0 = 1e10 g/cm3`.

```text
collapsar_ndaf_like
```

Literature-guided semi-analytic collapsar/NDAF-like preset. Do not refer to all
three simply as `fiducial`.

Outputs:

```text
output/collapsar_ndaf_like/collapsar_statistics.md
output/collapsar_ndaf_like/collapsar_mass_estimate.md
output/collapsar_ndaf_like/comparison_summary.md
output/collapsar_ndaf_like/uhe_opacity_unit_audit.md
output/collapsar_ndaf_like/fiducial_mev_density_statistics.md
plots/collapsar_ndaf_like/mev_luminosity_comparison.png
plots/collapsar_ndaf_like/uhe_opacity_comparison.png
plots/collapsar_ndaf_like/neutrinosphere_comparison.png
```

The preset uses `DENSITY_PROFILE=collapsar_ndaf_like`,
`MEV_THERMAL_PROFILE=collapsar_inner_hot`, and
`MEV_YE_PROFILE=collapsar_neutron_rich` in a controlled diagnostic comparison.
It does not change emissivity normalization.

## 6. Regenerate the Plot Dashboard

```bash
make dashboard PYTHON=python3
```

Open:

```text
dashboard/index.html
```

The dashboard only indexes existing outputs; it does not run simulations.

## Parameters

Thermal morphology:

```text
MEV_THERMAL_PROFILE=constant|inner_hot_torus|radial_powerlaw|torus_plus_cool_envelope|collapsar_inner_hot
MEV_T0_MEV=6.0
MEV_T_FLOOR_MEV=0.1
MEV_T_POWER=0.2
```

Electron-fraction morphology:

```text
MEV_YE_PROFILE=constant|neutron_rich_torus|funnel_proton_rich|torus_envelope_contrast|collapsar_neutron_rich
MEV_YE_TORUS=0.25
MEV_YE_FUNNEL=0.55
MEV_YE_ENVELOPE=0.45
MEV_YE_FLOOR=0.01
MEV_YE_CEIL=0.60
```

Spectral mode:

```text
MEV_SPECTRAL_MODE=monochromatic|fermi_dirac_band
MEV_E_MIN_MEV=3.0
MEV_E_MAX_MEV=50.0
MEV_N_BINS=8
```

Next-generation options:

```text
MEV_USE_DEGENERACY_CORRECTION=0
MEV_INCLUDE_ABS_N=1
MEV_INCLUDE_ABS_P=1
MEV_INCLUDE_SCAT_N=1
MEV_INCLUDE_SCAT_P=1
MEV_INCLUDE_SCAT_E=1
MEV_LUMINOSITY_NR=96
MEV_LUMINOSITY_NTH=72
MEV_LUMINOSITY_E_MIN_MEV=1.0
MEV_LUMINOSITY_E_MAX_MEV=80.0
MEV_LUMINOSITY_E_BINS=24
```

## What Can Be Claimed

- The CPU MeV module provides a physically motivated local emissivity and
  opacity proxy.
- The diagnostic images show how MeV morphology changes with temperature,
  composition, opacity, and energy band.
- The MeV and UHE opacity surfaces can be compared as structural diagnostics of
  different interaction regimes.
- Electron degeneracy maps identify where eta_e may be non-negligible under a
  transparent zero-temperature chemical-potential approximation.
- The luminosity integral is a useful morphology and scaling proxy.

## What Should Not Be Claimed Yet

- Calibrated absolute MeV neutrino luminosities.
- Self-consistent neutrino radiation hydrodynamics.
- A physical CUDA MeV result equivalent to the CPU module.
- Fully realistic collapsar thermodynamics.
- Precision beta-equilibrium or finite-temperature chemical potentials.

CUDA MeV remains legacy/toy and is not equivalent to the canonical CPU physical
MeV module.
