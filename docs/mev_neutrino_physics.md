# MeV Neutrino Physics Module

This module treats thermal MeV neutrinos in the torus as a separate
post-processing component from the UHE/DIS pipeline. The CPU implementation in
`include/mev_neutrino_physics.hpp` and `src/mev_neutrino_physics.cpp` is the
canonical reference implementation.

It is not full neutrino radiation hydrodynamics. It is not a calibrated
prediction of absolute neutrino luminosity. It is a controlled local emissivity
and opacity model used to compare thermal MeV structure with the UHE opacity
maps.

## Separation From UHE/DIS

The UHE channel uses DIS cross sections, UHE source morphologies, UHE spectral
models, and Kerr ray tracing. The MeV channel uses low-energy weak-interaction
processes:

- beta processes / URCA,
- electron-positron pair annihilation,
- nucleon-nucleon bremsstrahlung,
- charged-current absorption,
- neutral-current nucleon scattering.

The canonical CPU MeV physics is implemented in:

```text
include/mev_neutrino_physics.hpp
src/mev_neutrino_physics.cpp
```

The UHE/DIS cross sections and sigma tables are not modified by this module.

CUDA status:

```text
MEV_CPU_MODEL = physical or toy
MEV_CUDA_STATUS = legacy_toy_not_equivalent
```

The CUDA MeV path has not yet been ported to the canonical CPU physical module.
CUDA MeV outputs should therefore be treated as legacy/toy MeV diagnostics, not
as equivalent physical MeV calculations.

## Flavors

The supported flavor labels are:

```text
nu_e
anti_nu_e
nu_x
```

`nu_x` represents mu/tau neutrinos and antineutrinos as a grouped heavy-lepton
species.

## Emissivity Scalings

The physical MeV model uses approximate local spectral emissivity models. They
are physically motivated proxies, not full weak-interaction rates with
degeneracy or blocking corrections.

URCA / beta processes:

```text
j_urca proportional to rho T^6 f_flavor(E,T)
```

Electron capture contributes to `nu_e`; positron capture contributes to
`anti_nu_e`. The composition factors use `Ye` and `1-Ye`. `nu_x` receives no
URCA contribution.

Pair annihilation:

```text
j_pair proportional to T^9 f_flavor(E,T)
```

This is available for all flavors.

Nucleon-nucleon bremsstrahlung:

```text
j_brems proportional to rho^2 T^5 f_flavor(E,T)
```

This is available for all flavors.

The spectral shape is a simple thermal proxy:

```text
f(E,T) = E^2 / [exp(E/T) + 1] / T^3
```

The emitted quantity is stored as a local spectral emissivity proxy per MeV.
Do not interpret the normalization as a calibrated luminosity without a
separate calibration step.

The CPU interface exposes:

```text
mev_emissivity_urca(rho_gcm3, T_MeV, Ye, E_MeV, flavor)
mev_emissivity_pair(rho_gcm3, T_MeV, Ye, E_MeV, flavor)
mev_emissivity_brems(rho_gcm3, T_MeV, Ye, E_MeV, flavor)
mev_total_emissivity(rho_gcm3, T_MeV, Ye, E_MeV, params)
```

## Electron Degeneracy Diagnostic

The module includes a diagnostic estimate of the electron degeneracy parameter:

```text
eta_e = mu_e / T
```

with:

```text
n_e = rho Ye / m_u
p_F = hbar c (3 pi^2 n_e)^(1/3)
mu_e approximately sqrt(p_F^2 + m_e^2)
eta_e = mu_e / T
```

The approximation uses a zero-temperature Fermi momentum and ignores
finite-temperature chemical-potential corrections. It is useful for identifying
where electron degeneracy may matter, not for precision beta-equilibrium
calculations.

The CPU interface exposes:

```text
electron_number_density_cm3(rho_gcm3, Ye)
electron_fermi_momentum_MeV(rho_gcm3, Ye)
electron_chemical_potential_MeV(rho_gcm3, Ye, T_MeV)
electron_degeneracy_eta(rho_gcm3, Ye, T_MeV)
mev_urca_degeneracy_correction(rho_gcm3, Ye, T_MeV, flavor)
```

An optional bounded URCA-like correction can be enabled with:

```text
MEV_USE_DEGENERACY_CORRECTION=1
```

The default is:

```text
MEV_USE_DEGENERACY_CORRECTION=0
```

because the correction is approximate. The implemented factor is deliberately
bounded between 0.5 and 1.5 to avoid unphysical numerical spikes.

## Opacity

The low-energy MeV opacity uses the approximate weak cross-section scaling:

```text
sigma(E) = sigma0 (E_MeV / 1 MeV)^2
```

Absorption:

```text
alpha_abs = n_target sigma_abs(E)
```

Scattering:

```text
alpha_scat = n_b sigma_scat(E) composition_factor
```

The returned opacity unit is:

```text
cm^-1
```

The default constants are:

```text
sigma_abs0_cm2  = 9.6e-44
sigma_scat0_cm2 = 1.7e-44
```

The CPU interface exposes:

```text
mev_opacity_absorption_cm_inv(rho_gcm3, Ye, E_MeV, params)
mev_opacity_scattering_cm_inv(rho_gcm3, Ye, E_MeV, params)
mev_total_opacity_cm_inv(rho_gcm3, Ye, E_MeV, params)
```

In the current approximation, `nu_x` has no charged-current absorption, but it
does have neutral-current scattering.

The opacity is decomposed into named diagnostic components:

```text
mev_opacity_absorption_neutron_cm_inv(...)
mev_opacity_absorption_proton_cm_inv(...)
mev_opacity_scattering_neutron_cm_inv(...)
mev_opacity_scattering_proton_cm_inv(...)
mev_opacity_scattering_electron_cm_inv(...)
```

The compatibility wrappers `mev_opacity_absorption_cm_inv` and
`mev_opacity_scattering_cm_inv` call the decomposed components. Channels can be
enabled or disabled with:

```text
MEV_INCLUDE_ABS_N
MEV_INCLUDE_ABS_P
MEV_INCLUDE_SCAT_N
MEV_INCLUDE_SCAT_P
MEV_INCLUDE_SCAT_E
```

Defaults are `1`. Channels that are not physically allowed for a flavor return
zero; for example, `nu_x` has neutral-current scattering but no charged-current
absorption in this approximation.

## Radiative Transfer

The MeV channel is integrated along the cached Kerr rays with the local transfer
equation:

```text
dI/ds = j - alpha I
```

For a path segment:

```text
d_tau = alpha ds
I_out = I_in exp(-d_tau) + (j/alpha) [1 - exp(-d_tau)]
```

The optically thin limit is handled explicitly:

```text
if alpha -> 0:
    I_out += j ds
```

This replaces the previous MeV image path in the CPU post-processing integrator,
which used a phenomenological leakage factor as the primary attenuation model.

## Diagnostic Temperature and Ye Backgrounds

The physical MeV module can use semi-analytic diagnostic temperature and
electron-fraction morphologies tied to the existing density background. These
are controlled post-processing fields, not hydrodynamic thermodynamic
solutions.

Available temperature profiles:

```text
constant
inner_hot_torus
radial_powerlaw
torus_plus_cool_envelope
collapsar_inner_hot
```

Available electron-fraction profiles:

```text
constant
neutron_rich_torus
funnel_proton_rich
torus_envelope_contrast
collapsar_neutron_rich
```

Main parameters:

```text
MEV_THERMAL_PROFILE
MEV_YE_PROFILE
MEV_T0_MEV
MEV_T_FLOOR_MEV
MEV_T_POWER
MEV_YE_TORUS
MEV_YE_FUNNEL
MEV_YE_ENVELOPE
MEV_YE_FLOOR
MEV_YE_CEIL
```

The code clamps `T` to a positive floor and clamps `Ye` to the configured
physical interval. These safeguards prevent non-finite emissivities and
opacities in low-density regions.

## MeV Spectral Integration

The default remains monochromatic for backward compatibility:

```text
MEV_SPECTRAL_MODE=monochromatic
```

For diagnostic band images, the CPU integrator also supports:

```text
MEV_SPECTRAL_MODE=fermi_dirac_band
```

with a Fermi-Dirac-like weight:

```text
w(E,T) = E^2 / [exp(E/T) + 1]
```

The band is controlled by:

```text
MEV_E_MIN_MEV
MEV_E_MAX_MEV
MEV_N_BINS
```

This is a spectral weighting for morphology and attenuation diagnostics. It is
not a calibrated absolute neutrino luminosity model.

The multiband diagnostic image uses false colors:

```text
blue  = low-energy MeV band, 3-8 MeV
green = intermediate-energy MeV band, 8-20 MeV
red   = high-energy MeV band, 20-50 MeV
```

The colors encode MeV neutrino energy intervals, not physical photon colors.

## Diagnostic MeV Neutrinosphere

The target `make mev_neutrinosphere` extracts diagnostic MeV optical-depth
surfaces using the MeV weak absorption plus scattering opacity, not DIS. The
surface is defined by:

```text
tau_MeV(r_tau, theta, E) = tau_target
```

The current implementation is axisymmetric and integrates radially outward over
the semi-analytic diagnostic density, temperature, and `Ye` fields. It uses
interpolation between radial samples and marks directions with no crossing
instead of creating artificial surfaces.

Outputs:

```text
output/mev_neutrinosphere/mev_tau_surface_tau067_E10MeV.dat
output/mev_neutrinosphere/mev_tau_surface_tau1_E10MeV.dat
output/mev_neutrinosphere/mev_tau_surface_tau3_E10MeV.dat
plots/mev_physics/mev_neutrinosphere_tau067.png
plots/mev_physics/mev_neutrinosphere_energy_dependence.png
plots/mev_physics/mev_tau_phase_diagram.png
plots/mev_physics/mev_vs_uhe_opacity_surfaces.png
```

The MeV-vs-UHE comparison is a structural diagnostic only: the two surfaces are
computed from different interaction physics and different energy regimes.

## Diagnostic Luminosity Proxy

The target `make mev_luminosity` computes a volume-integrated diagnostic MeV
energy-emissivity proxy:

```text
L_nu_proxy = integral q_nu dV
```

for:

```text
nu_e
anti_nu_e
nu_x
```

The current volume element is axisymmetric:

```text
dV = 2 pi r^2 sin(theta) dr dtheta
```

with `r` converted to physical centimeters using `r_g = GM/c^2`. The energy
integration uses the local spectral emissivity proxy and multiplies by the
neutrino energy in erg.

Outputs:

```text
output/mev_luminosity/mev_luminosity_summary.csv
output/mev_luminosity/mev_luminosity_summary.md
plots/mev_physics/mev_luminosity_vs_temperature.png
plots/mev_physics/mev_luminosity_vs_density.png
plots/mev_physics/mev_luminosity_channel_breakdown.png
```

NDAF/collapsar literature often discusses neutrino luminosities around
`1e51-1e53 erg/s` for high accretion-rate disks. The values here should still
be treated as diagnostic proxies unless the emissivity normalization is
calibrated against a physical model. Do not force agreement with literature
bands by arbitrary normalization.

For a complete accounting of the current baseline luminosity proxy, run:

```bash
make audit_mev_luminosity NTHREADS=2 PYTHON=python3
```

The audit writes:

```text
output/mev_luminosity/luminosity_budget.md
output/mev_luminosity/unit_audit.md
output/mev_luminosity/realism_gap_analysis.md
```

The literature context is summarized in:

```text
docs/mev_literature_comparison.md
```

## Makefile Parameters

```text
MEV_MODEL
MEV_INCLUDE_URCA
MEV_INCLUDE_PAIR
MEV_INCLUDE_BREMS
MEV_INCLUDE_ABSORPTION
MEV_INCLUDE_SCATTERING
MEV_FLAVOR
MEV_ENERGY_MEV
MEV_ENU
MEV_NORM
MEV_THERMAL_PROFILE
MEV_YE_PROFILE
MEV_T0_MEV
MEV_T_FLOOR_MEV
MEV_T_POWER
MEV_YE_TORUS
MEV_YE_FUNNEL
MEV_YE_ENVELOPE
MEV_YE_FLOOR
MEV_YE_CEIL
MEV_SPECTRAL_MODE
MEV_E_MIN_MEV
MEV_E_MAX_MEV
MEV_N_BINS
MEV_USE_DEGENERACY_CORRECTION
MEV_INCLUDE_ABS_N
MEV_INCLUDE_ABS_P
MEV_INCLUDE_SCAT_N
MEV_INCLUDE_SCAT_P
MEV_INCLUDE_SCAT_E
MEV_LUMINOSITY_NR
MEV_LUMINOSITY_NTH
MEV_LUMINOSITY_E_MIN_MEV
MEV_LUMINOSITY_E_MAX_MEV
MEV_LUMINOSITY_E_BINS
```

Safe defaults:

```text
MEV_MODEL=physical
MEV_INCLUDE_URCA=1
MEV_INCLUDE_PAIR=1
MEV_INCLUDE_BREMS=1
MEV_INCLUDE_ABSORPTION=1
MEV_INCLUDE_SCATTERING=1
MEV_FLAVOR=anti_nu_e
MEV_ENERGY_MEV=10
MEV_THERMAL_PROFILE=inner_hot_torus
MEV_YE_PROFILE=neutron_rich_torus
MEV_SPECTRAL_MODE=monochromatic
MEV_USE_DEGENERACY_CORRECTION=0
MEV_INCLUDE_ABS_N=1
MEV_INCLUDE_ABS_P=1
MEV_INCLUDE_SCAT_N=1
MEV_INCLUDE_SCAT_P=1
MEV_INCLUDE_SCAT_E=1
```

`MEV_MODEL=toy` is retained for comparison with the previous phenomenological
emissivity.

## Metadata

Image outputs using the CPU post-processing path record:

```text
MEV_MODEL
MEV_CPU_MODEL
MEV_CUDA_STATUS
MEV_FLAVOR
MEV_ENERGY_MEV
MEV_INCLUDE_URCA
MEV_INCLUDE_PAIR
MEV_INCLUDE_BREMS
MEV_INCLUDE_ABSORPTION
MEV_INCLUDE_SCATTERING
MEV_THERMAL_PROFILE
MEV_YE_PROFILE
MEV_T0_MEV
MEV_T_FLOOR_MEV
MEV_T_POWER
MEV_YE_TORUS
MEV_YE_FUNNEL
MEV_YE_ENVELOPE
MEV_YE_FLOOR
MEV_YE_CEIL
MEV_SPECTRAL_MODE
MEV_E_MIN_MEV
MEV_E_MAX_MEV
MEV_N_BINS
MEV_USE_DEGENERACY_CORRECTION
MEV_INCLUDE_ABS_N
MEV_INCLUDE_ABS_P
MEV_INCLUDE_SCAT_N
MEV_INCLUDE_SCAT_P
MEV_INCLUDE_SCAT_E
mev_emissivity_units
mev_opacity_units
```

## Validation

Run:

```bash
make validate_mev_physics NTHREADS=2 PYTHON=python3
```

Outputs:

```text
output/validation/mev_physics_validation.txt
plots/mev_physics/mev_emissivity_channels.png
plots/mev_physics/mev_opacity_vs_energy.png
plots/mev_physics/mev_channel_dominance_map.png
plots/mev_physics/mev_flavor_comparison.png
plots/mev_physics/mev_transfer_limits.png
plots/mev_physics/mev_image_physical_vs_toy.png
```

Additional realism-upgrade diagnostics:

```bash
make mev_multiband_image NTHREADS=2 PYTHON=python3
make mev_neutrinosphere PYTHON=python3
make validate_mev_degeneracy PYTHON=python3
make validate_mev_opacity_components PYTHON=python3
make mev_luminosity NTHREADS=2 PYTHON=python3
make validate_mev_upgrades NTHREADS=2 PYTHON=python3
```

Outputs:

```text
output/validation/mev_realism_upgrade_validation.txt
output/mev_neutrinosphere/mev_multiband_flux.csv
plots/mev_physics/mev_multiband_false_color_image.png
plots/mev_physics/mev_temperature_profiles.png
plots/mev_physics/mev_ye_profiles.png
plots/mev_physics/mev_fermi_dirac_spectral_weights.png
plots/mev_physics/mev_neutrinosphere_tau067.png
plots/mev_physics/mev_neutrinosphere_energy_dependence.png
plots/mev_physics/mev_tau_phase_diagram.png
plots/mev_physics/mev_vs_uhe_opacity_surfaces.png
plots/mev_physics/mev_electron_degeneracy_map.png
plots/mev_physics/mev_eta_vs_density_temperature.png
plots/mev_physics/mev_opacity_components_vs_energy.png
plots/mev_physics/mev_opacity_components_rhoT_map.png
plots/mev_physics/mev_luminosity_vs_temperature.png
plots/mev_physics/mev_luminosity_vs_density.png
plots/mev_physics/mev_luminosity_channel_breakdown.png
```

The validation checks finite outputs, positive emissivity and opacity, pair
temperature scaling, bremsstrahlung density scaling, flavor behavior, `nu_x`
charged-current absorption suppression, the alpha-to-zero limit, the optically
thin transfer limit, the optically thick source-function limit, difference from
the toy image, URCA dominance in dense hot electron-flavor conditions, and steep
pair temperature growth.

## Limitations

The module currently omits:

- a full finite-temperature electron-degeneracy treatment,
- finite-temperature chemical-potential corrections,
- nucleon degeneracy,
- final-state blocking factors,
- detailed weak magnetism corrections,
- full neutrino radiation hydrodynamics,
- calibrated absolute MeV luminosities,
- flavor oscillations,
- self-consistent thermal/radiation feedback.

The CUDA MeV path still uses the previous toy/leakage implementation and is not
part of the physical MeV validation. Use the CPU post-processing path for the
physical MeV model until the CUDA MeV kernel is updated and validated
separately.
