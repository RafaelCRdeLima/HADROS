# Torus Regime Literature Review

This note collects representative values used to classify the current
semi-analytic torus. The purpose is not to tune the model toward a desired
luminosity, but to understand which physical regime it currently resembles.

## Comparison Table

| Class | Representative density | Representative temperature | Ye | Disk mass | Radius / size | Neutrino luminosity | Accretion rate | References |
|---|---:|---:|---:|---:|---|---:|---:|---|
| NDAF / hyperaccreting BH disk | high-density disk; representative audit band `1e10-1e13 g cm^-3` | hot disk; representative audit band `5-20 MeV` | inner neutronization can drive `Ye < 0.1`; broader disk values can approach `~0.5` | model-dependent, often `1e-3-1 Msun` in engine-scale estimates | inner tens of gravitational radii | commonly discussed around `1e51-1e53 erg/s` for luminous neutrino-cooled systems | `~0.001-10 Msun/s` | Liu, Gu & Zhang 2017; Liu et al. 2007 |
| Collapsar hyperaccretion / NDAF | massive-core fallback can exceed NDAF ignition rate during early accretion | MeV spectra peak around `10-20 MeV` in collapsar NDAF calculations | model-dependent; neutron-rich inner zones expected when electron degeneracy/neutronization matter | supplied by stellar fallback; engine disk/envelope problem | BH plus hyperaccretion disk in stellar core | expected NDAF-scale neutrino emission when above ignition | mass supply above NDAF ignition; review-scale `~0.001-10 Msun/s` | Wei, Liu & Song 2019; Liu, Gu & Zhang 2017 |
| BNS merger remnant / disk | remnant can reach nuclear-density regimes; disks can be massive and optically thick depending on remnant | remnant matter can reach several tens of MeV, with maxima near `~50 MeV` in some simulations | neutron rich; disk/remnant values depend on EOS and remnant type | disk/remnant masses are simulation-dependent; BH-NS example below gives `0.3 Msun` disk | compact remnant/disk | total neutrino luminosity `~3-8e53 erg/s` in one long-lived HMNS class | transient post-merger accretion | Sekiguchi et al. 2011; Perego, Bernuzzi & Radice 2019 |
| BH-NS merger disk | hot compact nuclear disk | initially `T ~ 6 MeV` in one high-spin BH-NS simulation | average disk `Ye` rises to `~0.2` then decreases | about `0.3 Msun` disk plus `0.08 Msun` ejecta in the cited case | compact accretion disk | initially `L_nu ~ 1e54 erg/s`, falling by about an order of magnitude over `50 ms` | transient hyperaccretion | Deaton et al. 2013 |

## Source Notes

- Liu, Gu & Zhang review NDAFs as GRB engines and quote hyperaccretion rates
  around `0.001-10 Msun/s`. They emphasize high density and temperature, thick
  geometry, and the role of neutrino luminosity/annihilation.
- Liu, Gu, Xue & Lu compute neutrino-cooled disk structure and show that
  electron degeneracy and neutronization can reduce `Ye` below `0.1` in the
  inner disk.
- Wei, Liu & Song study collapsar black-hole hyperaccretion and report
  time-integrated MeV neutrino spectra with peaks around `10-20 MeV`.
- Sekiguchi et al. report total neutrino luminosity of a long-lived HMNS
  merger remnant around `3-8e53 erg/s`.
- Perego, Bernuzzi & Radice summarize merger thermodynamics: remnant matter can
  reach densities up to several times nuclear saturation and temperatures of
  order tens of MeV; BH disks are smaller and often more transparent, whereas
  massive-neutron-star remnants can be optically thick.
- Deaton et al. report a high-spin BH-NS case with about `0.3 Msun` in a hot
  compact disk, average `Ye ~ 0.2`, initial `T ~ 6 MeV`, and initial neutrino
  luminosity around `1e54 erg/s`.

## References

- Tong Liu, Wei-Min Gu & Bing Zhang, "Neutrino-dominated accretion flows as the
  central engine of gamma-ray bursts", arXiv:1705.05516,
  https://arxiv.org/abs/1705.05516
- Tong Liu, Wei-Min Gu, Li Xue & Ju-Fu Lu, "Structure and Luminosity of
  Neutrino-cooled Accretion Disks", arXiv:astro-ph/0702186,
  https://arxiv.org/abs/astro-ph/0702186
- Yun-Feng Wei, Tong Liu & Cui-Ying Song, "Black hole hyperaccretion in
  collapsars. I. MeV neutrinos", arXiv:1905.04850,
  https://arxiv.org/abs/1905.04850
- Yuichiro Sekiguchi, Kenta Kiuchi, Koutarou Kyutoku & Masaru Shibata,
  "Gravitational waves and neutrino emission from the merger of binary neutron
  stars", arXiv:1105.2125, https://arxiv.org/abs/1105.2125
- Albino Perego, Sebastiano Bernuzzi & David Radice, "Thermodynamics conditions
  of matter in neutron star mergers", arXiv:1903.07898,
  https://arxiv.org/abs/1903.07898
- M. Brett Deaton et al., "Black Hole-Neutron Star Mergers with a Hot Nuclear
  Equation of State", arXiv:1304.3384, https://arxiv.org/abs/1304.3384

## Interpretation Rule

The current model should be compared against these values as a controlled
semi-analytic background. Agreement in one variable, such as peak density, is
not enough to claim an NDAF/collapsar/merger regime if temperature, mass,
emitting volume, accretion state, and calibrated luminosity are not also in the
same physical range.

## Semi-Analytic Collapsar/NDAF-like Preset

The project includes a second controlled preset:

```text
DENSITY_PROFILE=collapsar_ndaf_like
MEV_THERMAL_PROFILE=collapsar_inner_hot
MEV_YE_PROFILE=collapsar_neutron_rich
```

Regime names are deliberately explicit:

```text
fiducial_uhe_default
```

Historical low-density UHE opacity fiducial, with `rho0 = 1e-2 g cm^-3`.

```text
fiducial_mev_density
```

Weakly neutrino-cooled MeV luminosity-audit fiducial, with
`rho0 = 1e10 g cm^-3`.

```text
collapsar_ndaf_like
```

The literature-guided semi-analytic collapsar/NDAF-like preset. These regimes
should not all be called simply `fiducial`.

The diagnostic comparison target uses literature-guided parameters rather than
luminosity tuning:

```text
rho0 = 3e10 g cm^-3
r0 = 12 rg
sigma_r = 8 rg
H/R = 0.45
T0 = 18 MeV
Ye_torus = 0.20
Ye_funnel = 0.42
```

The resulting audited values are:

```text
M_torus = 7.63e-2 Msun
rho_max = 2.83e11 g cm^-3
rho_emissivity_weighted_median = 1.09e11 g cm^-3
T_max = 18.2 MeV
T_emissivity_weighted_median = 15.0 MeV
Ye_range = 0.15-0.42
```

This naturally moves the semi-analytic background into the thermodynamic range
associated with collapsar/NDAF-like disks in the literature review, while still
remaining a post-processing morphology and not a hydrodynamical simulation.
