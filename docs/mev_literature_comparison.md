# MeV Luminosity Literature Comparison

This note compares the current diagnostic MeV luminosity proxy with typical
values discussed for NDAF, collapsar, and compact-merger accretion systems.

The current code value is:

```text
L_proxy ~= 2.8e44 erg/s-proxy
```

This is not a calibrated luminosity. The current local emissivity coefficients
are physically motivated proxies, so the comparison below is a realism-gap
diagnostic rather than a claim of direct disagreement.

## Reference Regimes

| System | Typical thermodynamic regime | Typical MeV neutrino luminosity | Notes |
|---|---|---:|---|
| NDAF / collapsar hyperaccretion disk | Hot, dense neutrino-cooled disk; spectra peak around 10-20 MeV in collapsar NDAF calculations | often discussed around 1e51-1e53 erg/s depending on accretion state | Wei et al. discuss collapsar black-hole hyperaccretion and MeV neutrino spectra with peak energies around 10-20 MeV. |
| NDAF structure models | Density, temperature, electron fraction, optical depth, and electron degeneracy are coupled | neutrino cooling/annihilation depends strongly on degeneracy and Ye | Liu et al. emphasize that electron degeneracy and neutronization can lower Ye below 0.1 in inner disks and modify optical depths/luminosities. |
| Binary neutron-star merger remnant | Hot HMNS/disk environment with nuclear EOS and neutrino cooling | total neutrino luminosity about 3-8e53 erg/s in one GR simulation class | Sekiguchi et al. report this scale for long-lived HMNS merger remnants. |
| BH-NS merger disk | Hot compact disk, initially T of order several MeV | initially around 1e54 erg/s, decreasing by about an order of magnitude over tens of ms in one simulation | Deaton et al. report an initially hot disk and high neutrino luminosity for a high-spin case. |

## Sources

- Wei, Liu & Song, "Black hole hyperaccretion in collapsars. I. MeV neutrinos", arXiv:1905.04850, https://arxiv.org/abs/1905.04850
- Liu, Gu, Xue & Lu, "Structure and Luminosity of Neutrino-cooled Accretion Disks", arXiv:astro-ph/0702186, https://arxiv.org/abs/astro-ph/0702186
- Sekiguchi, Kiuchi, Kyutoku & Shibata, "Gravitational waves and neutrino emission from the merger of binary neutron stars", arXiv:1105.2125, https://arxiv.org/abs/1105.2125
- Deaton et al., "Black Hole-Neutron Star Mergers with a Hot Nuclear Equation of State", arXiv:1304.3384, https://arxiv.org/abs/1304.3384

## Comparison With Current Model

The current baseline differs from the literature regimes in four main ways:

1. The emissivity coefficients are diagnostic and not calibrated weak-interaction
   cooling rates.
2. The temperature and electron-fraction profiles are semi-analytic
   post-processing fields, not the result of thermal balance or beta
   equilibrium.
3. The density morphology is a controlled Gaussian/semi-analytic torus, not a
   hydrodynamic NDAF/collapsar or merger-disk solution.
4. The emitting volume is compact and the luminosity is dominated by a small
   region of the semi-analytic torus.

## Required Direction To Reach Literature Regimes

The controlled scans show that raising density and especially temperature
increases the proxy luminosity, as expected from URCA-like `T^6` and pair-like
`T^9` scaling. However, matching `1e51-1e53 erg/s` should not be done by
arbitrarily multiplying a normalization.

Physically meaningful steps would be:

- calibrate local emissivities against standard neutrino cooling formulae,
- use thermodynamic profiles from a published NDAF/collapsar/merger-disk model,
- include finite-temperature chemical potentials and blocking factors,
- add composition/nuclear state information beyond a prescribed `Ye`.

## Bottom Line

The current `~2.8e44 erg/s-proxy` value should be interpreted as a consequence
of the present diagnostic semi-analytic setup. It is not, by itself, evidence
that the ray tracing, volume conversion, or opacity code is broken.
