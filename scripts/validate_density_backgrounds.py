"""Validate controlled semi-analytic density morphologies.

This script plots radial and 2D density diagnostics for the density fields used
as controlled semi-analytic backgrounds in the UHE opacity study. These profiles
are not hydrodynamical solutions or equilibrium accretion disks.
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


OUT_DIR = Path("plots/validation_density_backgrounds")

RHO0 = 1.0e10
R0 = 10.0
SIGMA_R = 5.0
H_OVER_R = 0.25
RADIAL_POWER = 2.5
FUNNEL_DEPLETION = 0.95
FUNNEL_THETA = np.deg2rad(20.0)
ENVELOPE_RHO0 = 1.0e6
ENVELOPE_ALPHA = 2.5
R_MIN = 4.0
R_MAX = 80.0
RHO_FLOOR = 1.0e-99


def uses_powerlaw(profile):
    return profile in {
        "powerlaw",
        "powerlaw_funnel",
        "powerlaw_envelope",
        "powerlaw_funnel_envelope",
    }


def uses_funnel(profile):
    return profile in {
        "gaussian_funnel",
        "powerlaw_funnel",
        "powerlaw_funnel_envelope",
    }


def uses_envelope(profile):
    return profile in {
        "gaussian_envelope",
        "powerlaw_envelope",
        "powerlaw_funnel_envelope",
    }


def gaussian_shape(r, theta):
    delta = theta - 0.5 * np.pi
    shape = np.exp(-((r - R0) / SIGMA_R) ** 2) * np.exp(-(delta / H_OVER_R) ** 2)
    mask = (r > 4.0) & (r < 18.0) & (np.abs(delta) < 0.45)
    return np.where(mask, shape, 0.0)


def powerlaw_shape(r, theta):
    in_domain = (r >= R_MIN) & (r <= R_MAX) & (r > 0.0)
    vertical = np.exp(-(np.cos(theta) / H_OVER_R) ** 2)
    inner_taper = 1.0 - np.exp(-((r - R_MIN) / SIGMA_R) ** 2)
    outer_taper = np.exp(-(r / R_MAX) ** 4)
    radial = (r / R0) ** (-RADIAL_POWER)
    return np.where(in_domain, radial * vertical * inner_taper * outer_taper, 0.0)


def funnel_factor(profile, theta):
    if not uses_funnel(profile):
        return np.ones_like(theta)

    north = np.exp(-(theta / FUNNEL_THETA) ** 2)
    south = np.exp(-((np.pi - theta) / FUNNEL_THETA) ** 2)
    return np.clip(1.0 - FUNNEL_DEPLETION * (north + south), 0.0, 1.0)


def envelope_rho(profile, r):
    if not uses_envelope(profile):
        return np.zeros_like(r)

    mask = (r >= R0) & (r <= R_MAX)
    rho = ENVELOPE_RHO0 * (r / R0) ** (-ENVELOPE_ALPHA)
    return np.where(mask, rho, 0.0)


def rho(profile, r, theta):
    shape = powerlaw_shape(r, theta) if uses_powerlaw(profile) else gaussian_shape(r, theta)
    raw = RHO0 * shape * funnel_factor(profile, theta) + envelope_rho(profile, r)
    return np.maximum(raw, RHO_FLOOR)


def plot_profile(profile):
    r = np.linspace(1.0, R_MAX, 900)
    rho_eq = rho(profile, r, np.full_like(r, 0.5 * np.pi))
    rho_pol = rho(profile, r, np.zeros_like(r))

    fig, ax = plt.subplots(figsize=(6.2, 4.0))
    ax.semilogy(r, rho_eq, label=r"$\theta=\pi/2$")
    ax.semilogy(r, rho_pol, label=r"$\theta=0$")
    ax.set_xlabel(r"$r/r_g$")
    ax.set_ylabel(r"$\rho\;[\mathrm{g\,cm^{-3}}]$")
    ax.set_title(profile.replace("_", r"\_"))
    ax.grid(True, which="both", alpha=0.25)
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(OUT_DIR / f"{profile}_radial_profiles.png", dpi=220)
    plt.close(fig)


def plot_map(profile):
    r = np.linspace(1.0, R_MAX, 420)
    theta = np.linspace(0.0, np.pi, 300)
    rr, tt = np.meshgrid(r, theta)
    x = rr * np.sin(tt)
    z = rr * np.cos(tt)
    dens = rho(profile, rr, tt)

    fig, ax = plt.subplots(figsize=(6.2, 4.8))
    im = ax.pcolormesh(
        x,
        z,
        np.log10(dens),
        shading="auto",
        cmap="magma",
        vmin=-2,
        vmax=np.log10(RHO0),
    )
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel(r"$R/r_g$")
    ax.set_ylabel(r"$z/r_g$")
    ax.set_title(profile.replace("_", r"\_"))
    fig.colorbar(im, ax=ax, label=r"$\log_{10}\rho\;[\mathrm{g\,cm^{-3}}]$")
    fig.tight_layout()
    fig.savefig(OUT_DIR / f"{profile}_density_map.png", dpi=240)
    plt.close(fig)


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    profiles = [
        "gaussian",
        "powerlaw",
        "powerlaw_funnel",
        "powerlaw_funnel_envelope",
    ]

    for profile in profiles:
        plot_profile(profile)
        plot_map(profile)

    print(f"Saved validation plots to {OUT_DIR}")


if __name__ == "__main__":
    main()
