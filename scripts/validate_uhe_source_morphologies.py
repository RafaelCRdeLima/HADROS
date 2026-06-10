"""Generate diagnostic plots for phenomenological UHE source prescriptions."""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


PI = np.pi


def gaussian_density(r, theta, rho0=1.0e-2, r0=10.0, sigma=5.0, h=0.25):
    delta = theta - 0.5 * PI
    shape = np.exp(-((r - r0) / sigma) ** 2) * np.exp(-(delta / h) ** 2)
    mask = (r > 4.0) & (r < 18.0) & (np.abs(delta) < 0.45)
    return rho0 * shape * mask


def powerlaw_funnel_envelope_density(
    r,
    theta,
    rho0=1.0e-2,
    r0=10.0,
    sigma=5.0,
    h=0.25,
    radial_power=2.0,
    funnel_depletion=1.0,
    funnel_theta_deg=20.0,
    envelope_rho0=1.0e-4,
    envelope_alpha=2.5,
    r_min=4.0,
    r_max=60.0,
    rho_floor=1.0e-99,
):
    in_domain = (r >= r_min) & (r <= r_max)
    vertical = np.exp(-(np.cos(theta) / h) ** 2)
    inner_taper = 1.0 - np.exp(-((r - r_min) / sigma) ** 2)
    outer_taper = np.exp(-(r / r_max) ** 4)
    radial = np.power(np.maximum(r / r0, 1.0e-300), -radial_power)
    disk = rho0 * radial * vertical * inner_taper * outer_taper * in_domain

    theta_f = np.deg2rad(funnel_theta_deg)
    north = np.exp(-(theta / theta_f) ** 2)
    south = np.exp(-((PI - theta) / theta_f) ** 2)
    funnel = np.clip(1.0 - funnel_depletion * (north + south), 0.0, 1.0)

    envelope = np.where(
        (r >= r0) & (r <= r_max),
        envelope_rho0 * np.power(np.maximum(r / r0, 1.0e-300), -envelope_alpha),
        0.0,
    )

    return np.maximum(disk * funnel + envelope, rho_floor)


def bipolar_gaussian(theta, theta0, width):
    return np.exp(-((theta - theta0) / width) ** 2) + np.exp(
        -((theta - (PI - theta0)) / width) ** 2
    )


def source_weight(model, r, theta):
    source_r = 3.5
    sigma_r = 1.0
    theta_width = np.deg2rad(15.0)
    funnel_theta = np.deg2rad(20.0)
    density = powerlaw_funnel_envelope_density(r, theta)
    rho_norm = np.maximum(
        powerlaw_funnel_envelope_density(source_r, 0.5 * PI), 1.0e-300
    )

    if model == "inner_ring":
        return np.exp(-((r - source_r) / sigma_r) ** 2) * np.exp(
            -((theta - 0.5 * PI) / theta_width) ** 2
        )

    if model == "funnel_wall":
        return np.exp(-((r - source_r) / sigma_r) ** 2) * bipolar_gaussian(
            theta, funnel_theta, theta_width
        )

    if model == "jet_base":
        return np.exp(-((r / source_r) ** 2)) * bipolar_gaussian(
            theta, 0.0, theta_width
        )

    if model == "density_weighted":
        q = 1.0
        s = 2.0
        weight = np.power(np.maximum(density / rho_norm, 0.0), q) * np.power(
            np.maximum(r / source_r, 1.0e-300), -s
        )
        return np.clip(weight, 0.0, 1.0e2)

    if model == "shock_layer":
        dr = 0.1
        dtheta = np.deg2rad(1.0)
        rho_r_plus = powerlaw_funnel_envelope_density(r + dr, theta)
        rho_r_minus = powerlaw_funnel_envelope_density(np.maximum(r - dr, 1.0e-6), theta)
        rho_t_plus = powerlaw_funnel_envelope_density(r, np.clip(theta + dtheta, 0.0, PI))
        rho_t_minus = powerlaw_funnel_envelope_density(r, np.clip(theta - dtheta, 0.0, PI))
        d_r = (rho_r_plus - rho_r_minus) / (2.0 * dr)
        d_theta = (rho_t_plus - rho_t_minus) / (2.0 * dtheta * np.maximum(r, 1.0e-6))
        gradient = np.sqrt(d_r * d_r + d_theta * d_theta) / rho_norm
        return np.clip(gradient * np.exp(-((r / source_r) ** 2)), 0.0, 1.0e2)

    raise ValueError(model)


def tag(model):
    return model.replace("_", "-")


def main():
    outdir = Path("plots/validation_uhe_sources")
    outdir.mkdir(parents=True, exist_ok=True)

    r = np.linspace(1.5, 30.0, 360)
    theta = np.linspace(0.0, PI, 240)
    rr, tt = np.meshgrid(r, theta)
    xx = rr * np.sin(tt)
    zz = rr * np.cos(tt)

    models = [
        "inner_ring",
        "funnel_wall",
        "jet_base",
        "shock_layer",
        "density_weighted",
    ]

    for model in models:
        weight = source_weight(model, rr, tt)
        norm = np.max(weight)
        if norm > 0.0:
            weight = weight / norm

        fig, ax = plt.subplots(figsize=(6.2, 4.8), constrained_layout=True)
        im = ax.pcolormesh(xx, zz, np.log10(weight + 1.0e-12), shading="auto")
        ax.set_aspect("equal", adjustable="box")
        ax.set_xlabel(r"$R_{\rm cyl}/r_g$")
        ax.set_ylabel(r"$z/r_g$")
        ax.set_title(f"UHE source morphology: {model}")
        cbar = fig.colorbar(im, ax=ax)
        cbar.set_label(r"$\log_{10}(S/S_{\max} + 10^{-12})$")
        fig.savefig(outdir / f"uhe_source_{tag(model)}_map.png", dpi=180)
        plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.0), constrained_layout=True)
    r_line = np.linspace(1.5, 30.0, 500)
    theta_eq = np.full_like(r_line, 0.5 * PI)
    theta_pole = np.zeros_like(r_line)

    for model in models:
        eq = source_weight(model, r_line, theta_eq)
        pole = source_weight(model, r_line, theta_pole)
        if np.max(eq) > 0.0:
            eq = eq / np.max(eq)
        if np.max(pole) > 0.0:
            pole = pole / np.max(pole)
        axes[0].plot(r_line, eq, label=model)
        axes[1].plot(r_line, pole, label=model)

    axes[0].set_title("Equatorial source profile")
    axes[1].set_title("Polar source profile")
    for ax in axes:
        ax.set_xlabel(r"$r/r_g$")
        ax.set_ylabel(r"$S/S_{\max}$")
        ax.set_yscale("log")
        ax.set_ylim(1.0e-8, 2.0)
        ax.grid(alpha=0.25)
    axes[1].legend(fontsize=7)
    fig.savefig(outdir / "uhe_source_radial_profiles.png", dpi=180)
    plt.close(fig)

    print(f"Saved source diagnostics to {outdir}")


if __name__ == "__main__":
    main()
