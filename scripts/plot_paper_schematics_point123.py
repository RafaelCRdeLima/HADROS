"""Generate updated schematic Figures 1 and 2 for the paper.

The figures are conceptual summaries of the current controlled semi-analytic
model. They do not introduce new physics or replace the numerical ray-tracing
outputs.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
from matplotlib import patheffects
from matplotlib.colors import LogNorm
from matplotlib.patches import Circle, FancyArrowPatch, Polygon, Rectangle


ROOT = Path(__file__).resolve().parents[1]
PLOT_DIR = ROOT / "plots" / "paper_schematics"
PAPER_DIR = ROOT.parents[0] / "paper"


def density_powerlaw_funnel_envelope(r: np.ndarray, theta: np.ndarray) -> np.ndarray:
    rho0 = 1.0
    r0 = 10.0
    radial_power = 2.0
    h_over_r = 0.25
    r_min = 4.0
    r_max = 60.0
    rho_floor = 1.0e-8
    funnel_depletion = 0.98
    theta_f = np.deg2rad(20.0)
    envelope_rho0 = 1.0e-3
    envelope_alpha = 2.5

    theta_equator = np.pi / 2.0
    disk = rho0 * np.power(np.maximum(r, r_min) / r0, -radial_power)
    disk *= np.exp(-((np.cos(theta) / h_over_r) ** 2))
    disk = np.where((r >= r_min) & (r <= r_max), disk, 0.0)

    polar_angle = np.minimum(theta, np.pi - theta)
    funnel = 1.0 - funnel_depletion * np.exp(-((polar_angle / theta_f) ** 2))
    disk *= np.clip(funnel, 1.0e-4, 1.0)

    envelope = envelope_rho0 * np.power(np.maximum(r, r0) / r0, -envelope_alpha)
    envelope *= np.where(r >= r0, 1.0, 0.0)

    # Keep the equatorial disk visibly dominant while preserving the envelope.
    rho = disk + envelope + rho_floor
    return np.maximum(rho, rho_floor)


def label(ax, x, y, text, color="white", size=9, ha="center", va="center"):
    t = ax.text(x, y, text, color=color, fontsize=size, ha=ha, va=va)
    t.set_path_effects([patheffects.withStroke(linewidth=3, foreground="#0a0d12")])
    return t


def add_funnel_wall_band(ax, theta_deg: float, color: str) -> None:
    theta = np.deg2rad(theta_deg)
    width = np.deg2rad(3.0)
    r_inner = 4.2
    r_outer = 31.0

    for zsign in (-1, 1):
        for xsign in (-1, 1):
            th0 = theta - width
            th1 = theta + width
            points = []
            for r in np.linspace(r_inner, r_outer, 80):
                points.append((xsign * r * np.sin(th0), zsign * r * np.cos(th0)))
            for r in np.linspace(r_outer, r_inner, 80):
                points.append((xsign * r * np.sin(th1), zsign * r * np.cos(th1)))
            patch = Polygon(
                points,
                closed=True,
                facecolor=color,
                edgecolor="#f8ffff",
                linewidth=0.9,
                alpha=0.42,
                zorder=4,
            )
            ax.add_patch(patch)


def make_figure1() -> Path:
    x = np.linspace(-70.0, 70.0, 560)
    z = np.linspace(-70.0, 70.0, 560)
    X, Z = np.meshgrid(x, z)
    R = np.sqrt(X * X + Z * Z)
    theta = np.where(R > 0.0, np.arccos(np.clip(Z / np.maximum(R, 1.0e-30), -1.0, 1.0)), 0.0)

    rho = density_powerlaw_funnel_envelope(R, theta)
    fig, ax = plt.subplots(figsize=(7.2, 6.5), constrained_layout=True)
    fig.patch.set_facecolor("#f7f8fb")
    ax.set_facecolor("#080b10")

    im = ax.pcolormesh(
        X,
        Z,
        rho,
        shading="auto",
        cmap="magma",
        norm=LogNorm(vmin=1.0e-8, vmax=1.0),
    )

    # Funnel opening lines.
    theta_f = np.deg2rad(20.0)
    for sign in (-1, 1):
        ax.plot([0, sign * 70.0 * np.sin(theta_f)], [0, 70.0 * np.cos(theta_f)], color="#7fffd4", lw=1.2, ls="--", alpha=0.8)
        ax.plot([0, sign * 70.0 * np.sin(theta_f)], [0, -70.0 * np.cos(theta_f)], color="#7fffd4", lw=1.2, ls="--", alpha=0.8)

    add_funnel_wall_band(ax, 20.0, "#8fffea")

    # BH and characteristic radii.
    ax.add_patch(Circle((0, 0), 2.0, facecolor="black", edgecolor="white", lw=1.2, zorder=5))
    ax.add_patch(Circle((0, 0), 4.0, facecolor="none", edgecolor="#d6d6d6", lw=0.8, ls=":", alpha=0.9))
    ax.add_patch(Circle((0, 0), 10.0, facecolor="none", edgecolor="#ffd166", lw=1.0, ls="--", alpha=0.9))
    ax.add_patch(Circle((0, 0), 60.0, facecolor="none", edgecolor="#e5e7eb", lw=0.9, ls="--", alpha=0.55))

    label(ax, 0, 0, "BH", size=8)
    label(ax, 38, 9, "power-law disk / torus", color="#ffd166", size=9)
    label(ax, 19, 42, "depleted\npolar funnel", color="#7fffd4", size=9)
    label(ax, -50, -46, "external envelope", color="#f5c2ff", size=9)
    label(ax, -27, 19, "funnel-wall\nUHE source", color="#ccfff4", size=9)

    ax.annotate("", xy=(-14, 36), xytext=(-24, 23), arrowprops=dict(arrowstyle="->", color="#ccfff4", lw=1.3))
    ax.annotate("", xy=(36, 0), xytext=(39, 7), arrowprops=dict(arrowstyle="->", color="#ffd166", lw=1.3))
    ax.annotate("", xy=(-54, -31), xytext=(-50, -42), arrowprops=dict(arrowstyle="->", color="#f5c2ff", lw=1.3))

    cb = fig.colorbar(im, ax=ax, pad=0.02, shrink=0.88)
    cb.set_label(r"relative density $\rho/\rho_0$")

    ax.set_xlim(-70, 70)
    ax.set_ylim(-70, 70)
    ax.set_aspect("equal")
    ax.set_xlabel(r"$x/r_g$")
    ax.set_ylabel(r"$z/r_g$")
    ax.set_title("Controlled semi-analytic background sampled in a meridional plane")

    PLOT_DIR.mkdir(parents=True, exist_ok=True)
    output = PLOT_DIR / "paper_fig1_semi_analytic_background.png"
    fig.savefig(output, dpi=300)
    fig.savefig(PAPER_DIR / output.name, dpi=300)
    plt.close(fig)
    return output


def add_box(ax, xy, w, h, title, body, fc, ec):
    box = Rectangle(xy, w, h, facecolor=fc, edgecolor=ec, lw=1.3)
    ax.add_patch(box)
    ax.text(xy[0] + w / 2, xy[1] + h * 0.68, title, ha="center", va="center", fontsize=10, fontweight="bold", color="#111827")
    ax.text(xy[0] + w / 2, xy[1] + h * 0.34, body, ha="center", va="center", fontsize=8.3, color="#1f2937", linespacing=1.15)


def arrow(ax, start, end, color="#374151"):
    ax.add_patch(FancyArrowPatch(start, end, arrowstyle="-|>", mutation_scale=14, lw=1.3, color=color))


def make_figure2() -> Path:
    fig, ax = plt.subplots(figsize=(8.0, 5.4), constrained_layout=True)
    fig.patch.set_facecolor("#f8fafc")
    ax.set_facecolor("#f8fafc")
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 7)
    ax.axis("off")

    add_box(
        ax,
        (0.45, 4.65),
        2.1,
        1.05,
        "Kerr camera",
        "observer position\nspin, inclination\npixel grid",
        "#e0f2fe",
        "#0284c7",
    )
    add_box(
        ax,
        (3.0, 4.65),
        2.1,
        1.05,
        "Geodesic cache",
        "null rays\nproper path length\nredshift factors",
        "#e0e7ff",
        "#4f46e5",
    )
    add_box(
        ax,
        (5.55, 4.65),
        2.1,
        1.05,
        "Opacity integral",
        r"$d\tau=n_b\sigma_{\nu N}dl$" + "\nGBW / IIM tables",
        "#fee2e2",
        "#dc2626",
    )
    add_box(
        ax,
        (8.1, 4.65),
        1.45,
        1.05,
        "Image",
        "attenuated\nUHE map",
        "#fef3c7",
        "#d97706",
    )

    add_box(
        ax,
        (1.05, 1.2),
        2.25,
        1.4,
        "Density background",
        "power-law disk/torus\npolar funnel\nenvelope + floor",
        "#fce7f3",
        "#db2777",
    )
    add_box(
        ax,
        (3.85, 1.2),
        2.25,
        1.4,
        "UHE source",
        "inner ring\nfunnel wall\njet base / gradients",
        "#ecfccb",
        "#65a30d",
    )
    add_box(
        ax,
        (6.65, 1.2),
        2.25,
        1.4,
        "Diagnostics",
        r"$\tau$ maps, $P_{\rm surv}$" + "\nrobustness scans\nopacity surfaces",
        "#ede9fe",
        "#7c3aed",
    )

    arrow(ax, (2.55, 5.18), (3.0, 5.18))
    arrow(ax, (5.1, 5.18), (5.55, 5.18))
    arrow(ax, (7.65, 5.18), (8.1, 5.18))
    arrow(ax, (2.2, 2.6), (5.95, 4.65), "#db2777")
    arrow(ax, (4.95, 2.6), (8.25, 4.65), "#65a30d")
    arrow(ax, (6.7, 4.65), (7.75, 2.6), "#7c3aed")

    ax.text(
        5.0,
        6.45,
        "End-to-end synthetic UHE-neutrino image pipeline",
        ha="center",
        va="center",
        fontsize=13,
        fontweight="bold",
        color="#111827",
    )
    ax.text(
        5.0,
        0.50,
        "Opacity surfaces use only the medium and DIS cross section; image intensities also depend on source morphology.",
        ha="center",
        va="center",
        fontsize=8.7,
        color="#374151",
    )

    PLOT_DIR.mkdir(parents=True, exist_ok=True)
    output = PLOT_DIR / "paper_fig2_kerr_opacity_pipeline.png"
    fig.savefig(output, dpi=300)
    fig.savefig(PAPER_DIR / output.name, dpi=300)
    plt.close(fig)
    return output


def main() -> None:
    fig1 = make_figure1()
    fig2 = make_figure2()
    print(f"Saved: {fig1}")
    print(f"Saved: {PAPER_DIR / fig1.name}")
    print(f"Saved: {fig2}")
    print(f"Saved: {PAPER_DIR / fig2.name}")


if __name__ == "__main__":
    main()
