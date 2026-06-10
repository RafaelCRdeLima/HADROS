#!/usr/bin/env python3
"""Create lightweight schematic figures for the HADROS user manual."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, Rectangle


ROOT = Path(__file__).resolve().parents[1]
FIG_DIR = ROOT / "docs" / "figures"


def setup_ax(figsize=(10, 4)):
    fig, ax = plt.subplots(figsize=figsize)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    return fig, ax


def box(ax, xy, text, width=0.16, height=0.18, face="#eef4f8"):
    x, y = xy
    patch = Rectangle((x, y), width, height, ec="#1f3b57", fc=face, lw=1.4)
    ax.add_patch(patch)
    ax.text(x + width / 2, y + height / 2, text, ha="center", va="center", fontsize=9, wrap=True)
    return patch


def arrow(ax, start, end):
    ax.add_patch(
        FancyArrowPatch(
            start,
            end,
            arrowstyle="->",
            mutation_scale=14,
            lw=1.5,
            color="#38485a",
        )
    )


def manual_pipeline_scheme():
    fig, ax = setup_ax((11, 3.6))
    labels = [
        "UHE / MeV\nsource model",
        "Kerr camera\nand ray cache",
        "density, T, Ye\nbackground",
        "opacity and\nradiative transfer",
        "images, CSV,\nplots, dashboard",
    ]
    xs = [0.03, 0.23, 0.43, 0.63, 0.82]
    for x, label in zip(xs, labels):
        box(ax, (x, 0.42), label, width=0.15, height=0.2)
    for x1, x2 in zip(xs[:-1], xs[1:]):
        arrow(ax, (x1 + 0.15, 0.52), (x2, 0.52))
    ax.text(0.5, 0.82, "HADROS high-level workflow", ha="center", fontsize=13, weight="bold")
    fig.savefig(FIG_DIR / "manual_pipeline_scheme.png", dpi=180, bbox_inches="tight")
    plt.close(fig)


def manual_tau_concept():
    fig, ax = setup_ax((8, 4.2))
    ax.plot([0.1, 0.9], [0.25, 0.75], color="#1f3b57", lw=3)
    for t, rho in zip([0.2, 0.35, 0.5, 0.65, 0.8], [0.15, 0.35, 0.75, 0.45, 0.2]):
        x = 0.1 + 0.8 * t
        y = 0.25 + 0.5 * t
        ax.scatter([x], [y], s=800 * rho, color="#d75c37", alpha=0.35, edgecolor="none")
    arrow(ax, (0.12, 0.2), (0.88, 0.68))
    ax.text(0.5, 0.9, r"Optical depth along one ray", ha="center", fontsize=13, weight="bold")
    ax.text(0.5, 0.1, r"$\tau_{\rm DIS}(E)=\int n_b(s)\,\sigma_{\nu N}(E)\,ds$" "\n" r"$P_{\rm surv}=\exp(-\tau)$", ha="center", fontsize=13)
    ax.text(0.18, 0.62, "denser regions\ncontribute more", ha="center", fontsize=9)
    fig.savefig(FIG_DIR / "manual_tau_concept.png", dpi=180, bbox_inches="tight")
    plt.close(fig)


def manual_uhe_vs_mev_physics():
    fig, ax = setup_ax((10, 4.8))
    box(ax, (0.08, 0.58), "UHE neutrinos\nGeV to EeV\nDIS opacity", width=0.28, height=0.22, face="#e8f1fb")
    box(ax, (0.64, 0.58), "MeV neutrinos\nthermal diagnostic\nlocal j and alpha", width=0.28, height=0.22, face="#f8eee8")
    box(ax, (0.08, 0.18), r"$\tau=\int n_b\sigma_{\nu N}ds$" "\nsource morphologies\nspectral weights", width=0.28, height=0.24, face="#f5f9fc")
    box(ax, (0.64, 0.18), "URCA-like, pair,\nbremsstrahlung\n" r"$dI/ds=j-\alpha I$", width=0.28, height=0.24, face="#fcf7f4")
    arrow(ax, (0.22, 0.58), (0.22, 0.42))
    arrow(ax, (0.78, 0.58), (0.78, 0.42))
    ax.text(0.5, 0.88, "UHE and MeV modules are separate diagnostics", ha="center", fontsize=13, weight="bold")
    ax.text(0.5, 0.5, "same semi-analytic background\nnot the same opacity physics", ha="center", fontsize=10)
    fig.savefig(FIG_DIR / "manual_uhe_vs_mev_physics.png", dpi=180, bbox_inches="tight")
    plt.close(fig)


def manual_directory_map():
    fig, ax = setup_ax((10, 5.2))
    entries = [
        ("data/", "input sigma tables"),
        ("include/", "C++ headers"),
        ("src/", "core implementations"),
        ("apps/", "command-line executables"),
        ("scripts/", "plotting and validation"),
        ("docs/", "manuals and notes"),
        ("output/", "generated data products"),
        ("plots/", "generated figures"),
        ("dashboard/", "static HTML index"),
    ]
    ax.text(0.5, 0.94, "HADROS directory map", ha="center", fontsize=13, weight="bold")
    for idx, (name, desc) in enumerate(entries):
        row = idx // 3
        col = idx % 3
        x = 0.06 + col * 0.31
        y = 0.68 - row * 0.24
        box(ax, (x, y), f"{name}\n{desc}", width=0.25, height=0.16, face="#f2f5ef")
    fig.savefig(FIG_DIR / "manual_directory_map.png", dpi=180, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    manual_pipeline_scheme()
    manual_tau_concept()
    manual_uhe_vs_mev_physics()
    manual_directory_map()
    print(f"Wrote figures to {FIG_DIR.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
