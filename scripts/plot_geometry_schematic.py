import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Wedge


def arg_float(index, default):
    if len(sys.argv) > index:
        return float(sys.argv[index])
    return default


def tag_float(value):
    return f"{value:.2g}".replace("+", "_").replace("-", "m").replace(".", "p")


a_spin = arg_float(1, 0.0)
M_bh_msun = arg_float(2, 3.0)
torus_r0 = arg_float(3, 10.0)
torus_sigma = arg_float(4, 5.0)
torus_h_over_r = arg_float(5, 0.25)
source_r = arg_float(6, 3.5)
source_sigma = arg_float(7, 1.0)
source_theta_deg = arg_float(8, 15.0)
cam_theta_deg = arg_float(9, 75.0)

plot_dir = Path("plots")
plot_dir.mkdir(parents=True, exist_ok=True)

a_spin = np.clip(a_spin, 0.0, 0.999999)
r_h = 1.0 + np.sqrt(1.0 - a_spin * a_spin)

torus_inner = max(0.0, torus_r0 - torus_sigma)
torus_outer = torus_r0 + torus_sigma
source_inner = max(r_h, source_r - source_sigma)
source_outer = source_r + source_sigma

plot_rmax = max(1.25 * torus_outer, 1.25 * source_outer, 6.0)

fig, ax = plt.subplots(figsize=(8, 8))
ax.set_facecolor("#07090f")

torus = Wedge(
    (0.0, 0.0),
    torus_outer,
    0.0,
    360.0,
    width=max(torus_outer - torus_inner, 1.0e-6),
    facecolor="#1f77b4",
    edgecolor="#9dd7ff",
    alpha=0.32,
    linewidth=1.4,
)
ax.add_patch(torus)

torus_core = Circle(
    (0.0, 0.0),
    torus_r0,
    fill=False,
    linestyle="--",
    linewidth=1.2,
    edgecolor="#cfeeff",
    alpha=0.9,
)
ax.add_patch(torus_core)

source_ring = Wedge(
    (0.0, 0.0),
    source_outer,
    0.0,
    360.0,
    width=max(source_outer - source_inner, 1.0e-6),
    facecolor="#ffd84d",
    edgecolor="#fff6a8",
    alpha=0.72,
    linewidth=1.2,
)
ax.add_patch(source_ring)

source_core = Circle(
    (0.0, 0.0),
    source_r,
    fill=False,
    linestyle=":",
    linewidth=1.6,
    edgecolor="#fff8bf",
)
ax.add_patch(source_core)

horizon = Circle(
    (0.0, 0.0),
    r_h,
    facecolor="black",
    edgecolor="#ffffff",
    linewidth=1.3,
    zorder=5,
)
ax.add_patch(horizon)

ergosphere_eq = Circle(
    (0.0, 0.0),
    2.0,
    fill=False,
    linestyle="-.",
    linewidth=0.9,
    edgecolor="#999999",
    alpha=0.7,
)
ax.add_patch(ergosphere_eq)

def radial_arrow(radius, angle_deg, color, label, text_offset=0.0):
    angle = np.deg2rad(angle_deg)
    end = np.array([radius * np.cos(angle), radius * np.sin(angle)])
    ax.annotate(
        "",
        xy=end,
        xytext=(0.0, 0.0),
        arrowprops=dict(
            arrowstyle="->",
            color=color,
            lw=1.3,
            shrinkA=0,
            shrinkB=0,
        ),
    )
    text_pos = 1.06 * end
    normal = np.array([-np.sin(angle), np.cos(angle)])
    text_pos = text_pos + text_offset * normal
    ax.text(
        text_pos[0],
        text_pos[1],
        label,
        color=color,
        fontsize=9,
        ha="center",
        va="center",
        bbox=dict(facecolor="#07090f", edgecolor="none", alpha=0.75, pad=1.5),
    )


radial_arrow(source_r, 70.0, "#ffd84d", rf"$r_{{\rm UHE}}={source_r:.1f}\,r_g$")
radial_arrow(torus_r0, 135.0, "#8fd3ff", rf"$r_{{\rm torus}}={torus_r0:.1f}\,r_g$")
radial_arrow(torus_outer, 205.0, "#cfeeff", rf"$r_0+\sigma={torus_outer:.1f}\,r_g$")

ax.text(
    0.0,
    0.0,
    "BH",
    color="white",
    fontsize=13,
    fontweight="bold",
    ha="center",
    va="center",
    zorder=7,
)

ax.annotate(
    "UHE ring",
    xy=(source_r, 0.0),
    xytext=(0.45 * plot_rmax, 0.18 * plot_rmax),
    color="#ffd84d",
    fontsize=11,
    arrowprops=dict(arrowstyle="->", color="#ffd84d", lw=1.2),
    bbox=dict(facecolor="#07090f", edgecolor="#ffd84d", alpha=0.85, pad=3),
)

ax.annotate(
    "Dense torus",
    xy=(0.0, torus_r0),
    xytext=(-0.78 * plot_rmax, 0.58 * plot_rmax),
    color="#8fd3ff",
    fontsize=11,
    arrowprops=dict(arrowstyle="->", color="#8fd3ff", lw=1.2),
    bbox=dict(facecolor="#07090f", edgecolor="#8fd3ff", alpha=0.85, pad=3),
)

ax.annotate(
    "Kerr horizon",
    xy=(r_h / np.sqrt(2.0), r_h / np.sqrt(2.0)),
    xytext=(0.25 * plot_rmax, -0.42 * plot_rmax),
    color="white",
    fontsize=11,
    arrowprops=dict(arrowstyle="->", color="white", lw=1.2),
    bbox=dict(facecolor="#07090f", edgecolor="white", alpha=0.85, pad=3),
)

info = (
    rf"$a={a_spin:.3g}$" "\n"
    rf"$M={M_bh_msun:.2g}\,M_\odot$" "\n"
    rf"$H/R={torus_h_over_r:.2g}$" "\n"
    rf"$\sigma_{{\rm UHE}}={source_sigma:.1f}\,r_g$" "\n"
    rf"$\theta_{{\rm src}}={source_theta_deg:.1f}^\circ$" "\n"
    rf"$\theta_{{\rm cam}}={cam_theta_deg:.1f}^\circ$"
)
ax.text(
    0.96 * plot_rmax,
    -0.96 * plot_rmax,
    info,
    color="#d7dce8",
    fontsize=10,
    ha="right",
    va="bottom",
    bbox=dict(facecolor="#101522", edgecolor="#3e4a61", alpha=0.9, pad=6),
)

ax.set_xlim(-plot_rmax, plot_rmax)
ax.set_ylim(-plot_rmax, plot_rmax)
ax.set_aspect("equal", adjustable="box")
ax.set_xlabel(r"$x/r_g$")
ax.set_ylabel(r"$y/r_g$")
ax.grid(color="#2a3140", linestyle=":", linewidth=0.7, alpha=0.8)
ax.set_title("Equatorial-plane geometry schematic")

tag = (
    f"a_{tag_float(a_spin)}"
    f"_torusR_{tag_float(torus_r0)}"
    f"_sourceR_{tag_float(source_r)}"
    f"_camTheta_{tag_float(cam_theta_deg)}"
)
output = plot_dir / f"geometry_schematic_{tag}.png"

fig.tight_layout()
fig.savefig(output, dpi=250)
plt.close(fig)

print(f"Saved: {output}")
