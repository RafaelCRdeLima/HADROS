import sys
import struct
import os
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.cm import ScalarMappable

# ==========================================================
# Parameters from Makefile
# Usage:
#   python scripts/plot_kerr_image.py 1e9 10.0 1.0 72.0
# ==========================================================

if len(sys.argv) > 1:
    Enu_GeV = float(sys.argv[1])
else:
    Enu_GeV = 1.0e9

if len(sys.argv) > 2:
    MeV_Enu = float(sys.argv[2])
else:
    MeV_Enu = 10.0

if len(sys.argv) > 3:
    MeV_Norm = float(sys.argv[3])
else:
    MeV_Norm = 1.0

if len(sys.argv) > 4:
    Cam_Theta_Deg = float(sys.argv[4])
else:
    Cam_Theta_Deg = 0.0

if len(sys.argv) > 5:
    Torus_Rho0 = float(sys.argv[5])
else:
    Torus_Rho0 = float(os.environ.get("TORUS_RHO0", "1.0e-2"))


def tag_sci(value):
    return f"{value:.0e}".replace("+", "_").replace("-", "m").replace(".", "p")


def tag_fixed(value):
    return f"{value:.1f}".replace("+", "_").replace("-", "m").replace(".", "p")


def energy_folder_tag(value):
    mantissa, exponent = f"{value:.0e}".split("e")
    return f"{mantissa}e{int(exponent)}"


def rho_folder_tag(value):
    mantissa, exponent = f"{value:.0e}".split("e")
    return f"rho0_torus_{mantissa}e{int(exponent)}"


energy_tag = tag_sci(Enu_GeV)
mev_energy_tag = tag_sci(MeV_Enu)
mev_norm_tag = tag_sci(MeV_Norm)
cam_theta_tag = tag_fixed(Cam_Theta_Deg)
rho_tag = rho_folder_tag(Torus_Rho0)
density_profile = os.environ.get("DENSITY_PROFILE", "").strip()
source_model = os.environ.get("SOURCE_MODEL", "").strip()
spectral_model = os.environ.get("SPECTRAL_MODEL", "").strip()

run_tag = (
    f"Enu_{energy_tag}"
    f"_MeVEnu_{mev_energy_tag}"
    f"_MeVNorm_{mev_norm_tag}"
    f"_CamTheta_{cam_theta_tag}"
)

profile_source_tag = ""
if density_profile and density_profile != "gaussian":
    profile_source_tag += f"_Profile_{density_profile}"
if source_model and source_model != "inner_ring":
    profile_source_tag += f"_Source_{source_model}"

spectrum_tag = ""
if spectral_model and spectral_model != "monochromatic":
    spectrum_tag = f"_Spectrum_{spectral_model}"

KIMG_MAGIC = 0x4B494D47
KIMG_HEADER = struct.Struct("<7i4d")


def load_kerr_image_binary(path):
    with open(path, "rb") as f:
        header_bytes = f.read(KIMG_HEADER.size)
        if len(header_bytes) != KIMG_HEADER.size:
            raise ValueError(f"Invalid binary image header: {path}")

        (
            magic,
            version,
            backend,
            nx_file,
            ny_file,
            ncols,
            reserved,
            enu_file,
            mev_enu_file,
            mev_norm_file,
            cam_theta_file,
        ) = KIMG_HEADER.unpack(header_bytes)

        if magic != KIMG_MAGIC:
            raise ValueError(f"Invalid binary image magic: {path}")
        if version != 1:
            raise ValueError(f"Unsupported binary image version {version}: {path}")
        if ncols < 8:
            raise ValueError(f"Invalid binary image column count {ncols}: {path}")

        payload = np.fromfile(f, dtype="<f8")

    expected = nx_file * ny_file * ncols
    if payload.size != expected:
        raise ValueError(
            f"Invalid binary image payload size in {path}: "
            f"expected {expected}, found {payload.size}"
        )

    return payload.reshape((nx_file * ny_file, ncols))


def load_kerr_image(path):
    if path.suffix == ".bin":
        return load_kerr_image_binary(path)

    return np.loadtxt(path)


input_candidates = [
    Path(f"output/images/kerr_image_cuda_cache_GBW_{rho_tag}{profile_source_tag}{spectrum_tag}_{run_tag}.bin"),
    Path(f"output/images/kerr_image_cuda_cache_GBW_{rho_tag}{profile_source_tag}{spectrum_tag}_{run_tag}.dat"),
    Path(f"output/images/kerr_image_cuda_cache_{rho_tag}{profile_source_tag}{spectrum_tag}_{run_tag}.bin"),
    Path(f"output/images/kerr_image_cuda_cache_{rho_tag}{profile_source_tag}{spectrum_tag}_{run_tag}.dat"),
    Path(f"output/images/kerr_image_cuda_cache_GBW_{rho_tag}{profile_source_tag}_{run_tag}.bin"),
    Path(f"output/images/kerr_image_cuda_cache_GBW_{rho_tag}{profile_source_tag}_{run_tag}.dat"),
    Path(f"output/images/kerr_image_cuda_cache_{rho_tag}{profile_source_tag}_{run_tag}.bin"),
    Path(f"output/images/kerr_image_cuda_cache_{rho_tag}{profile_source_tag}_{run_tag}.dat"),
    Path(f"output/images/kerr_image_cuda_cache_{rho_tag}_{run_tag}.bin"),
    Path(f"output/images/kerr_image_cuda_{rho_tag}_{run_tag}.bin"),
    Path(f"output/images/kerr_image_{rho_tag}_{run_tag}.bin"),
    Path(f"output/images/kerr_image_cuda_cache_{run_tag}.bin"),
    Path(f"output/images/kerr_image_cuda_{run_tag}.bin"),
    Path(f"output/images/kerr_image_{run_tag}.bin"),
    Path(f"output/images/kerr_image_cuda_cache_{rho_tag}_{run_tag}.dat"),
    Path(f"output/images/kerr_image_cuda_{rho_tag}_{run_tag}.dat"),
    Path(f"output/images/kerr_image_{rho_tag}_{run_tag}.dat"),
    Path(f"output/images/kerr_image_cuda_cache_{run_tag}.dat"),
    Path(f"output/images/kerr_image_cuda_{run_tag}.dat"),
    Path(f"output/images/kerr_image_{run_tag}.dat"),
    Path(f"output/images/kerr_image_Enu_{energy_tag}.dat"),
]

input_file = None
for candidate in input_candidates:
    if candidate.exists():
        input_file = candidate
        break

if input_file is None:
    raise FileNotFoundError(
        "Could not find Kerr image file. Tried: "
        + ", ".join(str(p) for p in input_candidates)
    )

if input_file.name == f"kerr_image_Enu_{energy_tag}.dat":
    run_tag = f"Enu_{energy_tag}"
elif input_file.name.startswith("kerr_image_cuda_cache_"):
    run_tag = f"cuda_cache_{run_tag}"
elif input_file.name.startswith("kerr_image_cuda_"):
    run_tag = f"cuda_{run_tag}"

if spectral_model and spectral_model != "monochromatic":
    run_tag = f"{run_tag}_Spectrum_{spectral_model}"

plot_dir = (
    Path("plots")
    / rho_tag
    / f"ENU{energy_folder_tag(Enu_GeV)}"
)
plot_dir.mkdir(parents=True, exist_ok=True)

data = load_kerr_image(input_file)

alpha = data[:, 2]
beta = data[:, 3]

tau = data[:, 4]
P = data[:, 5]
Iobs = data[:, 6]
captured = data[:, 7]
has_mev_channel = data.shape[1] >= 11

pixel_i = data[:, 0].astype(int)
pixel_j = data[:, 1].astype(int)
nx = int(pixel_i.max()) + 1
ny = int(pixel_j.max()) + 1

# ==========================================================
# Normalized observed intensity
# ==========================================================

Iplot = Iobs.copy()

Imax = np.max(Iplot)
if Imax > 0:
    Iplot = Iplot / Imax

# ==========================================================
# Log intensity image
# ==========================================================

plt.figure(figsize=(7, 6))

plt.tricontourf(
    alpha,
    -beta,
    np.log10(Iplot + 1e-8),
    levels=120,
    vmin=-8,
    vmax=0
)

plt.colorbar(label=r"$\log_{10}(I_{\rm obs}/I_{\rm max})$")

plt.xlabel(r"Camera coordinate $\alpha$")
plt.ylabel(r"Camera coordinate $\beta$")

plt.gca().set_aspect("equal")
plt.tight_layout()

log_plot = plot_dir / f"kerr_uhe_neutrino_image_{run_tag}.png"

plt.savefig(log_plot, dpi=250)

# ==========================================================
# Enhanced contrast image
# ==========================================================

plt.figure(figsize=(7, 6))

plt.tricontourf(
    alpha,
    -beta,
    np.sqrt(Iplot),
    levels=120
)

plt.colorbar(label=r"$\sqrt{I_{\rm obs}/I_{\rm max}}$")

plt.xlabel(r"Camera coordinate $\alpha$")
plt.ylabel(r"Camera coordinate $\beta$")

plt.gca().set_aspect("equal")
plt.tight_layout()

sqrt_plot = plot_dir / f"kerr_uhe_neutrino_image_sqrt_{run_tag}.png"

plt.savefig(sqrt_plot, dpi=250)

if has_mev_channel:
    Imev = data[:, 10]
    Imev_plot = Imev.copy()
    Imev_max = np.max(Imev_plot)
    if Imev_max > 0:
        Imev_plot = Imev_plot / Imev_max

    plt.figure(figsize=(7, 6))

    plt.tricontourf(
        alpha,
        -beta,
        np.log10(Imev_plot + 1e-8),
        levels=120,
        vmin=-8,
        vmax=0
    )

    plt.colorbar(label=r"$\log_{10}(I_{\rm MeV}/I_{\rm MeV,max})$")

    plt.xlabel(r"Camera coordinate $\alpha$")
    plt.ylabel(r"Camera coordinate $\beta$")

    plt.gca().set_aspect("equal")
    plt.tight_layout()

    mev_plot = plot_dir / f"kerr_mev_thermal_neutrino_image_{run_tag}.png"
    plt.savefig(mev_plot, dpi=250)

    uhe_grid = np.zeros((ny, nx))
    mev_grid = np.zeros((ny, nx))

    image_row = (ny - 1) - pixel_j

    uhe_grid[image_row, pixel_i] = Iplot
    mev_grid[image_row, pixel_i] = Imev_plot

    uhe_rgb = np.sqrt(np.clip(uhe_grid, 0.0, 1.0))
    mev_rgb = np.sqrt(np.clip(mev_grid, 0.0, 1.0))

    rgb = np.zeros((ny, nx, 3))
    rgb[:, :, 0] = uhe_rgb
    rgb[:, :, 1] = uhe_rgb
    rgb[:, :, 2] = mev_rgb
    rgb = np.clip(rgb, 0.0, 1.0)

    alpha_extent = [alpha.min(), alpha.max()]
    beta_extent = [-beta.max(), -beta.min()]

    fig = plt.figure(figsize=(9.2, 6))
    ax = fig.add_axes([0.10, 0.13, 0.64, 0.76])
    ax.imshow(
        rgb,
        origin="lower",
        extent=[
            alpha_extent[0],
            alpha_extent[1],
            beta_extent[0],
            beta_extent[1],
        ],
        aspect="equal"
    )

    ax.set_xlabel(r"Camera coordinate $\alpha$")
    ax.set_ylabel(r"Camera coordinate $\beta$")

    uhe_cmap = LinearSegmentedColormap.from_list(
        "uhe_yellow",
        [(0.0, 0.0, 0.0), (1.0, 1.0, 0.0)]
    )
    mev_cmap = LinearSegmentedColormap.from_list(
        "mev_blue",
        [(0.0, 0.0, 0.0), (0.0, 0.35, 1.0)]
    )
    norm = Normalize(vmin=0.0, vmax=1.0)

    uhe_sm = ScalarMappable(norm=norm, cmap=uhe_cmap)
    mev_sm = ScalarMappable(norm=norm, cmap=mev_cmap)
    uhe_sm.set_array([])
    mev_sm.set_array([])

    cax_uhe = fig.add_axes([0.79, 0.13, 0.028, 0.76])
    cax_mev = fig.add_axes([0.90, 0.13, 0.028, 0.76])

    cbar_uhe = fig.colorbar(uhe_sm, cax=cax_uhe)
    cbar_uhe.set_label(r"UHE $\sqrt{I/I_{\max}}$")

    cbar_mev = fig.colorbar(mev_sm, cax=cax_mev)
    cbar_mev.set_label(r"MeV $\sqrt{I/I_{\max}}$")

    rgb_plot = plot_dir / f"kerr_rgb_uhe_mev_composite_{run_tag}.png"
    fig.savefig(rgb_plot, dpi=250)

    side_by_side_plot = plot_dir / f"combined_side_by_side_{run_tag}.png"
    fig_sbs, axes = plt.subplots(1, 2, figsize=(12, 5.4), constrained_layout=True)

    im_uhe = axes[0].imshow(
        uhe_rgb,
        origin="lower",
        extent=[
            alpha_extent[0],
            alpha_extent[1],
            beta_extent[0],
            beta_extent[1],
        ],
        aspect="equal",
        cmap=uhe_cmap,
        vmin=0.0,
        vmax=1.0
    )
    axes[0].set_title("Neutrino UHE type")
    axes[0].set_xlabel(r"Camera coordinate $\alpha$")
    axes[0].set_ylabel(r"Camera coordinate $\beta$")
    cbar = fig_sbs.colorbar(im_uhe, ax=axes[0], fraction=0.046, pad=0.04)

    im_mev = axes[1].imshow(
        mev_rgb,
        origin="lower",
        extent=[
            alpha_extent[0],
            alpha_extent[1],
            beta_extent[0],
            beta_extent[1],
        ],
        aspect="equal",
        cmap=mev_cmap,
        vmin=0.0,
        vmax=1.0
    )
    axes[1].set_title("Neutrino MeV type")
    axes[1].set_xlabel(r"Camera coordinate $\alpha$")
    axes[1].set_ylabel(r"Camera coordinate $\beta$")
    cbar = fig_sbs.colorbar(im_mev, ax=axes[1], fraction=0.046, pad=0.04)

    fig_sbs.savefig(side_by_side_plot, dpi=250)

    rgb_overlay_plot = plot_dir / f"combined_rgb_overlay_{run_tag}.png"
    fig.savefig(rgb_overlay_plot, dpi=250)

    x_grid = np.linspace(alpha_extent[0], alpha_extent[1], nx)
    y_grid = np.linspace(beta_extent[0], beta_extent[1], ny)
    X, Y = np.meshgrid(x_grid, y_grid)

    contour_plot = plot_dir / f"combined_contours_side_by_side_{run_tag}.png"
    fig_contour, contour_axes = plt.subplots(
        1,
        2,
        figsize=(12, 5.4),
        constrained_layout=True
    )

    contour_levels = np.array([0.10, 0.30, 0.55, 0.80])

    im_uhe_contour = contour_axes[0].imshow(
        uhe_rgb,
        origin="lower",
        extent=[
            alpha_extent[0],
            alpha_extent[1],
            beta_extent[0],
            beta_extent[1],
        ],
        aspect="equal",
        cmap=uhe_cmap,
        vmin=0.0,
        vmax=1.0
    )

    if np.nanmax(uhe_rgb) > 0.0:
        uhe_levels = contour_levels[contour_levels < np.nanmax(uhe_rgb)]

        uhe_contours = contour_axes[0].contour(
            X,
            Y,
            uhe_rgb,
            levels=uhe_levels,
            colors=["#fff8a8", "#ffe05c", "#ffb000", "#ffffff"][: len(uhe_levels)],
            linewidths=[0.8, 1.0, 1.3, 1.6][: len(uhe_levels)]
        )

        contour_axes[0].clabel(
            uhe_contours,
            inline=True,
            fontsize=8,
            fmt="%.1f"
        )

    contour_axes[0].set_title("Neutrino UHE type")
    contour_axes[0].set_xlabel(r"Camera coordinate $\alpha$")
    contour_axes[0].set_ylabel(r"Camera coordinate $\beta$")
    cbar = fig_contour.colorbar(
        im_uhe_contour,
        ax=contour_axes[0],
        fraction=0.046,
        pad=0.04
    )
    cbar.set_label(r"$\sqrt{I/I_{\max}}$")

    im_mev_contour = contour_axes[1].imshow(
        mev_rgb,
        origin="lower",
        extent=[
            alpha_extent[0],
            alpha_extent[1],
            beta_extent[0],
            beta_extent[1],
        ],
        aspect="equal",
        cmap=mev_cmap,
        vmin=0.0,
        vmax=1.0
    )

    if np.nanmax(mev_rgb) > 0.0:
        mev_levels = contour_levels[contour_levels < np.nanmax(mev_rgb)]

        mev_contours = contour_axes[1].contour(
            X,
            Y,
            mev_rgb,
            levels=mev_levels,
            colors=["#c9e7ff", "#87c6ff", "#4d96ff", "#ffffff"][: len(mev_levels)],
            linewidths=[0.7, 0.9, 1.1, 1.35][: len(mev_levels)]
        )

        contour_axes[1].clabel(
            mev_contours,
            inline=True,
            fontsize=8,
            fmt="%.1f"
        )

    contour_axes[1].set_title("Neutrino MeV type")
    contour_axes[1].set_xlabel(r"Camera coordinate $\alpha$")
    contour_axes[1].set_ylabel(r"Camera coordinate $\beta$")
    cbar = fig_contour.colorbar(
        im_mev_contour,
        ax=contour_axes[1],
        fraction=0.046,
        pad=0.04
    )
    cbar.set_label(r"$\sqrt{I/I_{\max}}$")

    fig_contour.savefig(contour_plot, dpi=250)

plt.close("all")

print(f"Saved: {log_plot}")
print(f"Saved: {sqrt_plot}")
if has_mev_channel:
    print(f"Saved: {mev_plot}")
    print(f"Saved: {rgb_plot}")
    print(f"Saved: {side_by_side_plot}")
    print(f"Saved: {rgb_overlay_plot}")
    print(f"Saved: {contour_plot}")
