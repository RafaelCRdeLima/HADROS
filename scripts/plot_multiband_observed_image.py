"""Create a multiband observed UHE image from monochromatic image slices."""

from __future__ import annotations

import os
import subprocess
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LinearSegmentedColormap


ROOT = Path(__file__).resolve().parents[1]
PLOTDIR = ROOT / "plots" / "spectra"
OUTDIR = ROOT / "output" / "spectra"


def env(name: str, default: str) -> str:
    return os.environ.get(name, default)


def sci_tag(value: float) -> str:
    return f"{value:.0e}".replace("+", "_").replace("-", "m").replace(".", "p")


def compact_tag(value: float) -> str:
    mantissa, exponent = f"{value:.0e}".split("e")
    return f"{mantissa}e{int(exponent)}"


def spectral_weight(E: np.ndarray, model: str, gamma: float, ecut: float) -> np.ndarray:
    if model == "powerlaw":
        return E ** (-gamma)
    if model == "powerlaw_cutoff":
        return E ** (-gamma) * np.exp(-E / ecut)
    if model == "monochromatic":
        center = float(env("ENU", env("SMALL_CACHE_ENU", "1e11")))
        width = 0.08
        return np.exp(-0.5 * (np.log10(E / center) / width) ** 2)
    raise ValueError(model)


def image_path(E: float, rho0: float, profile: str, source_model: str) -> Path:
    tag = ""
    if profile != "gaussian":
        tag += f"_Profile_{profile}"
    if source_model != "inner_ring":
        tag += f"_Source_{source_model}"
    name = (
        "kerr_image_cuda_cache_GBW_"
        f"rho0_torus_{compact_tag(rho0)}"
        f"{tag}"
        f"_Enu_{sci_tag(E)}"
        "_MeVEnu_1e_01"
        "_MeVNorm_1e_00"
        "_CamTheta_80p0.dat"
    )
    return ROOT / "output" / "images" / name


def ensure_image(E: float, source_emax: str) -> Path:
    rho0 = float(env("SMALL_CACHE_TORUS_RHO0", env("TORUS_RHO0", "1.0e0")))
    profile = env("DENSITY_PROFILE", "powerlaw_funnel_envelope")
    source_model = env("SOURCE_MODEL", "funnel_wall")
    path = image_path(E, rho0, profile, source_model)
    if path.exists():
        return path

    cmd = [
        "./compute_kerr_image_from_cache",
        f"{E:.12g}",
        env("ASPIN", "0.0001"),
        env("MBH_MSUN", "3.0"),
        f"{rho0:.12g}",
        env("TORUS_R0_RG", "10.0"),
        env("TORUS_SIGMA_RG", "5.0"),
        env("TORUS_H_OVER_R", "0.25"),
        env("SOURCE_R_RG", "3.5"),
        env("SOURCE_SIGMA_RG", "1.0"),
        env("SOURCE_THETA_DEG", "15.0"),
        env("SPECTRAL_GAMMA", "2.0"),
        source_emax,
        env("SOURCE_NORM", "1.0"),
        env("MEV_ENU", "10.0"),
        env("MEV_NORM", "1.0"),
        env("CAM_THETA_DEG", "80.0"),
        "data/sigma/sigma_nuN_CC_GBW.dat",
        profile,
        env("TORUS_RADIAL_POWER", "2.0"),
        env("FUNNEL_DEPLETION", "1.0"),
        env("FUNNEL_THETA_DEG", "20.0"),
        env("ENVELOPE_RHO0", "1.0e-2"),
        env("ENVELOPE_ALPHA", "2.5"),
        env("TORUS_R_MIN_RG", "4.0"),
        env("TORUS_R_MAX_RG", "60.0"),
        env("RHO_FLOOR", "1.0e-99"),
        env("CAM_R_OBS_RG", "60.0"),
        env("USE_F3", "1"),
        env("SMALL_CACHE_PATH", "output/rays/kerr_geodesics_small.bin"),
        source_model,
        env("SOURCE_FUNNEL_THETA_DEG", "20.0"),
        env("SOURCE_DENSITY_Q", "1.0"),
        env("SOURCE_RADIAL_S", "2.0"),
        env("SOURCE_GRADIENT_DR_RG", "0.1"),
        env("SOURCE_GRADIENT_DTHETA_DEG", "1.0"),
        env("SOURCE_RHO_REF", f"{rho0:.12g}"),
        env("SOURCE_CUTOFF_MIN", "0.0"),
        env("SOURCE_CUTOFF_MAX", "1.0e2"),
        "monochromatic",
        env("SPECTRAL_GAMMA", "2.0"),
        source_emax,
        f"{E:.12g}",
        f"{E:.12g}",
        "1",
    ]
    subprocess.run(cmd, cwd=ROOT, check=True)
    return path


def image_grid(data: np.ndarray, column: int) -> tuple[np.ndarray, tuple[float, float, float, float]]:
    pixel_i = data[:, 0].astype(int)
    pixel_j = data[:, 1].astype(int)
    nx = int(pixel_i.max()) + 1
    ny = int(pixel_j.max()) + 1
    row = (ny - 1) - pixel_j
    grid = np.zeros((ny, nx))
    grid[row, pixel_i] = data[:, column]
    alpha = data[:, 2]
    beta = data[:, 3]
    extent = (float(alpha.min()), float(alpha.max()), float(-beta.max()), float(-beta.min()))
    return grid, extent


def normalize(x: np.ndarray) -> np.ndarray:
    vmax = np.nanpercentile(x, 99.5)
    if vmax <= 0.0:
        vmax = np.nanmax(x)
    if vmax <= 0.0:
        return x
    return np.sqrt(np.clip(x / vmax, 0.0, 1.0))


def main() -> None:
    PLOTDIR.mkdir(parents=True, exist_ok=True)
    OUTDIR.mkdir(parents=True, exist_ok=True)

    model = env("SPECTRAL_MODEL", "powerlaw_cutoff")
    gamma = float(env("SPECTRAL_GAMMA", "2.0"))
    ecut = float(env("SPECTRAL_ECUT_GEV", "1.0e12"))
    e_min = float(env("SPECTRAL_E_MIN_GEV", "1.0e5"))
    e_max = float(env("SPECTRAL_E_MAX_GEV", "1.0e12"))
    n_bins = max(int(env("SPECTRAL_N_BINS", "8")), 3)
    source_emax = f"{ecut:.12g}" if model == "powerlaw_cutoff" else "1.0e99"

    edges = np.logspace(np.log10(e_min), np.log10(e_max), n_bins + 1)
    energies = np.sqrt(edges[:-1] * edges[1:])
    dE = edges[1:] - edges[:-1]
    weights = spectral_weight(energies, model, gamma, ecut) * dE

    bands = [
        ("low-energy UHE band", "blue", e_min, 1.0e7, np.array([0.15, 0.45, 1.00])),
        ("intermediate-energy UHE band", "green", 1.0e7, 1.0e10, np.array([0.05, 0.95, 0.55])),
        ("high-energy UHE band", "red", 1.0e10, e_max * 1.0001, np.array([1.00, 0.18, 0.12])),
    ]

    band_grids: dict[str, np.ndarray] = {}
    band_flux: dict[str, float] = {}
    extent = None

    for E, weight in zip(energies, weights):
        path = ensure_image(float(E), source_emax)
        data = np.loadtxt(path)
        grid, extent_i = image_grid(data, 6)
        extent = extent_i
        for name, _, lo, hi, _ in bands:
            if lo <= E < hi:
                band_grids[name] = band_grids.get(name, np.zeros_like(grid)) + weight * grid
                band_flux[name] = band_flux.get(name, 0.0) + float(weight * np.sum(grid))

    if extent is None:
        raise RuntimeError("No image slices were loaded.")

    rgb = np.zeros((*next(iter(band_grids.values())).shape, 3))
    for name, _, _, _, color in bands:
        grid = band_grids.get(name)
        if grid is None:
            continue
        rgb += normalize(grid)[:, :, None] * color[None, None, :]
    rgb = np.clip(rgb, 0.0, 1.0)

    fig, ax = plt.subplots(figsize=(7.0, 6.2), constrained_layout=True)
    ax.imshow(rgb, origin="lower", extent=extent, aspect="equal")
    ax.set_xlabel(r"Camera coordinate $\alpha$")
    ax.set_ylabel(r"Camera coordinate $\beta$")
    title = (
        "UHE energy-band composite image\n"
        f"blue: {e_min:.0e}-{1.0e7:.0e} GeV, "
        f"green: {1.0e7:.0e}-{1.0e10:.0e} GeV, "
        f"red: {1.0e10:.0e}-{e_max:.0e} GeV"
    )
    ax.set_title(title, fontsize=10.5)

    handles = []
    labels = []
    for name, false_color, lo, hi, color in bands:
        cmap = LinearSegmentedColormap.from_list(name, [(0, 0, 0), tuple(color)])
        handles.append(plt.Line2D([0], [0], color=cmap(1.0), lw=4))
        labels.append(f"{false_color}: {name}\n{lo:.0e}-{hi:.0e} GeV")
    ax.legend(handles, labels, loc="upper right", framealpha=0.8, fontsize=8)
    ax.text(
        0.5,
        -0.13,
        "False colors encode UHE neutrino energy intervals, not physical photon colors.",
        transform=ax.transAxes,
        ha="center",
        va="top",
        fontsize=8.2,
        color="#333333",
    )

    output = PLOTDIR / "observed_multiband_image.png"
    fig.savefig(output, dpi=220)
    plt.close(fig)

    rows = [
        "# product energy-band composite image",
        "# false_color_note colors encode UHE neutrino energy intervals, not physical photon colors",
        "# blue low-energy UHE band",
        "# green intermediate-energy UHE band",
        "# red high-energy UHE band",
        "band,false_color,E_min_GeV,E_max_GeV,weighted_observed_flux",
    ]
    for name, false_color, lo, hi, _ in bands:
        rows.append(f"{name},{false_color},{lo:.8e},{hi:.8e},{band_flux.get(name, 0.0):.8e}")
    (OUTDIR / "observed_multiband_flux.csv").write_text("\n".join(rows) + "\n")

    print(f"Saved: {output}")
    print(f"Saved: {OUTDIR / 'observed_multiband_flux.csv'}")


if __name__ == "__main__":
    main()
