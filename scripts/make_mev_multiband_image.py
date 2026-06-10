#!/usr/bin/env python3
"""Create a false-color MeV energy-band diagnostic image."""

from __future__ import annotations

import os
import subprocess
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
IMG_DIR = ROOT / "output" / "images"
OUT_DIR = ROOT / "output" / "mev_neutrinosphere"
PLOT_DIR = ROOT / "plots" / "mev_physics"
CACHE = os.environ.get("SMALL_CACHE_PATH", "output/rays/kerr_geodesics_small.bin")

BANDS = [
    ("low-energy MeV band", "blue", 3.0, 8.0),
    ("intermediate-energy MeV band", "green", 8.0, 20.0),
    ("high-energy MeV band", "red", 20.0, 50.0),
]


def env(name: str, default: str) -> str:
    return os.environ.get(name, default)


def sci_tag(value: float) -> str:
    return f"{value:.0e}".replace("+", "_").replace("-", "m").replace(".", "p")


def compact_tag(value: float) -> str:
    mantissa, exponent = f"{value:.0e}".split("e")
    return f"{mantissa}e{int(exponent)}"


def expected_path(rho0: float, enu: float, mev_enu: float, mev_norm: float, cam_theta: float, e_min: float, e_max: float) -> Path:
    profile = env("DENSITY_PROFILE", "gaussian")
    profile_part = "" if profile == "gaussian" else f"_Profile_{profile}"
    return IMG_DIR / (
        "kerr_image_cuda_cache_GBW"
        f"_rho0_torus_{compact_tag(rho0)}"
        f"{profile_part}"
        f"_Enu_{sci_tag(enu)}"
        f"_MeVEnu_{sci_tag(mev_enu)}"
        f"_MeVNorm_{sci_tag(mev_norm)}"
        f"_MeVSpectrum_fermi_dirac_band_E{sci_tag(e_min)}_{sci_tag(e_max)}"
        f"_CamTheta_{cam_theta:.1f}".replace(".", "p")
        + ".dat"
    )


def compute_band(e_min: float, e_max: float) -> Path:
    enu = float(env("SMALL_CACHE_ENU", "1e11"))
    rho0 = float(env("SMALL_CACHE_TORUS_RHO0", env("MEV_VALIDATION_TORUS_RHO0", "1.0e10")))
    mev_enu = float(env("MEV_ENU", env("MEV_ENERGY_MEV", "10.0")))
    mev_norm = float(env("MEV_NORM", "1.0"))
    cam_theta = float(env("CAM_THETA_DEG", "80.0"))
    cmd = [
        "./compute_kerr_image_from_cache",
        str(enu),
        env("ASPIN", "0.0001"),
        env("MBH_MSUN", "3.0"),
        str(rho0),
        env("TORUS_R0_RG", "10.0"),
        env("TORUS_SIGMA_RG", "5.0"),
        env("TORUS_H_OVER_R", "0.25"),
        env("SOURCE_R_RG", "3.5"),
        env("SOURCE_SIGMA_RG", "1.0"),
        env("SOURCE_THETA_DEG", "15.0"),
        env("SOURCE_POWERLAW", "2.0"),
        env("SOURCE_EMAX_GEV", "1.0e12"),
        env("SOURCE_NORM", "1.0"),
        str(mev_enu),
        str(mev_norm),
        str(cam_theta),
        "data/sigma/sigma_nuN_CC_GBW.dat",
        env("DENSITY_PROFILE", "gaussian"),
        env("TORUS_RADIAL_POWER", "2.0"),
        env("FUNNEL_DEPLETION", "0.0"),
        env("FUNNEL_THETA_DEG", "15.0"),
        env("ENVELOPE_RHO0", "0.0"),
        env("ENVELOPE_ALPHA", "2.5"),
        env("TORUS_R_MIN_RG", "4.0"),
        env("TORUS_R_MAX_RG", "60.0"),
        env("RHO_FLOOR", env("MEV_VALIDATION_RHO_FLOOR", "1.0e-20")),
        env("CAM_R_OBS_RG", "60.0"),
        env("USE_F3", "1"),
        CACHE,
        env("SOURCE_MODEL", "inner_ring"),
        env("SOURCE_FUNNEL_THETA_DEG", "20.0"),
        env("SOURCE_DENSITY_Q", "1.0"),
        env("SOURCE_RADIAL_S", "2.0"),
        env("SOURCE_GRADIENT_DR_RG", "0.1"),
        env("SOURCE_GRADIENT_DTHETA_DEG", "1.0"),
        env("SOURCE_RHO_REF", str(rho0)),
        env("SOURCE_CUTOFF_MIN", "0.0"),
        env("SOURCE_CUTOFF_MAX", "1.0e2"),
        env("SPECTRAL_MODEL", "monochromatic"),
        env("SPECTRAL_GAMMA", "2.0"),
        env("SPECTRAL_ECUT_GEV", "1.0e12"),
        env("SPECTRAL_E_MIN_GEV", "1.0e5"),
        env("SPECTRAL_E_MAX_GEV", "1.0e12"),
        env("SPECTRAL_N_BINS", "8"),
        env("MEV_MODEL", "physical"),
        env("MEV_FLAVOR", "anti_nu_e"),
        env("MEV_INCLUDE_URCA", "1"),
        env("MEV_INCLUDE_PAIR", "1"),
        env("MEV_INCLUDE_BREMS", "1"),
        env("MEV_INCLUDE_ABSORPTION", "1"),
        env("MEV_INCLUDE_SCATTERING", "1"),
        env("MEV_THERMAL_PROFILE", "inner_hot_torus"),
        env("MEV_YE_PROFILE", "neutron_rich_torus"),
        env("MEV_T0_MEV", "6.0"),
        env("MEV_T_FLOOR_MEV", "0.1"),
        env("MEV_T_POWER", "0.2"),
        env("MEV_YE_TORUS", "0.25"),
        env("MEV_YE_FUNNEL", "0.55"),
        env("MEV_YE_ENVELOPE", "0.45"),
        env("MEV_YE_FLOOR", "0.01"),
        env("MEV_YE_CEIL", "0.60"),
        "fermi_dirac_band",
        str(e_min),
        str(e_max),
        env("MEV_N_BINS", "8"),
        env("MEV_USE_DEGENERACY_CORRECTION", "0"),
        env("MEV_INCLUDE_ABS_N", "1"),
        env("MEV_INCLUDE_ABS_P", "1"),
        env("MEV_INCLUDE_SCAT_N", "1"),
        env("MEV_INCLUDE_SCAT_P", "1"),
        env("MEV_INCLUDE_SCAT_E", "1"),
    ]
    subprocess.run(cmd, cwd=ROOT, check=True)
    path = expected_path(rho0, enu, mev_enu, mev_norm, cam_theta, e_min, e_max)
    if not path.exists():
        raise FileNotFoundError(path)
    return path


def load_grid(path: Path) -> tuple[np.ndarray, np.ndarray]:
    data = np.loadtxt(path, comments="#")
    if data.ndim == 1:
        data = data[None, :]
    nx = int(np.max(data[:, 0])) + 1
    ny = int(np.max(data[:, 1])) + 1
    grid = np.zeros((ny, nx))
    for row in data:
        i = int(row[0])
        j = int(row[1])
        grid[ny - 1 - j, i] = row[10]
    return grid, data


def norm_channel(grid: np.ndarray) -> np.ndarray:
    g = np.clip(grid, 0.0, None)
    maxv = float(np.nanmax(g))
    if maxv <= 0.0:
        return np.zeros_like(g)
    return np.sqrt(g / maxv)


def main() -> None:
    PLOT_DIR.mkdir(parents=True, exist_ok=True)
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    channels = []
    rows = [
        "# product false-color MeV energy-band diagnostic image",
        "# blue low-energy MeV band 3-8 MeV",
        "# green intermediate-energy MeV band 8-20 MeV",
        "# red high-energy MeV band 20-50 MeV",
        "# note colors encode MeV neutrino energy intervals, not photon colors",
        "band,false_color,E_min_MeV,E_max_MeV,integrated_I_mev,data_path",
    ]
    for label, color, e_min, e_max in BANDS:
        path = compute_band(e_min, e_max)
        grid, data = load_grid(path)
        channels.append(norm_channel(grid))
        rows.append(f"{label},{color},{e_min:g},{e_max:g},{float(np.sum(data[:,10])):.8e},{path.relative_to(ROOT)}")

    blue, green, red = channels
    rgb = np.dstack([red, green, blue])
    plt.figure(figsize=(5.2, 5.0))
    plt.imshow(rgb, origin="lower")
    plt.title("False-color MeV energy-band diagnostic image\nB=3-8 MeV, G=8-20 MeV, R=20-50 MeV")
    plt.xlabel("camera pixel")
    plt.ylabel("camera pixel")
    plt.tight_layout()
    out = PLOT_DIR / "mev_multiband_false_color_image.png"
    plt.savefig(out, dpi=240)
    plt.close()
    csv_path = OUT_DIR / "mev_multiband_flux.csv"
    csv_path.write_text("\n".join(rows) + "\n")
    print(f"Saved: {out}")
    print(f"Saved: {csv_path}")


if __name__ == "__main__":
    main()
