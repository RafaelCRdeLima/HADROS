"""Plot emitted vs observed UHE spectrum using an existing geodesic cache."""

from __future__ import annotations

import os
import subprocess
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
OUTDIR = ROOT / "output" / "spectra"
PLOTDIR = ROOT / "plots" / "spectra"


def env(name: str, default: str) -> str:
    return os.environ.get(name, default)


def sci_tag(value: float) -> str:
    return f"{value:.0e}".replace("+", "_").replace("-", "m").replace(".", "p")


def compact_tag(value: float) -> str:
    mantissa, exponent = f"{value:.0e}".split("e")
    return f"{mantissa}e{int(exponent)}"


def spectral_weight(E: np.ndarray, model: str, gamma: float, ecut: float) -> np.ndarray:
    if model == "monochromatic":
        center = float(env("ENU", env("SMALL_CACHE_ENU", "1e11")))
        width = 0.08
        return np.exp(-0.5 * (np.log10(E / center) / width) ** 2)
    if model == "powerlaw":
        return E ** (-gamma)
    if model == "powerlaw_cutoff":
        return E ** (-gamma) * np.exp(-E / ecut)
    raise ValueError(f"Unknown spectral model: {model}")


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


def run_image(E: float, source_emax: str) -> Path:
    rho0 = env("SMALL_CACHE_TORUS_RHO0", env("TORUS_RHO0", "1.0e0"))
    profile = env("DENSITY_PROFILE", "powerlaw_funnel_envelope")
    source_model = env("SOURCE_MODEL", "funnel_wall")
    gamma = env("SPECTRAL_GAMMA", "2.0")

    cmd = [
        "./compute_kerr_image_from_cache",
        f"{E:.12g}",
        env("ASPIN", "0.0001"),
        env("MBH_MSUN", "3.0"),
        rho0,
        env("TORUS_R0_RG", "10.0"),
        env("TORUS_SIGMA_RG", "5.0"),
        env("TORUS_H_OVER_R", "0.25"),
        env("SOURCE_R_RG", "3.5"),
        env("SOURCE_SIGMA_RG", "1.0"),
        env("SOURCE_THETA_DEG", "15.0"),
        gamma,
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
        env("SOURCE_RHO_REF", rho0),
        env("SOURCE_CUTOFF_MIN", "0.0"),
        env("SOURCE_CUTOFF_MAX", "1.0e2"),
        "monochromatic",
        gamma,
        source_emax,
        f"{E:.12g}",
        f"{E:.12g}",
        "1",
    ]
    subprocess.run(cmd, cwd=ROOT, check=True)
    return image_path(E, float(rho0), profile, source_model)


def main() -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    PLOTDIR.mkdir(parents=True, exist_ok=True)

    model = env("SPECTRAL_MODEL", "powerlaw_cutoff")
    gamma = float(env("SPECTRAL_GAMMA", "2.0"))
    ecut = float(env("SPECTRAL_ECUT_GEV", "1.0e12"))
    e_min = float(env("SPECTRAL_E_MIN_GEV", "1.0e5"))
    e_max = float(env("SPECTRAL_E_MAX_GEV", "1.0e12"))
    n_bins = max(int(env("SPECTRAL_N_BINS", "8")), 2)

    if model == "powerlaw":
        source_emax = "1.0e99"
    elif model == "powerlaw_cutoff":
        source_emax = f"{ecut:.12g}"
    else:
        source_emax = env("SOURCE_EMAX_GEV", "1.0e12")

    energies = np.logspace(np.log10(e_min), np.log10(e_max), n_bins)
    emitted = spectral_weight(energies, model, gamma, ecut)
    observed = []
    mean_tau = []
    mean_psurv = []

    for E in energies:
        path = run_image(float(E), source_emax)
        data = np.loadtxt(path)
        observed.append(float(np.sum(data[:, 6])))
        mean_tau.append(float(np.mean(data[:, 4])))
        mean_psurv.append(float(np.mean(data[:, 5])))

    observed = np.array(observed)
    emitted = np.array(emitted)
    mean_tau = np.array(mean_tau)
    mean_psurv = np.array(mean_psurv)

    emitted_norm = emitted / max(np.max(emitted), 1.0e-300)
    observed_norm = observed / max(np.max(observed), 1.0e-300)

    table = np.column_stack(
        [energies, emitted, observed, emitted_norm, observed_norm, mean_tau, mean_psurv]
    )
    header = (
        "E_GeV emitted_weight observed_intensity "
        "emitted_norm observed_norm mean_tau mean_Psurv"
    )
    np.savetxt(OUTDIR / "emitted_vs_observed_spectrum.csv", table, header=header)

    fig, axes = plt.subplots(1, 2, figsize=(9.0, 4.2), constrained_layout=True)
    axes[0].plot(energies, emitted_norm, marker="o", label="emitted")
    axes[0].plot(energies, observed_norm, marker="s", label="observed")
    axes[0].set_xscale("log")
    axes[0].set_yscale("log")
    axes[0].set_xlabel(r"$E_\nu$ [GeV]")
    axes[0].set_ylabel("normalized spectrum")
    axes[0].set_title(f"{model}: emitted vs observed")
    axes[0].grid(alpha=0.25)
    axes[0].legend()

    ratio = observed_norm / np.maximum(emitted_norm, 1.0e-300)
    axes[1].plot(energies, ratio, marker="o", color="#dc2626")
    axes[1].set_xscale("log")
    axes[1].set_yscale("log")
    axes[1].set_xlabel(r"$E_\nu$ [GeV]")
    axes[1].set_ylabel("observed / emitted, normalized")
    axes[1].set_title("spectral attenuation")
    axes[1].grid(alpha=0.25)

    output = PLOTDIR / "emitted_vs_observed_spectrum.png"
    fig.savefig(output, dpi=190)
    plt.close(fig)

    print(f"Saved: {OUTDIR / 'emitted_vs_observed_spectrum.csv'}")
    print(f"Saved: {output}")


if __name__ == "__main__":
    main()
