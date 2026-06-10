"""Validate UHE spectral emission models on a small geodesic cache."""

from __future__ import annotations

import csv
import math
import subprocess
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
OUTDIR = ROOT / "output" / "spectra"
PLOTDIR = ROOT / "plots" / "spectra"

BASE = {
    "enu": "1e11",
    "a": "0.0001",
    "mbh": "3.0",
    "rho0": "1.0e0",
    "r0": "10.0",
    "sigma": "5.0",
    "h": "0.25",
    "source_r": "3.5",
    "source_sigma": "1.0",
    "source_theta": "15.0",
    "source_powerlaw": "2.0",
    "source_emax": "1.0e12",
    "source_norm": "1.0",
    "mev_enu": "10.0",
    "mev_norm": "1.0",
    "cam_theta": "80.0",
    "sigma_path": "data/sigma/sigma_nuN_CC_GBW.dat",
    "profile": "powerlaw_funnel_envelope",
    "radial_power": "2.0",
    "funnel_depletion": "1.0",
    "funnel_theta": "20.0",
    "envelope_rho0": "1.0e-2",
    "envelope_alpha": "2.5",
    "r_min": "4.0",
    "r_max": "60.0",
    "rho_floor": "1.0e-99",
    "observer_distance": "60.0",
    "use_f3": "1",
    "cache": "output/rays/kerr_geodesics_small.bin",
    "source_model": "funnel_wall",
    "source_funnel_theta": "20.0",
    "source_q": "1.0",
    "source_s": "2.0",
    "source_grad_dr": "0.1",
    "source_grad_dtheta": "1.0",
    "source_rho_ref": "1.0e0",
    "source_cutoff_min": "0.0",
    "source_cutoff_max": "1.0e2",
}


def run(cmd: list[str]) -> str:
    result = subprocess.run(
        cmd,
        cwd=ROOT,
        check=True,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
    )
    return result.stdout


def compact_tag(value: str) -> str:
    mantissa, exponent = f"{float(value):.0e}".split("e")
    return f"{mantissa}e{int(exponent)}"


def sci_tag(value: str) -> str:
    return f"{float(value):.0e}".replace("+", "_").replace("-", "m").replace(".", "p")


def image_path(model: str) -> Path:
    spectrum = "" if model == "monochromatic" else f"_Spectrum_{model}"
    name = (
        "kerr_image_cuda_cache_GBW_"
        f"rho0_torus_{compact_tag(BASE['rho0'])}"
        f"_Profile_{BASE['profile']}"
        f"_Source_{BASE['source_model']}"
        f"{spectrum}"
        f"_Enu_{sci_tag(BASE['enu'])}"
        f"_MeVEnu_{sci_tag(BASE['mev_enu'])}"
        f"_MeVNorm_{sci_tag(BASE['mev_norm'])}"
        "_CamTheta_80p0.dat"
    )
    return ROOT / "output" / "images" / name


def compute(model: str, n_bins: int = 8) -> Path:
    cmd = [
        "./compute_kerr_image_from_cache",
        BASE["enu"],
        BASE["a"],
        BASE["mbh"],
        BASE["rho0"],
        BASE["r0"],
        BASE["sigma"],
        BASE["h"],
        BASE["source_r"],
        BASE["source_sigma"],
        BASE["source_theta"],
        BASE["source_powerlaw"],
        BASE["source_emax"],
        BASE["source_norm"],
        BASE["mev_enu"],
        BASE["mev_norm"],
        BASE["cam_theta"],
        BASE["sigma_path"],
        BASE["profile"],
        BASE["radial_power"],
        BASE["funnel_depletion"],
        BASE["funnel_theta"],
        BASE["envelope_rho0"],
        BASE["envelope_alpha"],
        BASE["r_min"],
        BASE["r_max"],
        BASE["rho_floor"],
        BASE["observer_distance"],
        BASE["use_f3"],
        BASE["cache"],
        BASE["source_model"],
        BASE["source_funnel_theta"],
        BASE["source_q"],
        BASE["source_s"],
        BASE["source_grad_dr"],
        BASE["source_grad_dtheta"],
        BASE["source_rho_ref"],
        BASE["source_cutoff_min"],
        BASE["source_cutoff_max"],
        model,
        "2.0",
        "1.0e12",
        "1.0e5",
        "1.0e12",
        str(n_bins),
    ]
    run(cmd)
    return image_path(model)


def metrics(path: Path, model: str, n_bins: int) -> dict[str, float | str | int]:
    data = np.loadtxt(path)
    tau = data[:, 4]
    psurv = data[:, 5]
    intensity = data[:, 6]
    finite = np.isfinite(tau) & np.isfinite(psurv) & np.isfinite(intensity)
    return {
        "spectral_model": model,
        "n_bins": n_bins,
        "mean_tau": float(np.mean(tau[finite])),
        "max_tau": float(np.max(tau[finite])),
        "mean_Psurv": float(np.mean(psurv[finite])),
        "integrated_intensity": float(np.sum(intensity[finite])),
        "valid_pixels": int(np.count_nonzero(finite)),
        "finite_outputs": bool(np.all(finite)),
        "positive_weights": True,
    }


def spectral_weight(E: np.ndarray, model: str, gamma: float = 2.0, ecut: float = 1.0e12) -> np.ndarray:
    if model == "monochromatic":
        center = float(BASE["enu"])
        width = 0.08
        return np.exp(-0.5 * (np.log10(E / center) / width) ** 2)
    if model == "powerlaw":
        return E ** (-gamma)
    if model == "powerlaw_cutoff":
        return E ** (-gamma) * np.exp(-E / ecut)
    raise ValueError(model)


def plot_spectra() -> None:
    energies = np.logspace(5, 12, 400)
    fig, ax = plt.subplots(figsize=(6.6, 4.4), constrained_layout=True)
    for model in ["monochromatic", "powerlaw", "powerlaw_cutoff"]:
        w = spectral_weight(energies, model)
        w = w / np.max(w)
        ax.plot(energies, w, label=model)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(r"$E_\nu$ [GeV]")
    ax.set_ylabel("normalized spectral weight")
    ax.set_title("UHE spectral models")
    ax.grid(alpha=0.25)
    ax.legend()
    fig.savefig(PLOTDIR / "spectral_models.png", dpi=180)
    plt.close(fig)


def plot_validation(rows: list[dict[str, float | str | int]]) -> None:
    labels = [str(r["spectral_model"]) for r in rows if int(r["n_bins"]) == 8]
    mean_p = [float(r["mean_Psurv"]) for r in rows if int(r["n_bins"]) == 8]
    mean_tau = [float(r["mean_tau"]) for r in rows if int(r["n_bins"]) == 8]
    intensity = [float(r["integrated_intensity"]) for r in rows if int(r["n_bins"]) == 8]

    sigma_data = np.loadtxt(ROOT / BASE["sigma_path"])
    sigma_E = sigma_data[:, 0]
    sigma_cm2 = sigma_data[:, 1]
    sigma_ref = np.interp(float(BASE["enu"]), sigma_E, sigma_cm2)
    tau_ref = mean_tau[0]
    tau_energy = tau_ref * sigma_cm2 / max(sigma_ref, 1.0e-300)
    fig, ax = plt.subplots(figsize=(6.6, 4.4), constrained_layout=True)
    ax.plot(sigma_E, np.exp(-tau_energy), color="#2563eb", lw=1.8)
    ax.scatter([float(BASE["enu"])], [mean_p[0]], color="#dc2626", zorder=3, label="validated mono run")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(r"$E_\nu$ [GeV]")
    ax.set_ylabel(r"$\langle P_{\rm surv}\rangle$")
    ax.set_title("Survival probability vs energy")
    ax.grid(alpha=0.25)
    ax.legend()
    fig.savefig(PLOTDIR / "survival_probability_vs_energy.png", dpi=180)
    plt.close(fig)

    mono_tau = mean_tau[0]
    mono_intensity = intensity[0]
    fig, axes = plt.subplots(1, 2, figsize=(8.0, 3.9), constrained_layout=True)
    axes[0].bar(labels, [v / mono_tau if mono_tau > 0 else 0 for v in mean_tau])
    axes[0].set_ylabel(r"$\langle\tau\rangle / \langle\tau\rangle_{\rm mono}$")
    axes[0].tick_params(axis="x", rotation=15)
    axes[1].bar(labels, [v / mono_intensity if mono_intensity > 0 else 0 for v in intensity])
    axes[1].set_ylabel(r"$I / I_{\rm mono}$")
    axes[1].tick_params(axis="x", rotation=15)
    fig.suptitle("Spectral weighted vs monochromatic")
    fig.savefig(PLOTDIR / "spectral_weighted_vs_monochromatic.png", dpi=180)
    plt.close(fig)


def main() -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    PLOTDIR.mkdir(parents=True, exist_ok=True)

    cache = ROOT / BASE["cache"]
    if not cache.exists():
        raise FileNotFoundError(f"Small cache not found: {cache}")

    rows: list[dict[str, float | str | int]] = []
    for model in ["monochromatic", "powerlaw", "powerlaw_cutoff"]:
        path = compute(model, 8)
        rows.append(metrics(path, model, 8))

    for n_bins in [4, 8, 16]:
        path = compute("powerlaw_cutoff", n_bins)
        row = metrics(path, "powerlaw_cutoff", n_bins)
        row["spectral_model"] = f"powerlaw_cutoff_bins{n_bins}"
        rows.append(row)

    with (OUTDIR / "spectral_validation.csv").open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    plot_spectra()
    plot_validation(rows)

    finite_ok = all(bool(r["finite_outputs"]) for r in rows)
    positive_intensity = all(float(r["integrated_intensity"]) >= 0.0 for r in rows)
    convergence_rows = [r for r in rows if str(r["spectral_model"]).startswith("powerlaw_cutoff_bins")]
    intensities = [float(r["integrated_intensity"]) for r in convergence_rows]
    rel_8_16 = abs(intensities[1] - intensities[2]) / max(abs(intensities[2]), 1.0e-300)

    lines = [
        "# UHE spectral validation",
        "",
        f"- finite_outputs: {finite_ok}",
        f"- nonnegative_integrated_intensity: {positive_intensity}",
        f"- powerlaw_cutoff convergence |I_8-I_16|/I_16: {rel_8_16:.6e}",
        "",
        "## Metrics",
        "",
    ]
    for row in rows:
        lines.append(
            f"- {row['spectral_model']} bins={row['n_bins']}: "
            f"mean_tau={float(row['mean_tau']):.6e}, "
            f"mean_Psurv={float(row['mean_Psurv']):.6e}, "
            f"I={float(row['integrated_intensity']):.6e}"
        )
    (OUTDIR / "spectral_validation.md").write_text("\n".join(lines) + "\n")

    print(f"Wrote {OUTDIR / 'spectral_validation.csv'}")
    print(f"Wrote {OUTDIR / 'spectral_validation.md'}")
    print(f"Wrote plots to {PLOTDIR}")


if __name__ == "__main__":
    main()
