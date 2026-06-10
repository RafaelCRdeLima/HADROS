#!/usr/bin/env python3
"""Validate UHE observables against the DIS cross-section model choice.

This script deliberately keeps the astrophysical setup fixed inside each
regime and varies only the sigma_nuN table passed to the radiative-transfer
executable.
"""

from __future__ import annotations

import csv
import math
import os
import re
import subprocess
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = ROOT / "data" / "sigma"
OUTPUT_DIR = ROOT / "output" / "dis_validation"
PLOT_DIR = ROOT / "plots" / "dis_validation"
IMAGE_DIR = ROOT / "output" / "images"
RAY_DIR = ROOT / "output" / "rays"

GEV_MINUS2_TO_CM2 = 0.389379338e-27
POWERLAW_SCALE_A_CM2 = 5.53e-36
POWERLAW_SCALE_B = 0.363

CTW_NU_CC_TABLE = np.array(
    [
        [1.0e4, 0.48e-34],
        [2.5e4, 0.93e-34],
        [6.0e4, 0.16e-33],
        [1.0e5, 0.22e-33],
        [2.5e5, 0.36e-33],
        [6.0e5, 0.56e-33],
        [1.0e6, 0.72e-33],
        [2.5e6, 0.11e-32],
        [6.0e6, 0.16e-32],
        [1.0e7, 0.20e-32],
        [2.5e7, 0.29e-32],
        [6.0e7, 0.40e-32],
        [1.0e8, 0.48e-32],
        [2.5e8, 0.67e-32],
        [6.0e8, 0.91e-32],
        [1.0e9, 0.11e-31],
        [2.5e9, 0.14e-31],
        [6.0e9, 0.19e-31],
        [1.0e10, 0.22e-31],
        [2.5e10, 0.29e-31],
        [6.0e10, 0.37e-31],
        [1.0e11, 0.43e-31],
        [2.5e11, 0.56e-31],
        [6.0e11, 0.72e-31],
        [1.0e12, 0.83e-31],
    ],
    dtype=float,
)


@dataclass(frozen=True)
class Regime:
    name: str
    rho0: float
    r0: float
    sigma_r: float
    h_over_r: float
    profile: str
    radial_power: float
    funnel_depletion: float
    funnel_theta: float
    envelope_rho0: float
    envelope_alpha: float
    r_min: float
    r_max: float
    rho_floor: float


@dataclass(frozen=True)
class SigmaModel:
    name: str
    path: Path


REGIMES = [
    Regime(
        "fiducial_uhe_default",
        1.0e-2,
        10.0,
        5.0,
        0.25,
        "gaussian",
        2.0,
        0.0,
        15.0,
        0.0,
        2.5,
        4.0,
        60.0,
        1.0e-99,
    ),
    Regime(
        "collapsar_ndaf_like",
        3.0e10,
        12.0,
        8.0,
        0.45,
        "collapsar_ndaf_like",
        1.7,
        0.85,
        20.0,
        3.0e8,
        2.2,
        3.0,
        90.0,
        1.0e-20,
    ),
]


def ensure_dirs() -> None:
    for path in [OUTPUT_DIR, PLOT_DIR, IMAGE_DIR, RAY_DIR, DATA_DIR]:
        path.mkdir(parents=True, exist_ok=True)


def load_sigma_table(path: Path) -> np.ndarray:
    data = np.loadtxt(path, comments="#")
    if data.ndim != 2 or data.shape[1] < 3:
        raise ValueError(f"{path} is not a three-column sigma table")
    return data[:, :3]


def log_interp(x: np.ndarray, xp: np.ndarray, fp: np.ndarray) -> np.ndarray:
    return np.exp(np.interp(np.log(x), np.log(xp), np.log(fp)))


def write_literature_powerlaw_scale(grid: np.ndarray) -> Path:
    out = DATA_DIR / "sigma_nuN_CC_literature_powerlaw_scale.dat"
    sigma_cm2 = POWERLAW_SCALE_A_CM2 * np.power(grid, POWERLAW_SCALE_B)
    sigma_gev2 = sigma_cm2 / GEV_MINUS2_TO_CM2
    with out.open("w") as f:
        f.write("# model literature_powerlaw_scale\n")
        f.write("# channel charged_current_neutrino_nucleon\n")
        f.write("# reference Gandhi_Quigg_Reno_Sarcevic_style_powerlaw\n")
        f.write("# formula sigma_cm2 = A * (E_GeV)^B\n")
        f.write(f"# A_cm2 {POWERLAW_SCALE_A_CM2:.8e}\n")
        f.write(f"# B {POWERLAW_SCALE_B:.8e}\n")
        f.write("# notes Approximate literature-scale curve retained only as a scale check; not a PDF reference table.\n")
        f.write("# Enu_GeV sigma_GeV_minus2 sigma_cm2\n")
        for e, sg, sc in zip(grid, sigma_gev2, sigma_cm2):
            f.write(f"{e:.10e} {sg:.10e} {sc:.10e}\n")
    return out


def write_ctw_reference() -> Path:
    out = DATA_DIR / "sigma_nuN_CC_CTW_reference.dat"
    energies = CTW_NU_CC_TABLE[:, 0]
    sigma_cm2 = CTW_NU_CC_TABLE[:, 1]
    sigma_gev2 = sigma_cm2 / GEV_MINUS2_TO_CM2
    with out.open("w") as f:
        f.write("# model CTW_reference\n")
        f.write("# channel charged_current_neutrino_nucleon\n")
        f.write("# neutrino_status neutrino_not_antineutrino\n")
        f.write("# target isoscalar_nucleon\n")
        f.write("# source Connolly_Thorne_Waters_PhysRevD83_113009_Table_I\n")
        f.write("# arxiv https://arxiv.org/abs/1102.0691\n")
        f.write("# pdf_set MSTW_2008\n")
        f.write("# energy_range_GeV 1e4 1e12\n")
        f.write("# values table_published_central_values\n")
        f.write("# notes Table I nuN CC central cross sections; not NC and not antineutrino.\n")
        f.write("# Enu_GeV sigma_GeV_minus2 sigma_cm2\n")
        for e, sg, sc in zip(energies, sigma_gev2, sigma_cm2):
            f.write(f"{e:.10e} {sg:.10e} {sc:.10e}\n")
    return out


def audit_tables(models: list[SigmaModel]) -> None:
    lines = [
        "# DIS Table Audit",
        "",
        "Astrophysical model is frozen for this validation. The only quantity varied is the neutrino-nucleon charged-current cross-section table passed to `SigmaTable`.",
        "",
        "## Interpolation and extrapolation",
        "",
        "- Loader: `src/sigma_table.cpp`.",
        "- Columns read: `Enu_GeV`, `sigma_GeV_minus2`, `sigma_cm2`.",
        "- Interpolation: logarithmic in energy and logarithmic in cross section.",
        "- Extrapolation: disabled. Requests outside the table range throw `Requested energy outside sigma table range.`",
        "- The tables encode precomputed total CC values in one scalar column; they do not carry separate structure-function columns.",
        "- `use_F3` is recorded in image metadata, but xF3 inclusion cannot be inferred from the sigma table alone.",
        "",
        "## Tables",
        "",
        "| model | filename | rows | Emin [GeV] | Emax [GeV] | sigma(Emin) [cm2] | sigma(Emax) [cm2] | channel/status |",
        "|---|---|---:|---:|---:|---:|---:|---|",
    ]
    for model in models:
        data = load_sigma_table(model.path)
        lines.append(
            f"| {model.name} | `{model.path.relative_to(ROOT)}` | {len(data)} | "
            f"{data[0,0]:.6e} | {data[-1,0]:.6e} | {data[0,2]:.6e} | {data[-1,2]:.6e} | CC nuN scalar table |"
        )
    lines.extend(
        [
            "",
            "## Reference-model provenance",
            "",
            "`CTW_reference` is the published central charged-current nuN table from Connolly, Thorne & Waters, Phys. Rev. D 83, 113009, Table I:",
            "",
            "- URL: https://arxiv.org/abs/1102.0691",
            "- PDF set: MSTW 2008.",
            "- Status: neutrino, not antineutrino.",
            "- Channel: charged current only.",
            "- Energy range: 1e4 <= E_nu/GeV <= 1e12.",
            "- Units: cm^2, converted internally to GeV^-2 for the second table column.",
            "- Values: published table values, not digitized from a plot.",
            "",
            "`literature_powerlaw_scale` is retained only as an old scale-check curve:",
            "",
            "```text",
            f"sigma_CC(E) = {POWERLAW_SCALE_A_CM2:.3e} * (E_GeV)^{POWERLAW_SCALE_B:.3f} cm^2",
            "```",
            "",
            "It must not be called `PDF_reference` in paper text.",
            "",
            "Reference context:",
            "",
            "- Gandhi, Quigg, Reno & Sarcevic, Phys. Rev. D 58, 093009; UHE neutrino interactions over roughly 1e9-1e21 eV: https://arxiv.org/abs/hep-ph/9807264",
            "- Connolly, Thorne & Waters, Phys. Rev. D 83, 113009; MSTW PDF-based CC/NC cross sections and parametrizations over 1e4-1e12 GeV: https://arxiv.org/abs/1102.0691",
            "- Cooper-Sarkar, Mertsch & Sarkar, JHEP 08, 042; Standard Model high-energy neutrino cross sections and PDF uncertainties: https://arxiv.org/abs/1106.3723",
        ]
    )
    (OUTPUT_DIR / "dis_table_audit.md").write_text("\n".join(lines) + "\n")


def write_reference_model_audit(models: list[SigmaModel]) -> None:
    ctw = load_sigma_table(next(model.path for model in models if model.name == "CTW_reference"))
    lines = [
        "# Reference Model Audit",
        "",
        "## Implemented paper-level reference",
        "",
        "- Model name: `CTW_reference`.",
        "- Exact source: Connolly, Thorne & Waters, Phys. Rev. D 83, 113009, Table I.",
        "- Citation: Amy Connolly, Robert S. Thorne, David Waters, Calculation of High Energy Neutrino-Nucleon Cross Sections and Uncertainties Using the MSTW Parton Distribution Functions and Implications for Future Experiments.",
        "- URL: https://arxiv.org/abs/1102.0691",
        "- Values used: table values, not plot digitization.",
        "- Energy range: `1e4 <= E_nu/GeV <= 1e12`.",
        "- Units in source: `cm^2`.",
        "- Units in HADROS table: `Enu_GeV`, `sigma_GeV_minus2`, `sigma_cm2`.",
        "- Channel: charged current.",
        "- Particle: neutrino, not antineutrino.",
        "- Target: isoscalar nucleon.",
        "- PDF set/theory: MSTW 2008 PDFs, Standard Model DIS calculation.",
        "",
        "## Implemented data range",
        "",
        f"- Rows: {len(ctw)}.",
        f"- Minimum energy: {ctw[0,0]:.6e} GeV.",
        f"- Maximum energy: {ctw[-1,0]:.6e} GeV.",
        f"- sigma_CC(Emin): {ctw[0,2]:.6e} cm^2.",
        f"- sigma_CC(Emax): {ctw[-1,2]:.6e} cm^2.",
        "",
        "## Limitations",
        "",
        "- This is the CTW central-value table only; it does not include CTW PDF uncertainty bands.",
        "- It is CC-only; NC is not included in the opacity comparison.",
        "- It is neutrino-only; antineutrino values are tabulated separately in CTW Table II and are not used here.",
        "- HADROS uses log-log interpolation between the published table points and no extrapolation outside the table range.",
        "- Comparisons involving `CTW_reference` are restricted to the common energy range of all active tables.",
        "",
        "## Old approximate scale curve",
        "",
        "- Model name: `literature_powerlaw_scale`.",
        "- Formula: `sigma_CC(E) = 5.53e-36 * E_GeV^0.363 cm^2`.",
        "- Purpose: legacy scale check only.",
        "- It must not be described as `PDF_reference` or as a paper-level PDF validation.",
    ]
    (OUTPUT_DIR / "reference_model_audit.md").write_text("\n".join(lines) + "\n")


def plot_sigma_comparison(models: list[SigmaModel]) -> None:
    loaded = {model.name: load_sigma_table(model.path) for model in models}
    emin = max(data[0, 0] for data in loaded.values())
    emax = min(data[-1, 0] for data in loaded.values())
    energies = np.logspace(np.log10(emin), np.log10(emax), 240)
    sigma = {
        name: log_interp(energies, data[:, 0], data[:, 2])
        for name, data in loaded.items()
    }

    with (OUTPUT_DIR / "sigma_model_comparison.csv").open("w", newline="") as f:
        writer = csv.writer(f)
        model_names = [model.name for model in models]
        header = ["E_GeV"]
        header.extend(f"sigma_{name}_cm2" for name in model_names)
        header.extend(f"{name}_over_CTW_reference" for name in model_names if name != "CTW_reference")
        header.append("IIM_over_GBW")
        writer.writerow(header)
        for i, e in enumerate(energies):
            row = [f"{e:.10e}"]
            row.extend(f"{sigma[name][i]:.10e}" for name in model_names)
            row.extend(
                f"{sigma[name][i] / sigma['CTW_reference'][i]:.10e}"
                for name in model_names
                if name != "CTW_reference"
            )
            row.append(f"{sigma['IIM'][i] / sigma['GBW'][i]:.10e}")
            writer.writerow(row)

    plt.figure(figsize=(7.2, 5.0))
    for name, values in sigma.items():
        plt.loglog(energies, values, label=name)
    plt.xlabel(r"$E_\nu$ [GeV]")
    plt.ylabel(r"$\sigma_{\nu N}^{CC}$ [cm$^2$]")
    plt.title("DIS model cross-section comparison")
    plt.grid(True, which="both", alpha=0.25)
    plt.legend()
    plt.tight_layout()
    plt.savefig(PLOT_DIR / "sigma_nuN_model_comparison.png", dpi=180)
    plt.close()

    plt.figure(figsize=(7.2, 5.0))
    plt.semilogx(energies, sigma["GBW"] / sigma["CTW_reference"], label="GBW / CTW_reference")
    plt.semilogx(energies, sigma["IIM"] / sigma["CTW_reference"], label="IIM / CTW_reference")
    plt.semilogx(
        energies,
        sigma["literature_powerlaw_scale"] / sigma["CTW_reference"],
        label="literature_powerlaw_scale / CTW_reference",
    )
    plt.semilogx(energies, sigma["IIM"] / sigma["GBW"], label="IIM / GBW")
    plt.axhline(1.0, color="0.3", lw=0.8)
    plt.xlabel(r"$E_\nu$ [GeV]")
    plt.ylabel("cross-section ratio")
    plt.title("DIS model cross-section ratios")
    plt.grid(True, which="both", alpha=0.25)
    plt.legend()
    plt.tight_layout()
    plt.savefig(PLOT_DIR / "sigma_nuN_model_ratios.png", dpi=180)
    plt.close()


def run_command(args: list[str], env: dict[str, str]) -> str:
    proc = subprocess.run(
        args,
        cwd=ROOT,
        env=env,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        check=True,
    )
    return proc.stdout


def output_path_from_stdout(stdout: str) -> Path:
    match = re.search(r"Saved:\s*(\S+)", stdout)
    if not match:
        raise RuntimeError(f"Could not find output image path in stdout:\n{stdout}")
    return ROOT / match.group(1)


def make_geodesic_cache(env: dict[str, str]) -> Path:
    cache = RAY_DIR / "kerr_geodesics_dis_validation.bin"
    run_command(
        [
            "./compute_kerr_geodesics",
            "0.0001",
            "60.0",
            "80.0",
            "25.0",
            "32",
            "32",
            "80.0",
            "0.01",
            str(cache.relative_to(ROOT)),
        ],
        env,
    )
    return cache


def run_image(regime: Regime, model: SigmaModel, energy: float, cache: Path, env: dict[str, str]) -> Path:
    stdout = run_command(
        [
            "./compute_kerr_image_from_cache",
            f"{energy:.10e}",
            "0.0001",
            "3.0",
            f"{regime.rho0:.10e}",
            f"{regime.r0:.10e}",
            f"{regime.sigma_r:.10e}",
            f"{regime.h_over_r:.10e}",
            "3.5",
            "1.0",
            "15.0",
            "2.0",
            "1.0e12",
            "1.0",
            "10.0",
            "1.0",
            "80.0",
            str(model.path.relative_to(ROOT)),
            regime.profile,
            f"{regime.radial_power:.10e}",
            f"{regime.funnel_depletion:.10e}",
            f"{regime.funnel_theta:.10e}",
            f"{regime.envelope_rho0:.10e}",
            f"{regime.envelope_alpha:.10e}",
            f"{regime.r_min:.10e}",
            f"{regime.r_max:.10e}",
            f"{regime.rho_floor:.10e}",
            "60.0",
            "1",
            str(cache.relative_to(ROOT)),
            "inner_ring",
            "20.0",
            "1.0",
            "2.0",
            "0.1",
            "1.0",
            f"{regime.rho0:.10e}",
            "0.0",
            "1.0e2",
        ],
        env,
    )
    return output_path_from_stdout(stdout)


def load_image(path: Path) -> np.ndarray:
    data = np.loadtxt(path, comments="#")
    if data.ndim == 1:
        data = data.reshape(1, -1)
    return data


def metrics_from_image(path: Path) -> dict[str, float]:
    data = load_image(path)
    tau = data[:, 4]
    psurv = data[:, 5]
    intensity = data[:, 6]
    captured = data[:, 7]
    valid = captured < 0.5
    if np.any(valid):
        tau_valid = tau[valid]
        psurv_valid = psurv[valid]
    else:
        tau_valid = tau
        psurv_valid = psurv
    return {
        "mean_tau": float(np.mean(tau_valid)),
        "max_tau": float(np.max(tau_valid)),
        "mean_Psurv": float(np.mean(psurv_valid)),
        "total_intensity": float(np.sum(intensity)),
        "valid_rays": int(np.count_nonzero(valid)),
    }


def run_observable_suite(models: list[SigmaModel]) -> list[dict[str, str]]:
    env = os.environ.copy()
    env["OMP_NUM_THREADS"] = env.get("NTHREADS", "2")
    cache = make_geodesic_cache(env)
    rows: list[dict[str, str]] = []
    energy = 1.0e11
    for regime in REGIMES:
        for model in models:
            image = run_image(regime, model, energy, cache, env)
            metrics = metrics_from_image(image)
            row = {
                "regime": regime.name,
                "DIS_model": model.name,
                "E_GeV": f"{energy:.6e}",
                "image_path": str(image.relative_to(ROOT)),
                **{k: f"{v:.10e}" if isinstance(v, float) else str(v) for k, v in metrics.items()},
            }
            rows.append(row)
    with (OUTPUT_DIR / "dis_observable_comparison.csv").open("w", newline="") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=[
                "regime",
                "DIS_model",
                "E_GeV",
                "mean_tau",
                "max_tau",
                "mean_Psurv",
                "total_intensity",
                "valid_rays",
                "image_path",
            ],
        )
        writer.writeheader()
        writer.writerows(rows)
    return rows


def image_grid(path: Path, column: int = 6) -> np.ndarray:
    data = load_image(path)
    i = data[:, 0].astype(int)
    j = data[:, 1].astype(int)
    nx = int(i.max()) + 1
    ny = int(j.max()) + 1
    grid = np.full((ny, nx), np.nan)
    grid[j, i] = data[:, column]
    return grid


def plot_image_comparison(rows: list[dict[str, str]]) -> None:
    row_order = ["fiducial_uhe_default", "collapsar_ndaf_like"]
    col_order = ["GBW", "IIM", "literature_powerlaw_scale", "CTW_reference"]
    fig, axes = plt.subplots(len(row_order), len(col_order), figsize=(13.5, 6.8))
    for r, regime in enumerate(row_order):
        images = []
        for model in col_order:
            entry = next(row for row in rows if row["regime"] == regime and row["DIS_model"] == model)
            images.append(image_grid(ROOT / entry["image_path"], 6))
        vmax = max(float(np.nanmax(img)) for img in images)
        for c, (model, img) in enumerate(zip(col_order, images)):
            ax = axes[r, c]
            shown = np.sqrt(np.maximum(img, 0.0))
            norm_vmax = math.sqrt(vmax) if vmax > 0 else 1.0
            ax.imshow(shown, origin="lower", cmap="magma", vmin=0.0, vmax=norm_vmax)
            ax.set_title(f"{regime}\n{model}", fontsize=9)
            ax.set_xticks([])
            ax.set_yticks([])
    fig.suptitle("UHE image comparison across DIS models")
    fig.tight_layout()
    fig.savefig(PLOT_DIR / "uhe_image_dis_model_comparison.png", dpi=180)
    plt.close(fig)

    fig, axes = plt.subplots(len(row_order), 3, figsize=(10.5, 6.8))
    for r, regime in enumerate(row_order):
        ref_entry = next(row for row in rows if row["regime"] == regime and row["DIS_model"] == "CTW_reference")
        ref = image_grid(ROOT / ref_entry["image_path"], 6)
        for c, model in enumerate(["GBW", "IIM", "literature_powerlaw_scale"]):
            entry = next(row for row in rows if row["regime"] == regime and row["DIS_model"] == model)
            img = image_grid(ROOT / entry["image_path"], 6)
            ratio = np.divide(img, ref, out=np.ones_like(img), where=ref > 0)
            ax = axes[r, c]
            im = ax.imshow(np.clip(ratio, 0.0, 2.0), origin="lower", cmap="coolwarm", vmin=0.0, vmax=2.0)
            ax.set_title(f"{regime}\n{model} / CTW_reference", fontsize=9)
            ax.set_xticks([])
            ax.set_yticks([])
    fig.colorbar(im, ax=axes.ravel().tolist(), shrink=0.75, label="intensity ratio")
    fig.suptitle("UHE image ratios relative to CTW_reference")
    fig.savefig(PLOT_DIR / "uhe_image_dis_model_ratio.png", dpi=180)
    plt.close(fig)


def plot_spectral_propagation(rows: list[dict[str, str]], models: list[SigmaModel]) -> None:
    loaded = {model.name: load_sigma_table(model.path) for model in models}
    energies = np.logspace(5, 12, 64)
    phi = energies ** -2.0
    fig, axes = plt.subplots(1, 2, figsize=(11.0, 4.5), sharex=True)
    fig2, axes2 = plt.subplots(1, 2, figsize=(11.0, 4.5), sharex=True, sharey=True)
    for panel, regime in enumerate(["fiducial_uhe_default", "collapsar_ndaf_like"]):
        for model in models:
            entry = next(row for row in rows if row["regime"] == regime and row["DIS_model"] == model.name)
            data = load_image(ROOT / entry["image_path"])
            tau_ref = data[:, 4]
            intensity_ref = data[:, 6]
            sigma_ref = log_interp(np.array([float(entry["E_GeV"])]), loaded[model.name][:, 0], loaded[model.name][:, 2])[0]
            columns = tau_ref / sigma_ref
            emitted = np.divide(
                intensity_ref,
                np.exp(-np.clip(tau_ref, 0.0, 700.0)),
                out=np.zeros_like(intensity_ref),
                where=intensity_ref > 0,
            )
            emitted_sum = float(np.sum(emitted))
            observed = []
            attenuation = []
            for energy in energies:
                sigma_e = log_interp(np.array([energy]), loaded[model.name][:, 0], loaded[model.name][:, 2])[0]
                tau_e = columns * sigma_e
                obs = float(np.sum(emitted * np.exp(-np.clip(tau_e, 0.0, 700.0))))
                observed.append(phi[len(observed)] * obs)
                attenuation.append(obs / emitted_sum if emitted_sum > 0 else 0.0)
            observed_arr = np.asarray(observed)
            if np.any(observed_arr > 0):
                observed_plot = observed_arr
            else:
                observed_plot = np.full_like(observed_arr, 1.0e-300)
            axes[panel].loglog(energies, observed_plot, label=model.name)
            axes2[panel].semilogx(energies, attenuation, label=model.name)
        axes[panel].set_title(regime)
        axes[panel].set_xlabel(r"$E_\nu$ [GeV]")
        axes[panel].set_ylabel("observed spectrum proxy")
        axes[panel].grid(True, which="both", alpha=0.25)
        axes[panel].legend(fontsize=8)
        plotted_lines = [line.get_ydata() for line in axes[panel].lines]
        if plotted_lines and all(np.nanmax(y) <= 1.0e-299 for y in plotted_lines):
            axes[panel].text(
                0.05,
                0.08,
                "all observed values are zero;\nshown at numerical floor",
                transform=axes[panel].transAxes,
                fontsize=8,
                color="0.25",
            )
        axes2[panel].set_title(regime)
        axes2[panel].set_xlabel(r"$E_\nu$ [GeV]")
        axes2[panel].set_ylabel("observed / emitted")
        axes2[panel].grid(True, which="both", alpha=0.25)
        axes2[panel].legend(fontsize=8)
    fig.suptitle("Observed spectrum propagated through DIS model choices")
    fig.tight_layout()
    fig.savefig(PLOT_DIR / "observed_spectrum_dis_models.png", dpi=180)
    plt.close(fig)
    fig2.suptitle("Attenuation ratio across DIS model choices")
    fig2.tight_layout()
    fig2.savefig(PLOT_DIR / "attenuation_ratio_dis_models.png", dpi=180)
    plt.close(fig2)


def write_summary(rows: list[dict[str, str]]) -> None:
    by_regime = {}
    for row in rows:
        by_regime.setdefault(row["regime"], {})[row["DIS_model"]] = row
    lines = [
        "# DIS Model Validation Summary",
        "",
        "The astrophysical model was frozen inside each regime. Only the DIS cross-section table was changed.",
        "",
        "## Observable comparison at E = 1e11 GeV",
        "",
        "| regime | model | mean_tau | max_tau | mean_Psurv | total_intensity | valid_rays |",
        "|---|---|---:|---:|---:|---:|---:|",
    ]
    for row in rows:
        lines.append(
            f"| {row['regime']} | {row['DIS_model']} | {float(row['mean_tau']):.6e} | "
            f"{float(row['max_tau']):.6e} | {float(row['mean_Psurv']):.6e} | "
            f"{float(row['total_intensity']):.6e} | {row['valid_rays']} |"
        )
    lines.extend(["", "## Interpretation", ""])
    for regime, model_rows in by_regime.items():
        ctw_tau = float(model_rows["CTW_reference"]["mean_tau"])
        gbw_tau = float(model_rows["GBW"]["mean_tau"])
        iim_tau = float(model_rows["IIM"]["mean_tau"])
        scale_tau = float(model_rows["literature_powerlaw_scale"]["mean_tau"])
        lines.append(
            f"- `{regime}`: mean_tau ratios are GBW/CTW={gbw_tau / ctw_tau:.3e}, "
            f"IIM/CTW={iim_tau / ctw_tau:.3e}, and "
            f"literature_powerlaw_scale/CTW={scale_tau / ctw_tau:.3e}."
        )
    lines.extend(
        [
            "",
            "## Safe claims",
            "",
            "- The same ray cache, density background, source prescription, observer, spin, and energy were used inside each regime.",
            "- Differences in tau and P_surv in this suite are therefore DIS-table effects.",
            "- The CTW comparison can support statements about sensitivity to the assumed UHE neutrino-nucleon cross section relative to a published MSTW-2008 Standard Model DIS calculation.",
            "",
            "## Claims to avoid",
            "",
            "- Do not present `literature_powerlaw_scale` as a PDF reference.",
            "- Do not claim CTW PDF uncertainty-band coverage; this implementation uses the central Table I values only.",
            "- Do not claim CC+NC or neutrino+antineutrino coverage from these tables; the propagation files used here are charged-current scalar tables.",
            "- Do not claim xF3 status from the table alone; it is not represented as a separate column.",
        ]
    )
    (OUTPUT_DIR / "dis_model_validation_summary.md").write_text("\n".join(lines) + "\n")


def main() -> None:
    ensure_dirs()
    gbw = load_sigma_table(DATA_DIR / "sigma_nuN_CC_GBW.dat")
    powerlaw_path = write_literature_powerlaw_scale(gbw[:, 0])
    ctw_path = write_ctw_reference()
    models = [
        SigmaModel("GBW", DATA_DIR / "sigma_nuN_CC_GBW.dat"),
        SigmaModel("IIM", DATA_DIR / "sigma_nuN_CC_IIM.dat"),
        SigmaModel("literature_powerlaw_scale", powerlaw_path),
        SigmaModel("CTW_reference", ctw_path),
    ]
    audit_tables(models)
    write_reference_model_audit(models)
    plot_sigma_comparison(models)
    rows = run_observable_suite(models)
    plot_spectral_propagation(rows, models)
    plot_image_comparison(rows)
    write_summary(rows)
    print(f"Wrote {OUTPUT_DIR.relative_to(ROOT)}")
    print(f"Wrote {PLOT_DIR.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
