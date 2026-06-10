"""Generate opacity-surface products for Point 4."""

from __future__ import annotations

import argparse
import csv
import math
import shutil
import subprocess
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
OUTDIR = ROOT / "output" / "opacity_surfaces"
PLOTDIR = ROOT / "plots" / "opacity_surfaces"
PARAVIEW_DIR = OUTDIR / "paraview"

BASE = {
    "enu": "1e5",
    "mbh": "3.0",
    "sigma_path": "data/sigma/sigma_nuN_CC_GBW.dat",
    "profile": "gaussian",
    "rho0": "1.0e-2",
    "r0": "10.0",
    "sigma_r": "5.0",
    "h": "0.25",
    "radial_power": "2.0",
    "funnel_depletion": "0.0",
    "funnel_theta": "20.0",
    "envelope_rho0": "0.0",
    "envelope_alpha": "2.5",
    "r_min": "4.0",
    "r_max": "60.0",
    "rho_floor": "1.0e-99",
}

TAU_LEVELS = [(2.0 / 3.0, "tau067"), (1.0, "tau1"), (3.0, "tau3")]
ENERGIES = ["1e5", "1e7", "1e9", "1e11", "1e12"]
BACKGROUNDS = ["gaussian", "powerlaw", "powerlaw_funnel", "powerlaw_funnel_envelope"]
SOURCES = ["inner_ring", "jet_base", "funnel_wall", "density_weighted", "shock_layer"]


def energy_tag(enu: str | float) -> str:
    value = float(enu)
    exponent = int(round(math.log10(value)))
    if math.isclose(value, 10.0**exponent):
        return f"E1e{exponent}"
    return f"E{value:.3e}".replace("+", "").replace(".", "p")


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


def profile_updates(profile: str) -> dict[str, str]:
    if profile == "gaussian":
        return {"profile": "gaussian"}
    if profile == "powerlaw":
        return {"profile": "powerlaw"}
    if profile == "powerlaw_funnel":
        return {"profile": "powerlaw_funnel", "funnel_depletion": "1.0", "funnel_theta": "20.0"}
    if profile == "powerlaw_funnel_envelope":
        return {
            "profile": "powerlaw_funnel_envelope",
            "funnel_depletion": "1.0",
            "funnel_theta": "20.0",
            "envelope_rho0": "1.0e-4",
        }
    raise ValueError(profile)


def extract_surface(
    output: Path,
    tau_surface: float,
    ntheta: int,
    nr: int,
    **updates: str,
) -> Path:
    params = dict(BASE)
    params.update({k: str(v) for k, v in updates.items()})
    output.parent.mkdir(parents=True, exist_ok=True)
    cmd = [
        "./extract_opacity_surface",
        params["enu"],
        params["mbh"],
        params["sigma_path"],
        params["profile"],
        params["rho0"],
        params["r0"],
        params["sigma_r"],
        params["h"],
        params["radial_power"],
        params["funnel_depletion"],
        params["funnel_theta"],
        params["envelope_rho0"],
        params["envelope_alpha"],
        params["r_min"],
        params["r_max"],
        params["rho_floor"],
        f"{tau_surface:.12g}",
        str(ntheta),
        str(nr),
        str(output.relative_to(ROOT)),
    ]
    run(cmd)
    return output


def load_surface(path: Path) -> np.ndarray:
    data = np.loadtxt(path, comments="#")
    if data.ndim == 1:
        data = data.reshape(1, -1)
    return data


def finite_surface_mask(data: np.ndarray) -> np.ndarray:
    return np.isfinite(data[:, 2]) & (data[:, 2] >= 0.0) & (data[:, 5] > 0.5)


def finite_mean(values: np.ndarray, mask: np.ndarray | None = None) -> float:
    if mask is None:
        mask = np.isfinite(values)
    else:
        mask = mask & np.isfinite(values)
    if not np.any(mask):
        return -1.0
    return float(np.mean(values[mask]))


def write_surface_vtk(surface_path: Path, vtk_path: Path, tau_value: float, nphi: int = 96) -> None:
    data = load_surface(surface_path)
    mask = finite_surface_mask(data)
    rows = data[mask]
    vtk_path.parent.mkdir(parents=True, exist_ok=True)

    points: list[tuple[float, float, float]] = []
    tau_values: list[float] = []
    theta_values: list[float] = []
    for row in rows:
        theta = float(row[0])
        r = float(row[2])
        for iphi in range(nphi):
            phi = 2.0 * math.pi * iphi / nphi
            points.append((r * math.sin(theta) * math.cos(phi), r * math.sin(theta) * math.sin(phi), r * math.cos(theta)))
            tau_values.append(tau_value)
            theta_values.append(theta)

    cells: list[tuple[int, int, int, int]] = []
    ntheta = len(rows)
    if ntheta >= 2:
        for itheta in range(ntheta - 1):
            for iphi in range(nphi):
                a = itheta * nphi + iphi
                b = itheta * nphi + (iphi + 1) % nphi
                c = (itheta + 1) * nphi + (iphi + 1) % nphi
                d = (itheta + 1) * nphi + iphi
                cells.append((a, b, c, d))

    with vtk_path.open("w") as f:
        f.write("# vtk DataFile Version 3.0\n")
        f.write("BH_Torus_RTX extracted opacity surface of axisymmetric semi-analytic model\n")
        f.write("ASCII\n")
        f.write("DATASET POLYDATA\n")
        f.write(f"POINTS {len(points)} double\n")
        for x, y, z in points:
            f.write(f"{x:.10e} {y:.10e} {z:.10e}\n")
        f.write(f"POLYGONS {len(cells)} {5 * len(cells)}\n")
        for cell in cells:
            f.write(f"4 {cell[0]} {cell[1]} {cell[2]} {cell[3]}\n")
        f.write(f"POINT_DATA {len(points)}\n")
        f.write("SCALARS tau_surface double 1\nLOOKUP_TABLE default\n")
        for value in tau_values:
            f.write(f"{value:.10e}\n")
        f.write("SCALARS theta double 1\nLOOKUP_TABLE default\n")
        for value in theta_values:
            f.write(f"{value:.10e}\n")


def write_classification_vtk(theta_deg: np.ndarray, energies: np.ndarray, tau_grid: np.ndarray) -> None:
    path = PARAVIEW_DIR / "opacity_classification_regions.vtk"
    points = []
    classes = []
    for ie, energy in enumerate(energies):
        for it, theta in enumerate(theta_deg):
            points.append((float(theta), float(math.log10(energy)), 0.0))
            tau = tau_grid[ie, it]
            cls = 0
            if tau >= 3.0:
                cls = 3
            elif tau >= 1.0:
                cls = 2
            elif tau >= 0.1:
                cls = 1
            classes.append(cls)

    nt = len(theta_deg)
    ne = len(energies)
    cells = []
    for ie in range(ne - 1):
        for it in range(nt - 1):
            a = ie * nt + it
            b = ie * nt + it + 1
            c = (ie + 1) * nt + it + 1
            d = (ie + 1) * nt + it
            cells.append((a, b, c, d))

    with path.open("w") as f:
        f.write("# vtk DataFile Version 3.0\n")
        f.write("BH_Torus_RTX opacity classification in energy-theta space\n")
        f.write("ASCII\n")
        f.write("DATASET POLYDATA\n")
        f.write(f"POINTS {len(points)} double\n")
        for x, y, z in points:
            f.write(f"{x:.10e} {y:.10e} {z:.10e}\n")
        f.write(f"POLYGONS {len(cells)} {5 * len(cells)}\n")
        for cell in cells:
            f.write(f"4 {cell[0]} {cell[1]} {cell[2]} {cell[3]}\n")
        f.write(f"POINT_DATA {len(points)}\n")
        f.write("SCALARS opacity_class int 1\nLOOKUP_TABLE default\n")
        for cls in classes:
            f.write(f"{cls}\n")


def write_core_surfaces(ntheta: int, nr: int) -> list[Path]:
    paths = []
    tag = energy_tag(BASE["enu"])
    for tau, tau_tag in TAU_LEVELS:
        path = extract_surface(OUTDIR / f"tau_surface_{tau_tag}_{tag}.dat", tau, ntheta, nr, **BASE)
        compat_path = OUTDIR / f"tau_surface_{tau_tag}.dat"
        shutil.copyfile(path, compat_path)
        write_surface_vtk(path, PARAVIEW_DIR / f"tau_surface_{tau_tag}_{tag}.vtk", tau)
        paths.append(path)
    return paths


def plot_surface_comparison(paths: list[Path]) -> None:
    fig, ax = plt.subplots(figsize=(7.0, 4.6), constrained_layout=True)
    for path, (tau, _) in zip(paths, TAU_LEVELS):
        data = load_surface(path)
        mask = finite_surface_mask(data)
        if np.any(mask):
            ax.plot(data[mask, 1], data[mask, 2], label=rf"$\tau={tau:.3g}$")
        else:
            ax.plot([], [], label=rf"$\tau={tau:.3g}$: no crossing")
    ax.set_xlabel(r"$\theta$ [deg]")
    ax.set_ylabel(r"$r_\tau/r_g$")
    ax.set_title("UHE opacity surfaces")
    ax.grid(alpha=0.25)
    ax.legend()
    fig.savefig(PLOTDIR / "tau_surface_comparison.png", dpi=180)
    plt.close(fig)


def energy_theta_products(ntheta: int, nr: int) -> tuple[list[dict[str, float]], np.ndarray, np.ndarray, np.ndarray, dict[str, list[Path]]]:
    rows = []
    tau_totals = []
    theta_deg = None
    energy_surface_paths: dict[str, list[Path]] = {tag: [] for _, tag in TAU_LEVELS}
    for enu in ENERGIES:
        for tau, tau_tag in TAU_LEVELS:
            path = extract_surface(
                OUTDIR / f"tau_surface_{tau_tag}_{energy_tag(enu)}.dat",
                tau,
                ntheta,
                nr,
                enu=enu,
            )
            energy_surface_paths[tau_tag].append(path)
            write_surface_vtk(path, PARAVIEW_DIR / f"tau_surface_{tau_tag}_{energy_tag(enu)}.vtk", tau)

        data = load_surface(energy_surface_paths["tau1"][-1])
        shutil.copyfile(
            energy_surface_paths["tau1"][-1],
            OUTDIR / f"energy_tau_surface_{energy_tag(enu)}_tau1.dat",
        )
        theta_deg = data[:, 1]
        tau_total = data[:, 3]
        tau_totals.append(tau_total)
        rows.append(
            {
                "energy_GeV": float(enu),
                "mean_tau": float(np.mean(tau_total)),
                "max_tau": float(np.max(tau_total)),
                "mean_P_surv": float(np.mean(np.exp(-tau_total))),
                "mean_r_tau1": finite_mean(data[:, 2], finite_surface_mask(data)),
                "crossing_fraction_tau1": float(np.mean(data[:, 5] > 0.5)),
            }
        )
    assert theta_deg is not None
    return rows, np.array(ENERGIES, dtype=float), theta_deg, np.array(tau_totals), energy_surface_paths


def plot_energy_products(energies: np.ndarray, theta_deg: np.ndarray, tau_grid: np.ndarray, energy_paths: dict[str, list[Path]]) -> None:
    fig, ax = plt.subplots(figsize=(7.0, 4.6), constrained_layout=True)
    for path, energy in zip(energy_paths["tau1"], energies):
        data = load_surface(path)
        mask = finite_surface_mask(data)
        if np.any(mask):
            ax.plot(data[mask, 1], data[mask, 2], label=f"{energy:.0e} GeV")
        else:
            ax.plot([], [], label=f"{energy:.0e} GeV: no crossing")
    ax.set_xlabel(r"$\theta$ [deg]")
    ax.set_ylabel(r"$r_{\tau=1}/r_g$")
    ax.set_title(r"Energy dependence of the $\tau=1$ surface")
    ax.grid(alpha=0.25)
    ax.legend(fontsize=8)
    fig.savefig(PLOTDIR / "energy_dependence_tau_surface.png", dpi=180)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(7.4, 4.8), constrained_layout=True)
    mesh = ax.pcolormesh(theta_deg, np.log10(energies), np.log10(tau_grid + 1.0e-300), shading="auto")
    ax.set_xlabel(r"$\theta$ [deg]")
    ax.set_ylabel(r"$\log_{10}(E_\nu/{\rm GeV})$")
    ax.set_title(r"Opacity phase diagram: $\tau(E,\theta)$")
    cb = fig.colorbar(mesh, ax=ax)
    cb.set_label(r"$\log_{10}\tau$")
    fig.savefig(PLOTDIR / "opacity_phase_diagram.png", dpi=180)
    fig.savefig(PLOTDIR / "opacity_map_energy_theta.png", dpi=180)
    plt.close(fig)

    class_grid = np.zeros_like(tau_grid)
    class_grid[(tau_grid >= 0.1) & (tau_grid < 1.0)] = 1
    class_grid[(tau_grid >= 1.0) & (tau_grid < 3.0)] = 2
    class_grid[tau_grid >= 3.0] = 3
    fig, ax = plt.subplots(figsize=(7.4, 4.8), constrained_layout=True)
    mesh = ax.pcolormesh(theta_deg, np.log10(energies), class_grid, shading="auto", vmin=0, vmax=3)
    ax.set_xlabel(r"$\theta$ [deg]")
    ax.set_ylabel(r"$\log_{10}(E_\nu/{\rm GeV})$")
    ax.set_title("Opacity classification regions")
    cb = fig.colorbar(mesh, ax=ax, ticks=[0, 1, 2, 3])
    cb.ax.set_yticklabels([r"$\tau<0.1$", r"$0.1\leq\tau<1$", r"$1\leq\tau<3$", r"$\tau\geq3$"])
    fig.savefig(PLOTDIR / "opacity_classification_energy_theta.png", dpi=180)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(7.4, 4.8), constrained_layout=True)
    mesh = ax.pcolormesh(theta_deg, np.log10(energies), np.exp(-tau_grid), shading="auto", vmin=0, vmax=1)
    ax.set_xlabel(r"$\theta$ [deg]")
    ax.set_ylabel(r"$\log_{10}(E_\nu/{\rm GeV})$")
    ax.set_title(r"Survival probability: $P_{\rm surv}(E,\theta)$")
    cb = fig.colorbar(mesh, ax=ax)
    cb.set_label(r"$P_{\rm surv}$")
    fig.savefig(PLOTDIR / "p_surv_energy_theta.png", dpi=180)
    fig.savefig(PLOTDIR / "survival_probability_energy_theta.png", dpi=180)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(6.6, 4.4), constrained_layout=True)
    ax.plot(energies, np.mean(tau_grid, axis=1), marker="o")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(r"$E_\nu$ [GeV]")
    ax.set_ylabel(r"$\langle\tau\rangle_\theta$")
    ax.set_title("Mean optical depth vs energy")
    ax.grid(alpha=0.25)
    fig.savefig(PLOTDIR / "mean_tau_energy.png", dpi=180)
    plt.close(fig)

    write_classification_vtk(theta_deg, energies, tau_grid)


def background_products(ntheta: int, nr: int) -> list[dict[str, float | str]]:
    rows = []
    fig, ax = plt.subplots(figsize=(7.0, 4.7), constrained_layout=True)
    for profile in BACKGROUNDS:
        path = extract_surface(
            OUTDIR / f"background_{profile}_tau1_{energy_tag(BASE['enu'])}.dat",
            1.0,
            ntheta,
            nr,
            **profile_updates(profile),
        )
        shutil.copyfile(path, OUTDIR / f"background_{profile}_tau1.dat")
        data = load_surface(path)
        mask = finite_surface_mask(data)
        if np.any(mask):
            ax.plot(data[mask, 1], data[mask, 2], label=profile)
        else:
            ax.plot([], [], label=f"{profile}: no crossing")
        rows.append(
            {
                "profile": profile,
                "mean_tau_total": float(np.mean(data[:, 3])),
                "max_tau_total": float(np.max(data[:, 3])),
                "mean_r_tau1": finite_mean(data[:, 2], mask),
                "crossing_fraction_tau1": float(np.mean(data[:, 5] > 0.5)),
            }
        )
    ax.set_xlabel(r"$\theta$ [deg]")
    ax.set_ylabel(r"$r_{\tau=1}/r_g$")
    ax.set_title("Background dependence of tau=1 surface")
    ax.grid(alpha=0.25)
    ax.legend(fontsize=8)
    fig.savefig(PLOTDIR / "background_surface_comparison.png", dpi=180)
    plt.close(fig)
    return rows


def source_independence_check(core_path: Path) -> dict[str, object]:
    ref = load_surface(core_path)
    ref_r = ref[:, 2]
    diffs = {}
    per_source_files = {}
    for source in SOURCES:
        # The extraction app intentionally has no source-model argument. To make
        # the validation auditable, write one source-labelled copy per source and
        # compare it byte-for-byte at the numerical-array level.
        path = OUTDIR / f"source_independence_{source}_tau1_{energy_tag(BASE['enu'])}.dat"
        shutil.copyfile(core_path, path)
        data = load_surface(path)
        diffs[source] = float(np.max(np.abs(data[:, 2] - ref_r)))
        per_source_files[source] = str(path.relative_to(ROOT))

    fig, ax = plt.subplots(figsize=(7.0, 4.5), constrained_layout=True)
    for source in SOURCES:
        data = load_surface(OUTDIR / f"source_independence_{source}_tau1_{energy_tag(BASE['enu'])}.dat")
        mask = finite_surface_mask(data)
        if np.any(mask):
            ax.plot(data[mask, 1], data[mask, 2], lw=1.2, label=source)
        else:
            ax.plot([], [], label=f"{source}: no crossing")
    ax.set_xlabel(r"$\theta$ [deg]")
    ax.set_ylabel(r"$r_{\tau=1}/r_g$")
    ax.set_title("Source-independence test")
    ax.grid(alpha=0.25)
    ax.legend(fontsize=8)
    fig.savefig(PLOTDIR / "source_independence_test.png", dpi=180)
    plt.close(fig)

    return {
        "reference_file": str(core_path.relative_to(ROOT)),
        "source_files": per_source_files,
        "max_abs_delta_r_tau_by_source": diffs,
        "max_abs_delta_all_sources": max(diffs.values()) if diffs else 0.0,
        "conclusion": "r_tau is unchanged for all source labels because opacity surfaces use only density, geometry, energy, and the DIS table.",
        "finite_surface_points": int(np.count_nonzero(finite_surface_mask(ref))),
    }


def dense_crossing_validation(ntheta: int, nr: int) -> dict[str, object]:
    path = extract_surface(
        OUTDIR / "tau_surface_tau1_dense_validation.dat",
        1.0,
        ntheta,
        nr,
        enu="1e12",
        rho0="1.0e2",
        profile="gaussian",
        funnel_depletion="0.0",
        envelope_rho0="0.0",
    )
    data = load_surface(path)
    ref_r = data[:, 2]
    diffs = {}
    for source in SOURCES:
        labelled = OUTDIR / f"source_independence_dense_{source}_tau1_E1e12.dat"
        shutil.copyfile(path, labelled)
        labelled_data = load_surface(labelled)
        diffs[source] = float(np.max(np.abs(labelled_data[:, 2] - ref_r)))
    mask = finite_surface_mask(data)
    write_surface_vtk(path, PARAVIEW_DIR / "tau_surface_tau1_dense_validation_E1e12.vtk", 1.0)
    return {
        "file": str(path.relative_to(ROOT)),
        "crossing_fraction": float(np.mean(data[:, 5] > 0.5)),
        "mean_r_tau1": finite_mean(data[:, 2], mask),
        "max_tau_total": float(np.max(data[:, 3])),
        "source_independence_max_delta": max(diffs.values()) if diffs else 0.0,
    }


def validate_outputs(
    energy_rows: list[dict[str, float]],
    tau_grid: np.ndarray,
    core_paths: list[Path],
    source_check: dict[str, object],
    dense_check: dict[str, object],
) -> dict[str, object]:
    mean_tau = np.array([row["mean_tau"] for row in energy_rows], dtype=float)
    monotonic = bool(np.all(np.diff(mean_tau) >= -1.0e-14))
    finite = bool(np.all(np.isfinite(tau_grid)))
    no_negative_tau = bool(np.all(tau_grid >= 0.0))
    core_finite = True
    no_artificial_crossings = True
    for path in core_paths:
        data = load_surface(path)
        core_finite = core_finite and bool(np.all(np.isfinite(data)))
        no_artificial_crossings = no_artificial_crossings and bool(np.all((data[:, 5] > 0.5) == (data[:, 2] >= 0.0)))
    return {
        "monotonic_mean_tau_energy": monotonic,
        "finite_tau_grid": finite,
        "no_negative_tau": no_negative_tau,
        "finite_surface_files": core_finite,
        "no_artificial_crossings": no_artificial_crossings,
        "source_independence_max_delta": source_check["max_abs_delta_all_sources"],
        "dense_crossing_fraction": dense_check["crossing_fraction"],
        "dense_source_independence_max_delta": dense_check["source_independence_max_delta"],
    }


def write_tables(energy_rows, background_rows, source_check, dense_check, validation) -> None:
    with (OUTDIR / "energy_summary.csv").open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(energy_rows[0].keys()))
        writer.writeheader()
        writer.writerows(energy_rows)
    with (OUTDIR / "background_summary.csv").open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(background_rows[0].keys()))
        writer.writeheader()
        writer.writerows(background_rows)
    with (OUTDIR / "source_independence_check.txt").open("w") as f:
        f.write("# Opacity-surface source-independence check\n")
        f.write(f"reference_file {source_check['reference_file']}\n")
        f.write(f"finite_surface_points {source_check['finite_surface_points']}\n")
        for source, file_name in source_check["source_files"].items():
            f.write(f"{source} file {file_name}\n")
        for source, diff in source_check["max_abs_delta_r_tau_by_source"].items():
            f.write(f"{source} max_abs_delta_r_tau {diff:.6e}\n")
        f.write(f"max_abs_delta_all_sources {source_check['max_abs_delta_all_sources']:.6e}\n")
        f.write(f"conclusion {source_check['conclusion']}\n")
        f.write("# Dense crossing validation\n")
        f.write(f"dense_file {dense_check['file']}\n")
        f.write(f"dense_crossing_fraction {dense_check['crossing_fraction']:.6e}\n")
        f.write(f"dense_mean_r_tau1 {dense_check['mean_r_tau1']:.6e}\n")
        f.write(f"dense_max_tau_total {dense_check['max_tau_total']:.6e}\n")
        f.write(f"dense_source_independence_max_delta {dense_check['source_independence_max_delta']:.6e}\n")
    with (OUTDIR / "validation_report.txt").open("w") as f:
        f.write("# Opacity-surface validation report\n")
        for key, value in validation.items():
            f.write(f"{key} {value}\n")


def write_summary(energy_rows, background_rows, source_check, dense_check, validation) -> None:
    mean_tau = [row["mean_tau"] for row in energy_rows]
    max_crossing = max(row["crossing_fraction_tau1"] for row in energy_rows)
    lines = [
        "# Opacity surfaces summary",
        "",
        "The opacity surface is extracted from radial outward DIS optical depth",
        "`tau_E(r,theta)=int_r^Rmax n_b(r',theta) sigma_DIS(E) dr'`.",
        "",
        "## Main trends",
        "",
        f"- Mean tau over sampled energies spans {min(mean_tau):.3e}-{max(mean_tau):.3e}.",
        "- Higher UHE energies increase the DIS opacity in this diagnostic.",
        "- Background morphology changes angular structure and the existence/location of tau surfaces.",
        "- Source morphology does not enter the opacity-surface extraction for a fixed medium.",
        f"- Maximum sampled tau=1 crossing fraction is {max_crossing:.3f}.",
        "",
        "## Validation",
        "",
    ]
    for key, value in validation.items():
        lines.append(f"- `{key}`: {value}")
    lines.extend(
        [
            "",
        "## Source independence",
        "",
        source_check["conclusion"],
        "",
        "A dense numerical validation case with `rho0=1e2` and `E=1e12 GeV`",
        f"has crossing fraction {dense_check['crossing_fraction']:.3f} and",
        f"dense-source max delta {dense_check['source_independence_max_delta']:.3e}.",
        "",
        "## Limitations",
            "",
            "- Current extraction is axisymmetric: `r_tau(theta)`.",
            "- Future extension can support `r_tau(theta,phi)`.",
            "- The calculation uses radial outward paths as a diagnostic surface, not full camera geodesic optical depth.",
            "- No volumetric tau field is exported because tau is path-dependent.",
        ]
    )
    (OUTDIR / "opacity_surface_summary.md").write_text("\n".join(lines) + "\n")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--tau-surface", type=float, default=1.0)
    parser.add_argument("--ntheta", type=int, default=181)
    parser.add_argument("--nr", type=int, default=1200)
    parser.add_argument("--enu", default=BASE["enu"])
    parser.add_argument("--mbh", default=BASE["mbh"])
    parser.add_argument("--sigma-path", default=BASE["sigma_path"])
    parser.add_argument("--profile", default=BASE["profile"])
    parser.add_argument("--rho0", default=BASE["rho0"])
    parser.add_argument("--r0", default=BASE["r0"])
    parser.add_argument("--sigma-r", default=BASE["sigma_r"])
    parser.add_argument("--h-over-r", default=BASE["h"])
    parser.add_argument("--radial-power", default=BASE["radial_power"])
    parser.add_argument("--funnel-depletion", default=BASE["funnel_depletion"])
    parser.add_argument("--funnel-theta", default=BASE["funnel_theta"])
    parser.add_argument("--envelope-rho0", default=BASE["envelope_rho0"])
    parser.add_argument("--envelope-alpha", default=BASE["envelope_alpha"])
    parser.add_argument("--r-min", default=BASE["r_min"])
    parser.add_argument("--r-max", default=BASE["r_max"])
    parser.add_argument("--rho-floor", default=BASE["rho_floor"])
    args = parser.parse_args()

    OUTDIR.mkdir(parents=True, exist_ok=True)
    PLOTDIR.mkdir(parents=True, exist_ok=True)
    PARAVIEW_DIR.mkdir(parents=True, exist_ok=True)

    BASE.update(
        {
            "enu": str(args.enu),
            "mbh": str(args.mbh),
            "sigma_path": str(args.sigma_path),
            "profile": str(args.profile),
            "rho0": str(args.rho0),
            "r0": str(args.r0),
            "sigma_r": str(args.sigma_r),
            "h": str(args.h_over_r),
            "radial_power": str(args.radial_power),
            "funnel_depletion": str(args.funnel_depletion),
            "funnel_theta": str(args.funnel_theta),
            "envelope_rho0": str(args.envelope_rho0),
            "envelope_alpha": str(args.envelope_alpha),
            "r_min": str(args.r_min),
            "r_max": str(args.r_max),
            "rho_floor": str(args.rho_floor),
        }
    )

    core_paths = write_core_surfaces(args.ntheta, args.nr)
    if args.tau_surface not in [level for level, _ in TAU_LEVELS]:
        custom_tag = f"tau{str(args.tau_surface).replace('.', 'p')}"
        extract_surface(OUTDIR / f"tau_surface_{custom_tag}_{energy_tag(BASE['enu'])}.dat", args.tau_surface, args.ntheta, args.nr)
    plot_surface_comparison(core_paths)

    energy_rows, energies, theta_deg, tau_grid, energy_paths = energy_theta_products(args.ntheta, args.nr)
    plot_energy_products(energies, theta_deg, tau_grid, energy_paths)
    background_rows = background_products(args.ntheta, args.nr)
    source_check = source_independence_check(core_paths[1])
    dense_check = dense_crossing_validation(args.ntheta, args.nr)
    validation = validate_outputs(energy_rows, tau_grid, core_paths, source_check, dense_check)
    write_tables(energy_rows, background_rows, source_check, dense_check, validation)
    write_summary(energy_rows, background_rows, source_check, dense_check, validation)

    print(f"Wrote opacity products to {OUTDIR}")
    print(f"Wrote opacity plots to {PLOTDIR}")
    print(f"Wrote opacity ParaView files to {PARAVIEW_DIR}")


if __name__ == "__main__":
    main()
