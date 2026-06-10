"""Run moderate robustness scans for Point 3."""

from __future__ import annotations

import csv
import math
import subprocess
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
OUTDIR = ROOT / "output" / "scans"
PLOTDIR = ROOT / "plots" / "scans"


BASE = {
    "enu": "1e5",
    "aspin": "0.0001",
    "mbh": "3.0",
    "rho0": "1.0e-2",
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
    "density_profile": "gaussian",
    "radial_power": "2.0",
    "funnel_depletion": "0.0",
    "funnel_theta": "20.0",
    "envelope_rho0": "0.0",
    "envelope_alpha": "2.5",
    "r_min": "4.0",
    "r_max": "60.0",
    "rho_floor": "1.0e-99",
    "observer_distance": "60.0",
    "use_f3": "1",
    "source_model": "inner_ring",
    "source_funnel_theta": "20.0",
    "source_q": "1.0",
    "source_s": "2.0",
    "source_grad_dr": "0.1",
    "source_grad_dtheta": "1.0",
    "source_rho_ref": "1.0e-2",
    "source_cutoff_min": "0.0",
    "source_cutoff_max": "1.0e2",
}


def run(cmd: list[str], env: dict[str, str] | None = None) -> str:
    result = subprocess.run(
        cmd,
        cwd=ROOT,
        env=env,
        check=True,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
    )
    return result.stdout


def ensure_build() -> None:
    run(["make", "build", "NTHREADS=2"])


def cache_path(label: str) -> Path:
    return ROOT / "output" / "rays" / f"scan_{label}.bin"


def make_cache(label: str, spin: str, theta: str, nx: int = 16, ny: int = 16) -> Path:
    path = cache_path(label)
    if path.exists():
        return path
    run(
        [
            "./compute_kerr_geodesics",
            spin,
            BASE["observer_distance"],
            theta,
            "25.0",
            str(nx),
            str(ny),
            "60.0",
            "0.02",
            str(path.relative_to(ROOT)),
        ],
        env={"OMP_NUM_THREADS": "2"},
    )
    return path


def image_args(params: dict[str, str], cache: Path) -> list[str]:
    return [
        "./compute_kerr_image_from_cache",
        params["enu"],
        params["aspin"],
        params["mbh"],
        params["rho0"],
        params["r0"],
        params["sigma"],
        params["h"],
        params["source_r"],
        params["source_sigma"],
        params["source_theta"],
        params["source_powerlaw"],
        params["source_emax"],
        params["source_norm"],
        params["mev_enu"],
        params["mev_norm"],
        params["cam_theta"],
        params["sigma_path"],
        params["density_profile"],
        params["radial_power"],
        params["funnel_depletion"],
        params["funnel_theta"],
        params["envelope_rho0"],
        params["envelope_alpha"],
        params["r_min"],
        params["r_max"],
        params["rho_floor"],
        params["observer_distance"],
        params["use_f3"],
        str(cache.relative_to(ROOT)),
        params["source_model"],
        params["source_funnel_theta"],
        params["source_q"],
        params["source_s"],
        params["source_grad_dr"],
        params["source_grad_dtheta"],
        params["source_rho_ref"],
        params["source_cutoff_min"],
        params["source_cutoff_max"],
    ]


def run_case(scan: str, case: str, cache: Path, **updates: str) -> dict[str, object]:
    params = dict(BASE)
    params.update({k: str(v) for k, v in updates.items()})
    stdout = run(image_args(params, cache), env={"OMP_NUM_THREADS": "2"})
    saved = None
    for line in stdout.splitlines():
        if line.startswith("Saved: ") and line.endswith(".dat"):
            saved = ROOT / line.split("Saved: ", 1)[1]
    if saved is None:
        raise RuntimeError(f"No output file found for {scan}:{case}")
    metrics = image_metrics(saved)
    row = {
        "scan": scan,
        "case": case,
        "output_file": str(saved.relative_to(ROOT)),
        **params,
        **metrics,
    }
    return row


def image_metrics(path: Path) -> dict[str, float | int]:
    data = np.loadtxt(path, comments="#")
    if data.ndim == 1:
        data = data.reshape(1, -1)
    tau = data[:, 4]
    psurv = data[:, 5]
    intensity = data[:, 6]
    captured = data[:, 7]
    finite = np.isfinite(tau) & np.isfinite(psurv) & np.isfinite(intensity)
    return {
        "mean_tau": float(np.mean(tau[finite])),
        "max_tau": float(np.max(tau[finite])),
        "mean_P_surv": float(np.mean(psurv[finite])),
        "total_intensity": float(np.sum(intensity[finite])),
        "valid_rays": int(np.count_nonzero(finite & (captured == 0))),
    }


def run_scans() -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    base_cache = make_cache("base_a0001_th80", BASE["aspin"], BASE["cam_theta"])

    for model in ["inner_ring", "jet_base", "funnel_wall", "density_weighted", "shock_layer"]:
        rows.append(run_case("source_model", model, base_cache, source_model=model))

    density_cases = {
        "gaussian": {},
        "powerlaw": {"density_profile": "powerlaw"},
        "powerlaw_funnel": {
            "density_profile": "powerlaw_funnel",
            "funnel_depletion": "1.0",
            "funnel_theta": "20.0",
        },
        "powerlaw_funnel_envelope": {
            "density_profile": "powerlaw_funnel_envelope",
            "funnel_depletion": "1.0",
            "funnel_theta": "20.0",
            "envelope_rho0": "1.0e-4",
        },
    }
    for name, updates in density_cases.items():
        rows.append(run_case("density_background", name, base_cache, source_model="funnel_wall", **updates))

    for model in ["inner_ring", "jet_base", "funnel_wall", "density_weighted", "shock_layer"]:
        previous_tau = None
        for enu in ["1e4", "1e5", "1e6"]:
            row = run_case("energy", f"{model}_E{enu}", base_cache, source_model=model, enu=enu)
            if previous_tau is not None and row["mean_tau"] < previous_tau:
                row["anomaly"] = "mean_tau_decreased"
            else:
                row["anomaly"] = ""
            previous_tau = row["mean_tau"]
            rows.append(row)

    for theta in ["10.0", "30.0", "50.0", "70.0", "85.0"]:
        c = make_cache(f"inclination_th{theta.replace('.', 'p')}", BASE["aspin"], theta)
        rows.append(run_case("inclination", f"theta_{theta}", c, cam_theta=theta, source_model="funnel_wall"))

    for spin in ["0.0", "0.5", "0.9", "0.99"]:
        c = make_cache(f"spin_a{spin.replace('.', 'p')}", spin, BASE["cam_theta"])
        rows.append(run_case("spin", f"a_{spin}", c, aspin=spin, source_model="funnel_wall"))

    for q in ["0.0", "0.5", "1.0", "2.0"]:
        rows.append(
            run_case(
                "parameter_stability",
                f"density_weighted_q{q}",
                base_cache,
                source_model="density_weighted",
                source_q=q,
            )
        )
    for s in ["0.0", "1.0", "2.0", "3.0"]:
        rows.append(
            run_case(
                "parameter_stability",
                f"density_weighted_s{s}",
                base_cache,
                source_model="density_weighted",
                source_s=s,
            )
        )
    for wall in ["10.0", "20.0", "30.0"]:
        rows.append(
            run_case(
                "parameter_stability",
                f"funnel_wall_theta{wall}",
                base_cache,
                source_model="funnel_wall",
                source_funnel_theta=wall,
            )
        )
    for width in ["8.0", "15.0", "25.0"]:
        rows.append(
            run_case(
                "parameter_stability",
                f"funnel_wall_width{width}",
                base_cache,
                source_model="funnel_wall",
                source_theta=width,
            )
        )

    return rows


def write_csv(rows: list[dict[str, object]]) -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    keys = sorted({key for row in rows for key in row})
    with (OUTDIR / "scan_summary.csv").open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=keys)
        writer.writeheader()
        writer.writerows(rows)


def rows_for(rows, scan):
    return [r for r in rows if r["scan"] == scan]


def bar_plot(rows, title, filename):
    labels = [r["case"] for r in rows]
    x = np.arange(len(labels))
    fig, axes = plt.subplots(2, 2, figsize=(11, 7), constrained_layout=True)
    for ax, metric in zip(
        axes.flat,
        ["mean_tau", "max_tau", "mean_P_surv", "total_intensity"],
    ):
        ax.bar(x, [float(r[metric]) for r in rows])
        ax.set_title(metric)
        ax.set_xticks(x)
        ax.set_xticklabels(labels, rotation=35, ha="right")
        ax.grid(axis="y", alpha=0.25)
    fig.suptitle(title)
    fig.savefig(PLOTDIR / filename, dpi=170)
    plt.close(fig)


def line_plot(rows, group_key, x_getter, metrics, title, filename):
    fig, axes = plt.subplots(1, len(metrics), figsize=(5.3 * len(metrics), 4.2), constrained_layout=True)
    if len(metrics) == 1:
        axes = [axes]
    groups = {}
    for r in rows:
        groups.setdefault(r[group_key], []).append(r)
    for label, group in groups.items():
        group = sorted(group, key=x_getter)
        xs = [x_getter(r) for r in group]
        for ax, metric in zip(axes, metrics):
            ax.plot(xs, [float(r[metric]) for r in group], marker="o", label=label)
    for ax, metric in zip(axes, metrics):
        ax.set_title(metric)
        ax.grid(alpha=0.25)
        ax.legend(fontsize=7)
    fig.suptitle(title)
    fig.savefig(PLOTDIR / filename, dpi=170)
    plt.close(fig)


def write_plots(rows: list[dict[str, object]]) -> None:
    PLOTDIR.mkdir(parents=True, exist_ok=True)
    bar_plot(rows_for(rows, "source_model"), "Source model dependence", "source_model_comparison.png")
    bar_plot(rows_for(rows, "density_background"), "Density background dependence", "density_background_comparison.png")

    energy_rows = rows_for(rows, "energy")
    for r in energy_rows:
        r["energy_group"] = str(r["case"]).split("_E", 1)[0]
    line_plot(
        energy_rows,
        "energy_group",
        lambda r: float(r["enu"]),
        ["mean_tau", "mean_P_surv"],
        "Energy dependence",
        "energy_dependence.png",
    )

    line_plot(
        rows_for(rows, "inclination"),
        "scan",
        lambda r: float(r["cam_theta"]),
        ["mean_tau", "total_intensity"],
        "Observer inclination dependence",
        "inclination_dependence.png",
    )

    line_plot(
        rows_for(rows, "spin"),
        "scan",
        lambda r: float(r["aspin"]),
        ["mean_tau", "total_intensity"],
        "Black-hole spin dependence",
        "spin_dependence.png",
    )

    bar_plot(rows_for(rows, "parameter_stability"), "Parameter stability", "parameter_stability.png")


def write_markdown(rows: list[dict[str, object]]) -> None:
    energy_anomalies = [r for r in rows if r.get("anomaly")]
    lines = [
        "# Robustness scan summary",
        "",
        "Moderate-resolution scans were run for Point 3. These are validation-scale",
        "parameter studies, not production-resolution figures.",
        "",
        "## Key observations",
        "",
    ]

    source = rows_for(rows, "source_model")
    if source:
        tau_values = [float(r["mean_tau"]) for r in source]
        lines.append(
            f"- Source morphology changes total intensity, but mean tau stays in "
            f"the range {min(tau_values):.3e}-{max(tau_values):.3e} for the fixed Gaussian background."
        )

    density = rows_for(rows, "density_background")
    if density:
        tau_values = [float(r["mean_tau"]) for r in density]
        lines.append(
            f"- Density morphology changes mean tau over {min(tau_values):.3e}-{max(tau_values):.3e}."
        )

    if energy_anomalies:
        lines.append(f"- Energy scan anomalies flagged: {len(energy_anomalies)} cases.")
    else:
        lines.append("- No monotonicity anomalies were flagged in mean tau over the sampled UHE energies.")

    spin = rows_for(rows, "spin")
    if spin:
        intensities = [float(r["total_intensity"]) for r in spin]
        lines.append(
            f"- Spin modifies image intensity over {min(intensities):.3e}-{max(intensities):.3e} "
            "in this validation grid."
        )

    lines.extend(
        [
            "",
            "## Outputs",
            "",
            "- `output/scans/scan_summary.csv`",
            "- `plots/scans/source_model_comparison.png`",
            "- `plots/scans/density_background_comparison.png`",
            "- `plots/scans/energy_dependence.png`",
            "- `plots/scans/inclination_dependence.png`",
            "- `plots/scans/spin_dependence.png`",
            "- `plots/scans/parameter_stability.png`",
            "",
            "## Interpretation guardrails",
            "",
            "These scans support robustness of qualitative trends at moderate resolution.",
            "They should not be presented as production-resolution convergence tests.",
        ]
    )
    (OUTDIR / "scan_summary.md").write_text("\n".join(lines) + "\n")


def main() -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    PLOTDIR.mkdir(parents=True, exist_ok=True)
    ensure_build()
    rows = run_scans()
    write_csv(rows)
    write_markdown(rows)
    write_plots(rows)
    print(f"Wrote {OUTDIR / 'scan_summary.csv'}")
    print(f"Wrote {OUTDIR / 'scan_summary.md'}")
    print(f"Wrote plots to {PLOTDIR}")


if __name__ == "__main__":
    main()
