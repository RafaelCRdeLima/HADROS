#!/usr/bin/env python3
"""Compare named UHE/MeV/collapsar semi-analytic regimes."""

from __future__ import annotations

import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
OUT_DIR = ROOT / "output" / "collapsar_ndaf_like"
PLOT_DIR = ROOT / "plots" / "collapsar_ndaf_like"
SIGMA_PATH = ROOT / "data" / "sigma" / "sigma_nuN_CC_GBW.dat"

PI = math.pi
M_U_G = 1.66053906660e-24
MSUN_G = 1.98847e33
MEV_TO_ERG = 1.602176634e-6


def rg_cm(m_msun=3.0):
    return 6.67430e-8 * m_msun * 1.98847e33 / (2.99792458e10**2)


PRESETS = {
    "fiducial_mev_density": {
        "profile": "gaussian",
        "rho0": 1.0e10,
        "r0": 10.0,
        "sigma": 5.0,
        "H": 0.25,
        "rmin": 4.0,
        "rmax": 60.0,
        "T0": 6.0,
        "Tfloor": 0.1,
        "Ye_torus": 0.25,
        "Ye_funnel": 0.55,
    },
    "collapsar_ndaf_like": {
        "profile": "collapsar_ndaf_like",
        "rho0": 3.0e10,
        "r0": 12.0,
        "sigma": 8.0,
        "H": 0.45,
        "rmin": 3.0,
        "rmax": 80.0,
        "T0": 18.0,
        "Tfloor": 0.2,
        "Ye_torus": 0.20,
        "Ye_funnel": 0.42,
    },
}

UHE_COMPARISON_PRESETS = {
    "fiducial_uhe_default": {
        **PRESETS["fiducial_mev_density"],
        "rho0": 1.0e-2,
    },
    "fiducial_mev_density": PRESETS["fiducial_mev_density"],
    "collapsar_ndaf_like": PRESETS["collapsar_ndaf_like"],
}


def make_grid(nr=220, nth=160, rmax=90.0):
    r_edges = np.linspace(2.5, rmax, nr + 1)
    th_edges = np.linspace(1e-4, PI - 1e-4, nth + 1)
    r = 0.5 * (r_edges[:-1] + r_edges[1:])
    th = 0.5 * (th_edges[:-1] + th_edges[1:])
    dr = np.diff(r_edges) * rg_cm()
    dth = np.diff(th_edges)
    rr, tt = np.meshgrid(r, th, indexing="ij")
    drr, dtt = np.meshgrid(dr, dth, indexing="ij")
    rcm = rr * rg_cm()
    dV = 2 * PI * rcm * rcm * np.sin(tt) * drr * dtt
    return rr, tt, dV


def density(p, r, th):
    delta = th - 0.5 * PI
    if p["profile"] == "gaussian":
        g = np.exp(-((r - p["r0"]) / p["sigma"]) ** 2) * np.exp(-(delta / p["H"]) ** 2)
        g = np.where((r <= 4.0) | (r >= 18.0) | (np.abs(delta) >= 0.45), 0.0, g)
        return np.maximum(p["rho0"] * g, 1e-20)
    mu = np.cos(th)
    vertical = np.exp(-(mu / p["H"]) ** 2)
    inner = 1.0 / (1.0 + np.exp(-(r - p["rmin"]) / (0.35 * p["sigma"])))
    outer = 1.0 / (1.0 + np.exp((r - p["rmax"]) / (0.60 * p["sigma"])))
    core = 0.45 * np.exp(-((r - 0.75 * p["r0"]) / (0.75 * p["sigma"])) ** 2)
    radial = np.maximum(r / p["r0"], 1e-300) ** -2.0
    tail = 0.08 * np.exp(-((r - p["r0"]) / (1.8 * p["sigma"])) ** 2)
    shape = (radial + core + tail) * vertical * inner * outer
    return np.maximum(p["rho0"] * shape, 1e-20)


def temperature(p, r, th, rho):
    shape = np.clip(rho / max(p["rho0"], 1e-300), 0.0, 1.0)
    equatorial = np.exp(-((th - 0.5 * PI) / 0.45) ** 2)
    if p["profile"] == "gaussian":
        return np.maximum(p["Tfloor"] + p["T0"] * shape**0.2, p["Tfloor"])
    inner_hot = np.exp(-((r - p["rmin"]) / (1.15 * p["r0"])) ** 2)
    radial_cooling = (1 + np.maximum(r - p["rmin"], 0) / p["r0"]) ** -0.75
    T = (
        p["Tfloor"]
        + p["T0"]
        * (0.30 + 0.70 * shape**0.18)
        * (0.45 + 0.55 * equatorial)
        * (0.55 + 0.45 * inner_hot)
        * radial_cooling
    )
    return np.maximum(T, p["Tfloor"])


def ye_profile(p, r, th, rho):
    shape = np.clip(rho / max(p["rho0"], 1e-300), 0.0, 1.0)
    polar = np.exp(-(th / 0.35) ** 2) + np.exp(-((PI - th) / 0.35) ** 2)
    polar = np.clip(polar, 0, 1)
    if p["profile"] == "gaussian":
        ye = p["Ye_torus"] + 0.05 * np.exp(-((r - p["r0"]) / 4.0) ** 2)
        return np.clip(ye, 0.01, 0.60)
    radial_neutron = np.exp(-((r - p["r0"]) / (1.2 * p["r0"])) ** 2)
    torus_ye = p["Ye_torus"] - 0.05 * shape**0.25 * radial_neutron
    ye = torus_ye * (1 - polar) + p["Ye_funnel"] * polar
    return np.clip(ye, 0.05, 0.60)


def thermal_shape(E, T):
    x = E / np.maximum(T, 1e-300)
    return np.where(x > 120, 0.0, E * E / (np.exp(np.minimum(x, 120)) + 1) / np.maximum(T**3, 1e-300))


def luminosity_proxy(rho, T, Ye, dV):
    e_edges = np.logspace(0, np.log10(80), 32)
    total = np.zeros_like(rho)
    for e1, e2 in zip(e_edges[:-1], e_edges[1:]):
        e = math.sqrt(e1 * e2)
        de = e2 - e1
        rho10 = rho / 1e10
        t10 = T / 10
        spec = thermal_shape(e, T)
        urca = 1e30 * rho10 * t10**6 * spec
        pair = 3e28 * t10**9 * 2.4 * spec
        brems = 1e27 * rho10**2 * t10**5 * 2.6 * spec
        total += (urca + pair + brems) * e * MEV_TO_ERG * de
    return float(np.sum(total * dV)), total * dV


def stats(vals, mask, weights=None):
    v = vals[mask]
    if weights is None:
        return np.min(v), np.max(v), np.mean(v), np.median(v)
    w = weights[mask]
    good = w > 0
    v, w = v[good], w[good]
    order = np.argsort(v)
    c = np.cumsum(w[order]) / np.sum(w)
    return np.min(v), np.max(v), np.sum(v * w) / np.sum(w), np.interp(0.5, c, v[order])


def log10_sigma_cm2(E):
    data = np.loadtxt(SIGMA_PATH, comments="#")
    return float(np.interp(np.log10(E), np.log10(data[:, 0]), np.log10(data[:, 2])))


def radial_opacity_diagnostic(p, E_GeV, tau_target=1.0):
    theta = np.linspace(1e-4, PI - 1e-4, 181)
    r = np.linspace(p["rmin"], p["rmax"], 1200)
    sig = 10 ** log10_sigma_cm2(E_GeV)
    out = np.full_like(theta, np.nan)
    tau_inner = np.zeros_like(theta)
    mass_column = np.zeros_like(theta)
    baryon_column = np.zeros_like(theta)
    for i, th in enumerate(theta):
        rr = r
        rho = density(p, rr, np.full_like(rr, th))
        ds = np.gradient(rr) * rg_cm()
        alpha = rho / M_U_G * sig
        tau = np.cumsum((alpha * ds)[::-1])[::-1]
        mass_column[i] = float(np.sum(rho * ds))
        baryon_column[i] = float(np.sum(rho / M_U_G * ds))
        tau_inner[i] = tau[0]
        if np.max(tau) < tau_target:
            continue
        idx = np.where(tau >= tau_target)[0][-1]
        if idx >= len(r) - 1:
            out[i] = r[idx]
        else:
            t1, t2 = tau[idx], tau[idx + 1]
            r1, r2 = r[idx], r[idx + 1]
            out[i] = r2 + (tau_target - t2) / max(t1 - t2, 1e-300) * (r1 - r2)
    return {
        "theta": theta,
        "surface_r": out,
        "tau_dimensionless": tau_inner,
        "mass_column_gcm2": mass_column,
        "baryon_column_cm2": baryon_column,
        "sigma_cm2": sig,
    }


def run_case(name, p):
    rr, tt, dV = make_grid(rmax=max(90.0, p["rmax"]))
    rho = density(p, rr, tt)
    T = temperature(p, rr, tt, rho)
    Ye = ye_profile(p, rr, tt, rho)
    mask = rho > max(1e-3 * p["rho0"], 1e-18)
    mass = float(np.sum(rho[mask] * dV[mask]))
    L, Lcell = luminosity_proxy(rho, T, Ye, dV)
    result = {"mass": mass, "L": L}
    for label, vals in [("rho", rho), ("T", T), ("Ye", Ye)]:
        result[label] = stats(vals, mask)
        result[label + "_w"] = stats(vals, mask, Lcell)
    return result


def run_uhe_case(name, p, E_GeV=1e11):
    diag = radial_opacity_diagnostic({**p, "rmax": p["rmax"]}, E_GeV, 1.0)
    tau = diag["tau_dimensionless"]
    result = {
        "theta": diag["theta"],
        "r_uhe": diag["surface_r"],
        "sigma_cm2": diag["sigma_cm2"],
        "mean_tau_dimensionless": float(np.mean(tau)),
        "max_tau_dimensionless": float(np.max(tau)),
        "mean_Psurv": float(np.mean(np.exp(-np.minimum(tau, 700)))),
        "mean_baryon_column_cm2": float(np.mean(diag["baryon_column_cm2"])),
        "max_baryon_column_cm2": float(np.max(diag["baryon_column_cm2"])),
        "mean_mass_column_gcm2": float(np.mean(diag["mass_column_gcm2"])),
        "max_mass_column_gcm2": float(np.max(diag["mass_column_gcm2"])),
    }
    return result


def run_mev_surface_case(p):
    diag = radial_opacity_diagnostic({**p, "rmax": p["rmax"]}, 10.0, 2.0 / 3.0)
    return diag["theta"], diag["surface_r"]


def write_case_md(name, result, path):
    lines = [f"# {name} Statistics", ""]
    lines.append(f"- M_torus: `{result['mass']:.8e}` g = `{result['mass']/MSUN_G:.8e}` Msun")
    lines.append(f"- L_proxy: `{result['L']:.8e}` erg/s-proxy")
    for key in ["rho", "T", "Ye", "rho_w", "T_w", "Ye_w"]:
        mn, mx, mean, med = result[key]
        lines += [f"## {key}", f"- min: `{mn:.8e}`", f"- max: `{mx:.8e}`", f"- mean: `{mean:.8e}`", f"- median: `{med:.8e}`", ""]
    lines += [
        "## UHE radial opacity diagnostic",
        "",
        "UHE opacity diagnostics are written separately in `uhe_opacity_unit_audit.md`.",
        "This avoids mixing `fiducial_mev_density` with `fiducial_uhe_default`.",
        "",
    ]
    path.write_text("\n".join(lines))


def write_uhe_audit(uhe_results):
    lines = [
        "# UHE Opacity Unit Audit",
        "",
        "The audited quantity is true dimensionless optical depth:",
        "",
        "```text",
        "tau(E, theta) = integral n_b(r,theta) sigma_nuN(E) ds",
        "```",
        "",
        "with:",
        "",
        "```text",
        "n_b = rho / m_u              [cm^-3]",
        "sigma_nuN(E)                [cm^2]",
        "ds                          [cm]",
        "tau                         dimensionless",
        "baryon_column = integral n_b ds     [cm^-2]",
        "mass_column = integral rho ds       [g cm^-2]",
        "```",
        "",
        "The earlier `mean_tau` label in the collapsar comparison was dimensionless tau, not a column density. However, it was misleading because the entry called only `fiducial` used `fiducial_mev_density` (`rho0=1e10 g/cm3`), not `fiducial_uhe_default` (`rho0=1e-2 g/cm3`).",
        "",
        "| Case | sigma_cm2 | mean_tau_dimensionless | max_tau_dimensionless | mean_baryon_column_cm2 | mean_mass_column_gcm2 | mean_Psurv |",
        "|---|---:|---:|---:|---:|---:|---:|",
    ]
    for name, res in uhe_results.items():
        lines.append(
            f"| {name} | {res['sigma_cm2']:.8e} | {res['mean_tau_dimensionless']:.8e} | "
            f"{res['max_tau_dimensionless']:.8e} | {res['mean_baryon_column_cm2']:.8e} | "
            f"{res['mean_mass_column_gcm2']:.8e} | {res['mean_Psurv']:.8e} |"
        )
    lines += [
        "",
        "Paper-use status: do not use the previous `uhe_opacity_comparison.png`; use the regenerated plot with explicit `tau_dimensionless` labeling.",
        "",
    ]
    (OUT_DIR / "uhe_opacity_unit_audit.md").write_text("\n".join(lines))


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    PLOT_DIR.mkdir(parents=True, exist_ok=True)
    results = {name: run_case(name, p) for name, p in PRESETS.items()}
    uhe_results = {name: run_uhe_case(name, p) for name, p in UHE_COMPARISON_PRESETS.items()}
    write_uhe_audit(uhe_results)
    write_case_md("collapsar_ndaf_like", results["collapsar_ndaf_like"], OUT_DIR / "collapsar_statistics.md")
    write_case_md("fiducial_mev_density", results["fiducial_mev_density"], OUT_DIR / "fiducial_mev_density_statistics.md")
    (OUT_DIR / "collapsar_mass_estimate.md").write_text(
        f"# Collapsar NDAF-like Mass Estimate\n\n"
        f"- M_torus: `{results['collapsar_ndaf_like']['mass']:.8e}` g = "
        f"`{results['collapsar_ndaf_like']['mass']/MSUN_G:.8e}` Msun\n"
        f"- fiducial_mev_density_M_torus: `{results['fiducial_mev_density']['mass']/MSUN_G:.8e}` Msun\n"
        f"- mass_ratio_vs_fiducial_mev_density: `{results['collapsar_ndaf_like']['mass']/results['fiducial_mev_density']['mass']:.8e}`\n"
    )
    ratio = results["collapsar_ndaf_like"]["L"] / results["fiducial_mev_density"]["L"]
    (OUT_DIR / "comparison_summary.md").write_text(
        "\n".join([
            "# Collapsar NDAF-like Preset Comparison",
            "",
            f"- fiducial_mev_density_L_proxy: `{results['fiducial_mev_density']['L']:.8e}`",
            f"- collapsar_ndaf_like_L_proxy: `{results['collapsar_ndaf_like']['L']:.8e}`",
            f"- luminosity_ratio: `{ratio:.8e}`",
            f"- fiducial_mev_density_mass_Msun: `{results['fiducial_mev_density']['mass']/MSUN_G:.8e}`",
            f"- collapsar_mass_Msun: `{results['collapsar_ndaf_like']['mass']/MSUN_G:.8e}`",
            f"- collapsar_T_weighted_median_MeV: `{results['collapsar_ndaf_like']['T_w'][3]:.8e}`",
            f"- collapsar_rho_weighted_median_gcm3: `{results['collapsar_ndaf_like']['rho_w'][3]:.8e}`",
            "",
            "The preset changes thermodynamic regime without changing emissivity normalization.",
            "",
            "UHE opacity units and the fiducial-density distinction are audited in `uhe_opacity_unit_audit.md`.",
            "",
        ])
    )
    plt.figure(figsize=(5.8, 4.2))
    plt.bar(["fiducial_mev\ndensity", "collapsar\nndaf-like"], [results["fiducial_mev_density"]["L"], results["collapsar_ndaf_like"]["L"]])
    plt.yscale("log")
    plt.ylabel(r"$L_{\nu,\rm proxy}$ [erg s$^{-1}$ proxy]")
    plt.title(f"MeV luminosity proxy comparison\nratio={ratio:.2e}")
    plt.tight_layout()
    plt.savefig(PLOT_DIR / "mev_luminosity_comparison.png", dpi=240)
    plt.close()
    metrics = ["mean_tau_dimensionless", "max_tau_dimensionless", "mean_Psurv"]
    x = np.arange(len(metrics))
    width = 0.36
    plt.figure(figsize=(7.0, 4.2))
    plt.bar(x - width / 2, [uhe_results["fiducial_uhe_default"][m] for m in metrics], width, label="fiducial_uhe_default")
    plt.bar(x + width / 2, [uhe_results["collapsar_ndaf_like"][m] for m in metrics], width, label="collapsar_ndaf_like")
    plt.yscale("log")
    plt.xticks(x, ["mean tau\ndimensionless", "max tau\ndimensionless", "mean Psurv"])
    plt.title("Representative UHE opacity comparison, E=1e11 GeV\nfiducial uses historical UHE density rho0=1e-2")
    plt.legend(frameon=False)
    plt.tight_layout()
    plt.savefig(PLOT_DIR / "uhe_opacity_comparison.png", dpi=240)
    plt.close()
    plt.figure(figsize=(7.0, 4.4))
    for name, p in PRESETS.items():
        theta, r_mev = run_mev_surface_case(p)
        theta_deg = theta * 180 / PI
        plt.plot(theta_deg, r_mev, label=f"{name} MeV tau=2/3")
    for name in ["fiducial_uhe_default", "collapsar_ndaf_like"]:
        res = uhe_results[name]
        theta_deg = res["theta"] * 180 / PI
        plt.plot(theta_deg, res["r_uhe"], "--", label=f"{name} UHE tau=1")
    plt.xlabel(r"$\theta$ [deg]")
    plt.ylabel(r"$r_\tau/r_g$")
    plt.title("MeV neutrinosphere and UHE opacity-surface comparison")
    plt.legend(frameon=False, fontsize=7)
    plt.tight_layout()
    plt.savefig(PLOT_DIR / "neutrinosphere_comparison.png", dpi=240)
    plt.close()
    print(f"Saved: {OUT_DIR / 'comparison_summary.md'}")
    print(f"Saved: {PLOT_DIR / 'mev_luminosity_comparison.png'}")


if __name__ == "__main__":
    main()
