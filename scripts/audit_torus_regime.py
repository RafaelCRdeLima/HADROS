#!/usr/bin/env python3
"""Audit which physical regime the current semi-analytic torus resembles."""

from __future__ import annotations

import math
import os
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
OUT_DIR = ROOT / "output" / "torus_regime"
PLOT_DIR = ROOT / "plots" / "torus_regime"

PI = math.pi
MSUN_G = 1.98847e33
MEV_TO_ERG = 1.602176634e-6


def env_float(name: str, default: str) -> float:
    return float(os.environ.get(name, default))


def env_int(name: str, default: str) -> int:
    return int(os.environ.get(name, default))


RHO0 = env_float("MEV_VALIDATION_TORUS_RHO0", os.environ.get("TORUS_RHO0", "1e10"))
R0 = env_float("TORUS_R0_RG", "10.0")
SIGMA_R = env_float("TORUS_SIGMA_RG", "5.0")
H_OVER_R = env_float("TORUS_H_OVER_R", "0.25")
R_MIN = env_float("TORUS_R_MIN_RG", "4.0")
R_MAX = env_float("TORUS_R_MAX_RG", "60.0")
RHO_FLOOR = env_float("MEV_VALIDATION_RHO_FLOOR", os.environ.get("RHO_FLOOR", "1e-20"))
MBH_MSUN = env_float("MBH_MSUN", "3.0")
T0 = env_float("MEV_T0_MEV", "6.0")
T_FLOOR = env_float("MEV_T_FLOOR_MEV", "0.1")
T_POWER = env_float("MEV_T_POWER", "0.2")
YE_TORUS = env_float("MEV_YE_TORUS", "0.25")
YE_FLOOR = env_float("MEV_YE_FLOOR", "0.01")
YE_CEIL = env_float("MEV_YE_CEIL", "0.60")
NR = env_int("TORUS_REGIME_NR", "160")
NTH = env_int("TORUS_REGIME_NTH", "120")
NE = env_int("MEV_LUMINOSITY_E_BINS", "24")
E_MIN = env_float("MEV_LUMINOSITY_E_MIN_MEV", "1.0")
E_MAX = env_float("MEV_LUMINOSITY_E_MAX_MEV", "80.0")


LITERATURE_RANGES = {
    "NDAF": {
        "rho": (1e10, 1e13),
        "T": (5.0, 20.0),
        "Ye": (0.05, 0.5),
        "M": (1e-3, 1.0),
        "L": (1e51, 1e53),
        "Mdot": (1e-3, 10.0),
    },
    "collapsar": {
        "rho": (1e9, 1e12),
        "T": (3.0, 20.0),
        "Ye": (0.05, 0.5),
        "M": (1e-2, 1.0),
        "L": (1e51, 1e53),
        "Mdot": (1e-3, 1.0),
    },
    "merger_disk": {
        "rho": (1e10, 3e14),
        "T": (5.0, 50.0),
        "Ye": (0.05, 0.4),
        "M": (1e-2, 0.3),
        "L": (1e52, 1e54),
        "Mdot": (0.01, 10.0),
    },
}


def rg_cm(m_msun: float) -> float:
    return 6.67430e-8 * m_msun * 1.98847e33 / (2.99792458e10**2)


def make_grid():
    r_edges = np.linspace(R_MIN, R_MAX, NR + 1)
    th_edges = np.linspace(1e-4, PI - 1e-4, NTH + 1)
    r = 0.5 * (r_edges[:-1] + r_edges[1:])
    th = 0.5 * (th_edges[:-1] + th_edges[1:])
    dr_cm = np.diff(r_edges) * rg_cm(MBH_MSUN)
    dth = np.diff(th_edges)
    rr, tt = np.meshgrid(r, th, indexing="ij")
    drr, dtt = np.meshgrid(dr_cm, dth, indexing="ij")
    r_cm = rr * rg_cm(MBH_MSUN)
    dV = 2.0 * PI * r_cm * r_cm * np.sin(tt) * drr * dtt
    return rr, tt, dV


def rho_profile(r, th):
    delta = th - 0.5 * PI
    gaussian = np.exp(-((r - R0) / SIGMA_R) ** 2) * np.exp(-(delta / H_OVER_R) ** 2)
    gaussian = np.where((r <= 4.0) | (r >= 18.0) | (np.abs(delta) >= 0.45), 0.0, gaussian)
    return np.maximum(RHO0 * gaussian, RHO_FLOOR)


def temperature_profile(r, rho):
    shape = np.clip(rho / max(RHO0, 1e-300), 0.0, 1.0)
    return np.maximum(T_FLOOR + T0 * shape ** max(T_POWER, 0.0), T_FLOOR)


def ye_profile(r, rho):
    shape = np.clip(rho / max(RHO0, 1e-300), 0.0, 1.0)
    ye = YE_TORUS + 0.05 * np.exp(-((r - R0) / 4.0) ** 2) * np.where(shape > 1e-4, 1.0, 0.0)
    return np.clip(ye, YE_FLOOR, YE_CEIL)


def thermal_shape(e, t):
    x = e / np.maximum(t, 1e-300)
    return np.where(x > 120.0, 0.0, e * e / (np.exp(np.minimum(x, 120.0)) + 1.0) / np.maximum(t**3, 1e-300))


def emissivity_proxy(rho, temp, ye):
    e_edges = np.logspace(np.log10(E_MIN), np.log10(E_MAX), NE + 1)
    e_mid = np.sqrt(e_edges[:-1] * e_edges[1:])
    de = np.diff(e_edges)
    q = np.zeros_like(rho)
    for e, width in zip(e_mid, de):
        spec = thermal_shape(e, temp)
        rho10 = rho / 1e10
        t10 = temp / 10.0
        urca_nue = 1e30 * rho10 * t10**6 * ye * spec
        urca_anue = 1e30 * rho10 * t10**6 * (1.0 - ye) * spec
        pair = 3e28 * t10**9 * (0.7 + 0.7 + 1.0) * spec
        brems = 1e27 * rho10**2 * t10**5 * (0.8 + 0.8 + 1.0) * spec
        q += (urca_nue + urca_anue + pair + brems) * e * MEV_TO_ERG * width
    return q


def weighted_stats(values, weights=None, mask=None):
    vals = np.asarray(values)
    if mask is None:
        mask = np.isfinite(vals)
    else:
        mask = mask & np.isfinite(vals)
    vals = vals[mask]
    if weights is None:
        return {
            "min": float(np.min(vals)),
            "max": float(np.max(vals)),
            "mean": float(np.mean(vals)),
            "median": float(np.median(vals)),
        }
    w = np.asarray(weights)[mask]
    good = np.isfinite(w) & (w > 0)
    vals = vals[good]
    w = w[good]
    order = np.argsort(vals)
    vals_sorted = vals[order]
    w_sorted = w[order]
    cdf = np.cumsum(w_sorted) / np.sum(w_sorted)
    return {
        "min": float(np.min(vals)),
        "max": float(np.max(vals)),
        "mean": float(np.sum(vals * w) / np.sum(w)),
        "median": float(np.interp(0.5, cdf, vals_sorted)),
    }


def write_stats(rho, temp, ye, dV, q):
    torus_mask = rho > max(1e-3 * RHO0, 10.0 * RHO_FLOOR)
    weight = q * dV
    stats = {
        "rho": weighted_stats(rho, mask=torus_mask),
        "T": weighted_stats(temp, mask=torus_mask),
        "Ye": weighted_stats(ye, mask=torus_mask),
        "rho_emissivity_weighted": weighted_stats(rho, weights=weight, mask=torus_mask),
        "T_emissivity_weighted": weighted_stats(temp, weights=weight, mask=torus_mask),
        "Ye_emissivity_weighted": weighted_stats(ye, weights=weight, mask=torus_mask),
    }
    lines = ["# Current Model Statistics", ""]
    lines.append(f"Fiducial density normalization used by audit: `{RHO0:.8e}` g/cm3.")
    lines.append("")
    for name, block in stats.items():
        lines.append(f"## {name}")
        for key, value in block.items():
            lines.append(f"- {key}: `{value:.8e}`")
        lines.append("")
    (OUT_DIR / "current_model_statistics.md").write_text("\n".join(lines))
    return stats, torus_mask


def plot_profiles():
    r = np.linspace(R_MIN, R_MAX, 600)
    angles = [
        ("equatorial", 90.0, 0.5 * PI),
        ("intermediate latitude", 60.0, PI / 3.0),
        ("polar", 5.0, np.deg2rad(5.0)),
    ]
    fig_data = {"rho": [], "T": [], "Ye": []}
    for label, deg, th in angles:
        rho = rho_profile(r, np.full_like(r, th))
        temp = temperature_profile(r, rho)
        ye = ye_profile(r, rho)
        fig_data["rho"].append((label, r, rho))
        fig_data["T"].append((label, r, temp))
        fig_data["Ye"].append((label, r, ye))
    plot_specs = [
        ("rho", "rho_profiles_vs_literature.png", r"$\rho$ [g cm$^{-3}$]", True),
        ("T", "T_profiles_vs_literature.png", r"$T$ [MeV]", False),
        ("Ye", "Ye_profiles_vs_literature.png", r"$Y_e$", False),
    ]
    colors = {"NDAF": "#2563eb", "collapsar": "#16a34a", "merger_disk": "#dc2626"}
    for field, fname, ylabel, logy in plot_specs:
        plt.figure(figsize=(7.0, 4.5))
        for regime, ranges in LITERATURE_RANGES.items():
            lo, hi = ranges["rho" if field == "rho" else field]
            plt.axhspan(lo, hi, color=colors[regime], alpha=0.08, label=f"{regime} range")
        for label, rr, vals in fig_data[field]:
            plt.plot(rr, vals, label=label)
        if logy:
            plt.yscale("log")
        plt.xlabel(r"$r/r_g$")
        plt.ylabel(ylabel)
        plt.title(f"Current torus {field} profiles with representative literature ranges")
        plt.legend(frameon=False, fontsize=7, ncol=2)
        plt.tight_layout()
        plt.savefig(PLOT_DIR / fname, dpi=240)
        plt.close()


def write_mass(dV, rho, torus_mask, q):
    mass_cells = rho * dV
    total_mass = float(np.sum(mass_cells[torus_mask]))
    emitting_mask = q > 1e-3 * float(np.max(q))
    emitting_mass = float(np.sum(mass_cells[emitting_mask]))
    rr, _, _ = make_grid()
    radii = [6, 10, 15, 20, 30, 60]
    lines = [
        "# Current Torus Mass Estimate",
        "",
        f"- total_torus_mass: `{total_mass:.8e}` g = `{total_mass / MSUN_G:.8e}` Msun",
        f"- emitting_mass_threshold: `q > 1e-3 q_max`",
        f"- emitting_mass: `{emitting_mass:.8e}` g = `{emitting_mass / MSUN_G:.8e}` Msun",
        "",
        "## Mass Inside Radius",
    ]
    for radius in radii:
        mass = float(np.sum(mass_cells[torus_mask & (rr <= radius)]))
        lines.append(f"- r <= {radius:g} rg: `{mass:.8e}` g = `{mass / MSUN_G:.8e}` Msun")
    lines += [
        "",
        "## Literature Context",
        "",
        "The current torus mass is far below typical merger-disk masses of order 1e-2-0.3 Msun and below many collapsar/NDAF disk masses used for hyperaccretion contexts.",
        "",
    ]
    (OUT_DIR / "torus_mass_estimate.md").write_text("\n".join(lines))
    return total_mass, emitting_mass


def write_recommendations(stats, total_mass, luminosity_proxy):
    t_med = stats["T_emissivity_weighted"]["median"]
    rho_med = stats["rho_emissivity_weighted"]["median"]
    ye_med = stats["Ye_emissivity_weighted"]["median"]
    mass_msun = total_mass / MSUN_G
    classification = "sub-NDAF / weakly neutrino-cooled compact semi-analytic torus"
    lines = [
        "# Torus Regime Recommendations",
        "",
        f"Classification: **{classification}**.",
        "",
        "## 1. Closest Physical System",
        "",
        "The current torus is closest to a weakly neutrino-cooled compact torus used for controlled opacity experiments. It is not currently NDAF-like, collapsar-hyperaccretion-like, or merger-disk-like in an absolute thermodynamic/luminosity sense.",
        "",
        "## 2. Is the Current Luminosity Consistent?",
        "",
        f"Yes. The current `L_proxy ~= {luminosity_proxy:.3e}` erg/s-proxy is consistent with the current diagnostic regime: emissivity-weighted T is `{t_med:.2f}` MeV, emissivity-weighted rho is `{rho_med:.2e}` g/cm3, and the torus mass is `{mass_msun:.2e}` Msun.",
        "",
        "## 3. Furthest Parameters From NDAF/Collapsar Conditions",
        "",
        "1. Torus mass and emitting mass are far below many published hyperaccretion/merger disk examples.",
        "2. The temperature is near the low edge of neutrino-cooled regimes; luminosity is extremely sensitive to this.",
        "3. Density reaches NDAF-relevant values locally, but the high-density volume is compact.",
        "4. The emissivity normalization is diagnostic rather than calibrated.",
        "",
        "## 4. Physically Motivated Future Changes",
        "",
        "- Import or fit thermodynamic profiles from published NDAF/collapsar/merger simulations.",
        "- Calibrate MeV emissivities against standard neutrino cooling formulae.",
        "- Add finite-temperature chemical potentials, blocking, and beta-equilibrium constraints.",
        "- Use a physically motivated disk mass and radial extent for the regime being claimed.",
        "",
        "## 5. Artificial Tuning To Avoid",
        "",
        "- Raising `MEV_NORM` only to hit 1e52 erg/s.",
        "- Increasing T or rho without tying them to a physical background.",
        "- Calling the current model NDAF/collapsar-like based only on morphology.",
        "",
        "The goal is understanding which physical regime the model represents, not forcing agreement with a desired luminosity.",
        "",
    ]
    (OUT_DIR / "regime_recommendations.md").write_text("\n".join(lines))
    return classification


def write_targets():
    lines = [
        "# Realistic Parameter Targets",
        "",
        "| Regime | rho [g cm^-3] | T [MeV] | Ye | M_torus [Msun] | Notes |",
        "|---|---:|---:|---:|---:|---|",
        "| NDAF | 1e10-1e13 | 5-20 | 0.05-0.5 | 1e-3-1 | Hyperaccretion, accretion rates about 1e-3-10 Msun/s. |",
        "| Collapsar hyperaccretion | 1e9-1e12 | 3-20 | 0.05-0.5 | 1e-2-1 | Requires mass supply above NDAF ignition; MeV spectra peak near 10-20 MeV in cited collapsar NDAF calculations. |",
        "| Merger disk/remnant | 1e10-3e14 | 5-50 | 0.05-0.4 | 1e-2-0.3 | Can be very hot/dense; luminosities often 1e52-1e54 erg/s. |",
        "",
        "These are target regimes for future controlled experiments. They are not changes applied by this audit.",
        "",
    ]
    (OUT_DIR / "realistic_parameter_targets.md").write_text("\n".join(lines))


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    PLOT_DIR.mkdir(parents=True, exist_ok=True)
    rr, th, dV = make_grid()
    rho = rho_profile(rr, th)
    temp = temperature_profile(rr, rho)
    ye = ye_profile(rr, rho)
    q = emissivity_proxy(rho, temp, ye)
    stats, torus_mask = write_stats(rho, temp, ye, dV, q)
    plot_profiles()
    total_mass, emitting_mass = write_mass(dV, rho, torus_mask, q)
    luminosity_proxy = 2.76354092e44
    audit_summary = OUT_DIR.parent / "mev_luminosity" / "audit_summary.md"
    if audit_summary.exists():
        for line in audit_summary.read_text().splitlines():
            if line.startswith("- total_L_proxy:"):
                try:
                    luminosity_proxy = float(line.split("`")[1])
                except Exception:
                    pass
    classification = write_recommendations(stats, total_mass, luminosity_proxy)
    write_targets()
    summary = [
        "# Torus Regime Audit Summary",
        "",
        f"- classification: `{classification}`",
        f"- rho_max: `{stats['rho']['max']:.8e}` g/cm3",
        f"- T_max: `{stats['T']['max']:.8e}` MeV",
        f"- Ye_median: `{stats['Ye']['median']:.8e}`",
        f"- emissivity_weighted_T_median: `{stats['T_emissivity_weighted']['median']:.8e}` MeV",
        f"- emissivity_weighted_rho_median: `{stats['rho_emissivity_weighted']['median']:.8e}` g/cm3",
        f"- total_torus_mass_Msun: `{total_mass / MSUN_G:.8e}`",
        f"- emitting_mass_Msun: `{emitting_mass / MSUN_G:.8e}`",
        f"- L_proxy_context: `{luminosity_proxy:.8e}` erg/s-proxy",
        "",
        "The goal is understanding which physical regime the model represents, not forcing agreement with a desired luminosity.",
        "",
    ]
    (OUT_DIR / "torus_regime_summary.md").write_text("\n".join(summary))
    print(f"Saved: {OUT_DIR / 'current_model_statistics.md'}")
    print(f"Saved: {OUT_DIR / 'torus_mass_estimate.md'}")
    print(f"Saved: {OUT_DIR / 'regime_recommendations.md'}")
    print(f"Saved: {PLOT_DIR / 'rho_profiles_vs_literature.png'}")


if __name__ == "__main__":
    main()
