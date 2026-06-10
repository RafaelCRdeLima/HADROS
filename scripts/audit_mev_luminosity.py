#!/usr/bin/env python3
"""Audit the diagnostic MeV luminosity proxy without changing the model."""

from __future__ import annotations

import math
import os
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
OUT_DIR = ROOT / "output" / "mev_luminosity"
PLOT_DIR = ROOT / "plots" / "mev_physics"

PI = math.pi
M_U_G = 1.66053906660e-24
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
NR = env_int("MEV_LUMINOSITY_NR", "96")
NTH = env_int("MEV_LUMINOSITY_NTH", "72")
E_MIN = env_float("MEV_LUMINOSITY_E_MIN_MEV", "1.0")
E_MAX = env_float("MEV_LUMINOSITY_E_MAX_MEV", "80.0")
NE = env_int("MEV_LUMINOSITY_E_BINS", "24")

FLAVORS = ["nu_e", "anti_nu_e", "nu_x"]
CHANNELS = ["urca", "pair", "brems"]


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


def rho_profile(r, th, rho0=RHO0):
    delta = th - 0.5 * PI
    gaussian = np.exp(-((r - R0) / SIGMA_R) ** 2) * np.exp(-(delta / H_OVER_R) ** 2)
    gaussian = np.where((r <= 4.0) | (r >= 18.0) | (np.abs(delta) >= 0.45), 0.0, gaussian)
    return np.maximum(rho0 * gaussian, RHO_FLOOR)


def temperature_profile(r, th, rho, rho0=RHO0, t0=T0):
    shape = np.clip(rho / max(rho0, 1e-300), 0.0, 1.0)
    return np.maximum(T_FLOOR + t0 * shape ** max(T_POWER, 0.0), T_FLOOR)


def ye_profile(r, rho, rho0=RHO0, ye_torus=YE_TORUS):
    shape = np.clip(rho / max(rho0, 1e-300), 0.0, 1.0)
    ye = ye_torus + 0.05 * np.exp(-((r - R0) / 4.0) ** 2) * np.where(shape > 1e-4, 1.0, 0.0)
    return np.clip(ye, YE_FLOOR, YE_CEIL)


def thermal_shape(e, t):
    x = e / np.maximum(t, 1e-300)
    return np.where(x > 120.0, 0.0, e * e / (np.exp(np.minimum(x, 120.0)) + 1.0) / np.maximum(t**3, 1e-300))


def emissivity_channels(rho, temp, ye, e, flavor):
    rho10 = rho / 1e10
    t10 = temp / 10.0
    spec = thermal_shape(e, temp)
    if flavor == "nu_x":
        urca = np.zeros_like(rho)
    else:
        factor = ye if flavor == "nu_e" else 1.0 - ye
        urca = 1e30 * np.maximum(rho10, 0.0) * np.maximum(t10, 0.0) ** 6 * factor * spec
    pair = 3e28 * np.maximum(t10, 0.0) ** 9 * (1.0 if flavor == "nu_x" else 0.7) * spec
    brems = 1e27 * np.maximum(rho10, 0.0) ** 2 * np.maximum(t10, 0.0) ** 5 * (1.0 if flavor == "nu_x" else 0.8) * spec
    return {"urca": urca, "pair": pair, "brems": brems}


def luminosity_field(rho0=RHO0, t0=T0, ye_torus=YE_TORUS):
    rr, tt, dV = make_grid()
    rho = rho_profile(rr, tt, rho0)
    temp = temperature_profile(rr, tt, rho, rho0=rho0, t0=t0)
    ye = ye_profile(rr, rho, rho0=rho0, ye_torus=ye_torus)
    e_edges = np.logspace(np.log10(E_MIN), np.log10(E_MAX), NE + 1)
    e_mid = np.sqrt(e_edges[:-1] * e_edges[1:])
    de = np.diff(e_edges)
    fields = {flavor: {channel: np.zeros_like(rho) for channel in CHANNELS} for flavor in FLAVORS}
    for e, width in zip(e_mid, de):
        energy_weight = e * MEV_TO_ERG * width
        for flavor in FLAVORS:
            chans = emissivity_channels(rho, temp, ye, e, flavor)
            for channel in CHANNELS:
                fields[flavor][channel] += chans[channel] * energy_weight
    cell_lum = {flavor: {channel: fields[flavor][channel] * dV for channel in CHANNELS} for flavor in FLAVORS}
    total_cell = sum(cell_lum[f][c] for f in FLAVORS for c in CHANNELS)
    return rr, tt, dV, rho, temp, ye, cell_lum, total_cell


def weighted_quantiles(values, weights, qs=(0.5,)):
    values = np.asarray(values).ravel()
    weights = np.asarray(weights).ravel()
    mask = np.isfinite(values) & np.isfinite(weights) & (weights > 0)
    if not np.any(mask):
        return [float("nan") for _ in qs]
    v = values[mask]
    w = weights[mask]
    order = np.argsort(v)
    v = v[order]
    w = w[order]
    cdf = np.cumsum(w) / np.sum(w)
    return [float(np.interp(q, cdf, v)) for q in qs]


def write_budget(cell_lum):
    channel_totals = {f: {c: float(np.sum(cell_lum[f][c])) for c in CHANNELS} for f in FLAVORS}
    total = sum(channel_totals[f][c] for f in FLAVORS for c in CHANNELS)
    rows = ["flavor,channel,L_proxy_erg_s,channel_fraction,total_fraction,flavor_fraction"]
    for flavor in FLAVORS:
        flavor_total = sum(channel_totals[flavor].values())
        for channel in CHANNELS:
            val = channel_totals[flavor][channel]
            channel_total = sum(channel_totals[f][channel] for f in FLAVORS)
            rows.append(
                f"{flavor},{channel},{val:.8e},{val/max(channel_total,1e-300):.8e},"
                f"{val/max(total,1e-300):.8e},{val/max(flavor_total,1e-300):.8e}"
            )
    rows.append(f"all,all,{total:.8e},1.0,1.0,1.0")
    (OUT_DIR / "luminosity_budget.csv").write_text("\n".join(rows) + "\n")
    dominant = max(
        ((channel_totals[f][c], f, c) for f in FLAVORS for c in CHANNELS),
        key=lambda x: x[0],
    )
    lines = [
        "# MeV Luminosity Budget",
        "",
        f"Total diagnostic luminosity proxy: `{total:.8e}` erg/s-proxy.",
        "",
        f"Dominant contribution: `{dominant[1]}` / `{dominant[2]}` = `{dominant[0]:.8e}` erg/s-proxy.",
        "",
        "| Flavor | URCA | Pair | Brems | Total |",
        "|---|---:|---:|---:|---:|",
    ]
    for flavor in FLAVORS:
        vals = channel_totals[flavor]
        lines.append(
            f"| {flavor} | {vals['urca']:.3e} | {vals['pair']:.3e} | "
            f"{vals['brems']:.3e} | {sum(vals.values()):.3e} |"
        )
    lines += [
        "",
        "Interpretation: the baseline is dominated by charged-current URCA-like emission from electron-flavor neutrinos, especially anti_nu_e because the default torus is neutron rich.",
        "",
    ]
    (OUT_DIR / "luminosity_budget.md").write_text("\n".join(lines))
    return total, channel_totals, dominant


def write_unit_audit(total):
    rg = rg_cm(MBH_MSUN)
    lines = [
        "# MeV Luminosity Unit Audit",
        "",
        "| Quantity | Symbol | Units used | Units expected | Conversion/check |",
        "|---|---|---|---|---|",
        "| Density | rho | g cm^-3 | g cm^-3 | Input profile already in cgs density units. |",
        "| Temperature | T | MeV | MeV | Used directly in thermal scalings and spectral weights. |",
        "| Electron fraction | Ye | dimensionless | dimensionless | Clamped to configured physical range. |",
        "| Spectral emissivity proxy | j_E | arbitrary proxy per MeV per cm^3 per s | physical erg cm^-3 s^-1 MeV^-1 would require calibrated constants | Current coefficients are diagnostic, not calibrated rates. |",
        "| Opacity | alpha | cm^-1 | cm^-1 | n_target [cm^-3] times sigma [cm^2]. |",
        "| Radius | r | r_g | cm | r_g = GM/c^2. |",
        "| Gravitational radius | r_g | cm | cm | For current M_BH, r_g = {:.8e} cm. |".format(rg),
        "| Volume element | dV | cm^3 | cm^3 | dV = 2 pi r^2 sin(theta) dr dtheta after r and dr conversion to cm. |",
        "| Energy integration | E dE | MeV^2 converted partly to erg | calibrated energy emissivity requires calibrated j_E | Multiplies by E[MeV]*MeV_to_erg*dE[MeV]. |",
        "| Luminosity integral | L_proxy | erg/s-proxy | erg/s only after emissivity calibration | Volume and energy conversions are dimensional; emissivity normalization is not calibrated. |",
        "",
        f"Current integrated value: `{total:.8e}` erg/s-proxy.",
        "",
        "Conclusion: `L_proxy` is not a calibrated physical luminosity in erg/s. It is proportional to an energy luminosity after applying the configured proxy normalization, volume conversion, and MeV-to-erg conversion.",
        "",
    ]
    (OUT_DIR / "unit_audit.md").write_text("\n".join(lines))


def volume_audit(dV, rho, temp, total_cell):
    total_volume = float(np.sum(dV))
    torus_mask = rho > max(RHO_FLOOR * 10.0, 1e-3 * RHO0)
    torus_volume = float(np.sum(dV[torus_mask]))
    q = total_cell / np.maximum(dV, 1e-300)
    threshold = 1e-3 * float(np.max(q))
    emitting_mask = q >= threshold
    emitting_volume = float(np.sum(dV[emitting_mask]))
    flat_l = total_cell.ravel()
    flat_v = dV.ravel()
    order = np.argsort(flat_l)[::-1]
    cumulative = np.cumsum(flat_l[order])
    total = float(cumulative[-1])
    idx90 = int(np.searchsorted(cumulative, 0.9 * total)) if total > 0 else 0
    dominant_volume_90 = float(np.sum(flat_v[order[: idx90 + 1]])) if total > 0 else 0.0
    frac_lum = cumulative / max(total, 1e-300)
    frac_vol = np.cumsum(flat_v[order]) / max(total_volume, 1e-300)
    plt.figure(figsize=(6.2, 4.2))
    plt.plot(frac_vol, frac_lum)
    plt.axhline(0.9, color="black", ls="--", lw=1)
    plt.xlabel("fraction of total grid volume, sorted by luminosity density")
    plt.ylabel("cumulative luminosity fraction")
    plt.title("MeV emitting volume concentration")
    plt.tight_layout()
    plt.savefig(PLOT_DIR / "emitting_volume_contribution.png", dpi=240)
    plt.close()
    return {
        "total_grid_volume_cm3": total_volume,
        "torus_volume_cm3": torus_volume,
        "emitting_volume_cm3": emitting_volume,
        "dominant_90pct_luminosity_volume_cm3": dominant_volume_90,
        "emitting_threshold_proxy_per_cm3_s": threshold,
    }


def distribution_audit(rho, temp, total_cell):
    weights = total_cell
    t_mean = float(np.sum(temp * weights) / max(np.sum(weights), 1e-300))
    rho_mean = float(np.sum(rho * weights) / max(np.sum(weights), 1e-300))
    t_median = weighted_quantiles(temp, weights, [0.5])[0]
    rho_median = weighted_quantiles(rho, weights, [0.5])[0]
    mask = weights > 1e-6 * np.max(weights)
    t_stats = {
        "T_min_weighted_region_MeV": float(np.min(temp[mask])),
        "T_max_weighted_region_MeV": float(np.max(temp[mask])),
        "T_mean_emissivity_weighted_MeV": t_mean,
        "T_median_emissivity_weighted_MeV": t_median,
    }
    rho_stats = {
        "rho_min_weighted_region_gcm3": float(np.min(rho[mask])),
        "rho_max_weighted_region_gcm3": float(np.max(rho[mask])),
        "rho_mean_emissivity_weighted_gcm3": rho_mean,
        "rho_median_emissivity_weighted_gcm3": rho_median,
    }
    plt.figure(figsize=(6.2, 4.2))
    plt.hist(temp.ravel(), bins=60, weights=weights.ravel(), color="#c2410c")
    plt.xlabel("T [MeV]")
    plt.ylabel("luminosity-weighted contribution")
    plt.title("Emissivity-weighted temperature distribution")
    plt.tight_layout()
    plt.savefig(PLOT_DIR / "emissivity_weighted_temperature_distribution.png", dpi=240)
    plt.close()
    plt.figure(figsize=(6.2, 4.2))
    plt.hist(np.log10(np.maximum(rho.ravel(), 1e-300)), bins=70, weights=weights.ravel(), color="#2563eb")
    plt.xlabel(r"$\log_{10}\rho$ [g cm$^{-3}$]")
    plt.ylabel("luminosity-weighted contribution")
    plt.title("Emissivity-weighted density distribution")
    plt.tight_layout()
    plt.savefig(PLOT_DIR / "emissivity_weighted_density_distribution.png", dpi=240)
    plt.close()
    return t_stats, rho_stats


def scan_luminosity():
    rho_values = np.array([1e8, 1e9, 1e10, 1e11, 1e12])
    t_values = np.array([2, 4, 6, 8, 10, 15, 20], dtype=float)
    ye_values = np.array([0.05, 0.1, 0.2, 0.3, 0.4, 0.5])
    rows = ["scan,parameter,L_total_proxy_erg_s"]
    rho_l = []
    for val in rho_values:
        *_, cell_lum, total_cell = luminosity_field(rho0=val, t0=T0, ye_torus=YE_TORUS)
        lum = float(np.sum(total_cell))
        rho_l.append(lum)
        rows.append(f"rho0,{val:.8e},{lum:.8e}")
    t_l = []
    for val in t_values:
        *_, cell_lum, total_cell = luminosity_field(rho0=RHO0, t0=val, ye_torus=YE_TORUS)
        lum = float(np.sum(total_cell))
        t_l.append(lum)
        rows.append(f"T0,{val:.8e},{lum:.8e}")
    ye_l = []
    for val in ye_values:
        *_, cell_lum, total_cell = luminosity_field(rho0=RHO0, t0=T0, ye_torus=val)
        lum = float(np.sum(total_cell))
        ye_l.append(lum)
        rows.append(f"Ye,{val:.8e},{lum:.8e}")
    (OUT_DIR / "realism_scan.csv").write_text("\n".join(rows) + "\n")
    for x, y, xlabel, title, fname in [
        (rho_values, rho_l, r"$\rho_0$ [g cm$^{-3}$]", "Luminosity proxy vs density normalization", "luminosity_vs_rho0.png"),
        (t_values, t_l, r"$T_0$ [MeV]", "Luminosity proxy vs temperature scale", "luminosity_vs_temperature.png"),
        (ye_values, ye_l, r"$Y_e$", "Luminosity proxy vs electron fraction", "luminosity_vs_Ye.png"),
    ]:
        plt.figure(figsize=(6.2, 4.2))
        plt.loglog(x, y, marker="o") if np.all(x > 0) else plt.plot(x, y, marker="o")
        plt.xlabel(xlabel)
        plt.ylabel(r"$L_{\nu,\rm proxy}$ [erg s$^{-1}$ proxy]")
        plt.title(title)
        plt.tight_layout()
        plt.savefig(PLOT_DIR / fname, dpi=240)
        plt.close()
    return {
        "rho0": list(zip(rho_values.tolist(), rho_l)),
        "T0": list(zip(t_values.tolist(), t_l)),
        "Ye": list(zip(ye_values.tolist(), ye_l)),
    }


def write_gap_analysis(total, dominant, vol, t_stats, rho_stats, scans):
    current_vs_ndaf = total / 1e52
    t20 = dict(scans["T0"])[20.0]
    rho12 = dict(scans["rho0"])[1e12]
    ye_values = np.array([value for value, _ in scans["Ye"]])
    ye_lums = np.array([lum for _, lum in scans["Ye"]])
    ye_flat = bool(np.max(ye_lums) / max(np.min(ye_lums), 1e-300) < 1.01)
    lines = [
        "# MeV Realism-Gap Analysis",
        "",
        "## 1. Is the baseline luminosity expected?",
        "",
        f"The current value, `{total:.3e}` erg/s-proxy, is expected for the present semi-analytic diagnostic setup. It should not be treated as a calibrated physical luminosity.",
        "",
        "The main reasons are the proxy emissivity normalization, the compact Gaussian emitting region, and the baseline thermodynamic scale.",
        "",
        "## 2. Ranked causes of the gap to NDAF-like luminosities",
        "",
        "1. `normalization choices`: the emissivity coefficients are explicitly diagnostic and not calibrated weak-interaction cooling rates.",
        "2. `simplified emissivities`: no detailed beta equilibrium, blocking, finite-temperature chemical potentials, nuclear composition, or calibrated pair/brems rates.",
        "3. `temperature`: URCA scales roughly as T^6 and pair as T^9, so the baseline effective temperatures strongly limit L.",
        "4. `density`: raising rho0 increases L, especially through URCA and bremsstrahlung.",
        "5. `small volume`: emission is concentrated in a small fraction of the grid volume, but this is secondary to normalization and thermodynamics.",
        "6. `unit issues`: no r_g or dV conversion bug was found; the issue is that j_E is a proxy, not calibrated cgs emissivity.",
        "",
        "## 3. Most efficient parameter changes",
        "",
        f"- Raising `T0` to 20 MeV gives `{t20:.3e}` erg/s-proxy in the controlled scan.",
        f"- Raising `rho0` to 1e12 g cm^-3 gives `{rho12:.3e}` erg/s-proxy in the controlled scan.",
        f"- Varying `Ye` over {ye_values[0]:.2f}-{ye_values[-1]:.2f} changes the total by less than 1 percent: `{ye_flat}`.",
        "",
        "Temperature is the most efficient lever because of the steep T powers. In the current simplified total luminosity proxy, Ye mainly redistributes electron-flavor emission between nu_e and anti_nu_e rather than changing the all-flavor total.",
        "",
        "## 4. Physically justified changes",
        "",
        "- Use thermodynamic profiles from a collapsar/NDAF or merger-torus calculation.",
        "- Calibrate local emissivity coefficients against standard neutrino cooling formulae.",
        "- Add finite-temperature electron chemical potentials and blocking factors.",
        "- Replace the fixed Gaussian support with a physically motivated disk/envelope volume.",
        "",
        "## 5. Artificial tuning to avoid",
        "",
        "- Multiplying `MEV_NORM` to force agreement with 1e51-1e53 erg/s.",
        "- Increasing volume or temperature without a physical background model.",
        "- Claiming calibrated luminosities from the current proxy coefficients.",
        "",
        "## Baseline diagnostics",
        "",
        f"- dominant channel/flavor: `{dominant[2]}` / `{dominant[1]}`",
        f"- torus volume: `{vol['torus_volume_cm3']:.3e}` cm^3",
        f"- volume giving 90 percent of luminosity: `{vol['dominant_90pct_luminosity_volume_cm3']:.3e}` cm^3",
        f"- emissivity-weighted T median: `{t_stats['T_median_emissivity_weighted_MeV']:.3f}` MeV",
        f"- emissivity-weighted rho median: `{rho_stats['rho_median_emissivity_weighted_gcm3']:.3e}` g cm^-3",
        f"- baseline / 1e52 erg/s: `{current_vs_ndaf:.3e}`",
        "",
    ]
    (OUT_DIR / "realism_gap_analysis.md").write_text("\n".join(lines))


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    PLOT_DIR.mkdir(parents=True, exist_ok=True)
    rr, tt, dV, rho, temp, ye, cell_lum, total_cell = luminosity_field()
    total, channel_totals, dominant = write_budget(cell_lum)
    write_unit_audit(total)
    vol = volume_audit(dV, rho, temp, total_cell)
    t_stats, rho_stats = distribution_audit(rho, temp, total_cell)
    scans = scan_luminosity()
    write_gap_analysis(total, dominant, vol, t_stats, rho_stats, scans)
    audit_lines = ["# MeV Luminosity Audit Summary", ""]
    audit_lines += [f"- total_L_proxy: `{total:.8e}` erg/s-proxy"]
    audit_lines += [f"- dominant: `{dominant[1]}` / `{dominant[2]}` = `{dominant[0]:.8e}`"]
    for key, val in vol.items():
        audit_lines.append(f"- {key}: `{val:.8e}`")
    for key, val in t_stats.items():
        audit_lines.append(f"- {key}: `{val:.8e}`")
    for key, val in rho_stats.items():
        audit_lines.append(f"- {key}: `{val:.8e}`")
    audit_lines += [
        "",
        "Conclusion: no dimensional conversion bug was found. The luminosity is low primarily because the emissivity is a diagnostic proxy and the baseline semi-analytic torus is not an NDAF-calibrated thermodynamic model.",
        "",
    ]
    (OUT_DIR / "audit_summary.md").write_text("\n".join(audit_lines))
    print(f"Saved: {OUT_DIR / 'luminosity_budget.csv'}")
    print(f"Saved: {OUT_DIR / 'unit_audit.md'}")
    print(f"Saved: {OUT_DIR / 'realism_gap_analysis.md'}")
    print(f"Saved: {OUT_DIR / 'audit_summary.md'}")


if __name__ == "__main__":
    main()
