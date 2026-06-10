#!/usr/bin/env python3
"""Compute diagnostic volume-integrated MeV luminosity proxies."""

from __future__ import annotations

import math
import os
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
OUT_DIR = ROOT / "output" / "mev_luminosity"
PLOT_DIR = ROOT / "plots" / "mev_physics"
VAL_DIR = ROOT / "output" / "validation"

M_U_G = 1.66053906660e-24
MEV_TO_ERG = 1.602176634e-6
PI = math.pi


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
YE_FUNNEL = env_float("MEV_YE_FUNNEL", "0.55")
YE_ENVELOPE = env_float("MEV_YE_ENVELOPE", "0.45")
YE_FLOOR = env_float("MEV_YE_FLOOR", "0.01")
YE_CEIL = env_float("MEV_YE_CEIL", "0.60")
NR = env_int("MEV_LUMINOSITY_NR", "96")
NTH = env_int("MEV_LUMINOSITY_NTH", "72")
E_MIN = env_float("MEV_LUMINOSITY_E_MIN_MEV", "1.0")
E_MAX = env_float("MEV_LUMINOSITY_E_MAX_MEV", "80.0")
NE = env_int("MEV_LUMINOSITY_E_BINS", "24")


def rg_cm(m_msun: float) -> float:
    return 6.67430e-8 * m_msun * 1.98847e33 / (2.99792458e10**2)


def rho_profile(r, th, rho0=RHO0):
    delta = th - 0.5 * PI
    gaussian = np.exp(-((r - R0) / SIGMA_R) ** 2) * np.exp(-(delta / H_OVER_R) ** 2)
    gaussian = np.where((r <= 4.0) | (r >= 18.0) | (np.abs(delta) >= 0.45), 0.0, gaussian)
    return np.maximum(rho0 * gaussian, RHO_FLOOR)


def temperature_profile(r, th, rho, t_scale=1.0):
    shape = np.clip(rho / max(RHO0, 1e-300), 0.0, 1.0)
    return np.maximum(T_FLOOR + t_scale * T0 * shape ** max(T_POWER, 0.0), T_FLOOR)


def ye_profile(r, th, rho):
    shape = np.clip(rho / max(RHO0, 1e-300), 0.0, 1.0)
    ye = YE_TORUS + 0.05 * np.exp(-((r - R0) / 4.0) ** 2) * np.where(shape > 1e-4, 1.0, 0.0)
    return np.clip(ye, YE_FLOOR, YE_CEIL)


def thermal_shape(e, t):
    x = e / np.maximum(t, 1e-300)
    return np.where(x > 120.0, 0.0, e * e / (np.exp(np.minimum(x, 120.0)) + 1.0) / np.maximum(t**3, 1e-300))


def channels(rho, t, ye, e, flavor):
    rho10 = rho / 1e10
    t10 = t / 10.0
    spec = thermal_shape(e, t)
    if flavor == "nu_x":
        urca = np.zeros_like(rho)
    else:
        ff = ye if flavor == "nu_e" else 1.0 - ye
        urca = 1e30 * np.maximum(rho10, 0.0) * np.maximum(t10, 0.0) ** 6 * ff * spec
    pair = 3e28 * np.maximum(t10, 0.0) ** 9 * (1.0 if flavor == "nu_x" else 0.7) * spec
    brems = 1e27 * np.maximum(rho10, 0.0) ** 2 * np.maximum(t10, 0.0) ** 5 * (1.0 if flavor == "nu_x" else 0.8) * spec
    return urca, pair, brems


def integrate_luminosity(rho_scale=1.0, t_scale=1.0):
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
    rho = rho_profile(rr, tt, RHO0 * rho_scale)
    temp = temperature_profile(rr, tt, rho, t_scale)
    ye = ye_profile(rr, tt, rho)
    e_edges = np.logspace(np.log10(E_MIN), np.log10(E_MAX), NE + 1)
    e_mid = np.sqrt(e_edges[:-1] * e_edges[1:])
    de = np.diff(e_edges)
    out = {}
    for flavor in ["nu_e", "anti_nu_e", "nu_x"]:
        ch_sum = {"urca": 0.0, "pair": 0.0, "brems": 0.0}
        for e, width in zip(e_mid, de):
            urca, pair, brems = channels(rho, temp, ye, e, flavor)
            energy_weight = e * MEV_TO_ERG * width
            ch_sum["urca"] += float(np.sum(urca * energy_weight * dV))
            ch_sum["pair"] += float(np.sum(pair * energy_weight * dV))
            ch_sum["brems"] += float(np.sum(brems * energy_weight * dV))
        out[flavor] = ch_sum
    return out


def total_for(result, flavor):
    return sum(result[flavor].values())


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    PLOT_DIR.mkdir(parents=True, exist_ok=True)
    VAL_DIR.mkdir(parents=True, exist_ok=True)

    base = integrate_luminosity()
    rows = ["flavor,urca_proxy_erg_s,pair_proxy_erg_s,brems_proxy_erg_s,total_proxy_erg_s"]
    for flavor in ["nu_e", "anti_nu_e", "nu_x"]:
        vals = base[flavor]
        rows.append(f"{flavor},{vals['urca']:.8e},{vals['pair']:.8e},{vals['brems']:.8e},{total_for(base, flavor):.8e}")
    total = sum(total_for(base, f) for f in base)
    rows.append(f"all_flavors,,,,{total:.8e}")
    csv_path = OUT_DIR / "mev_luminosity_summary.csv"
    csv_path.write_text("\n".join(rows) + "\n")

    md_path = OUT_DIR / "mev_luminosity_summary.md"
    md_path.write_text(
        "\n".join(
            [
                "# Diagnostic MeV Luminosity Proxy",
                "",
                "This is a post-processing volume integral of the approximate local MeV energy emissivity proxy.",
                "It is not a calibrated absolute neutrino luminosity prediction.",
                "",
                f"- grid: NR={NR}, NTH={NTH}, NE={NE}",
                f"- energy range: {E_MIN:g}-{E_MAX:g} MeV",
                f"- total_L_proxy: {total:.8e} erg/s proxy",
                "",
                "NDAF/collapsar literature often discusses neutrino luminosities around 1e51-1e53 erg/s for high accretion-rate disks.",
                "Do not force this diagnostic proxy to match that range by arbitrary normalization.",
                "",
            ]
        )
    )

    t_scales = np.array([0.5, 0.75, 1.0, 1.5, 2.0])
    l_t = np.array([sum(total_for(integrate_luminosity(t_scale=s), f) for f in ["nu_e", "anti_nu_e", "nu_x"]) for s in t_scales])
    plt.figure(figsize=(6.2, 4.2))
    plt.loglog(T0 * t_scales, l_t, marker="o")
    plt.xlabel(r"$T_0$ scale [MeV]")
    plt.ylabel(r"$L_\nu$ proxy [erg s$^{-1}$]")
    plt.title("Diagnostic MeV luminosity proxy increases strongly with temperature")
    plt.tight_layout()
    plt.savefig(PLOT_DIR / "mev_luminosity_vs_temperature.png", dpi=240)
    plt.close()

    rho_scales = np.array([0.1, 0.3, 1.0, 3.0, 10.0])
    l_rho = np.array([sum(total_for(integrate_luminosity(rho_scale=s), f) for f in ["nu_e", "anti_nu_e", "nu_x"]) for s in rho_scales])
    plt.figure(figsize=(6.2, 4.2))
    plt.loglog(RHO0 * rho_scales, l_rho, marker="o")
    plt.xlabel(r"$\rho_0$ [g cm$^{-3}$]")
    plt.ylabel(r"$L_\nu$ proxy [erg s$^{-1}$]")
    plt.title("Diagnostic MeV luminosity proxy increases with density normalization")
    plt.tight_layout()
    plt.savefig(PLOT_DIR / "mev_luminosity_vs_density.png", dpi=240)
    plt.close()

    labels = []
    values = []
    for flavor in ["nu_e", "anti_nu_e", "nu_x"]:
        for channel in ["urca", "pair", "brems"]:
            labels.append(f"{flavor}\n{channel}")
            values.append(base[flavor][channel])
    plt.figure(figsize=(8.2, 4.4))
    plt.bar(np.arange(len(values)), values, color="#4776b4")
    plt.yscale("log")
    plt.xticks(np.arange(len(values)), labels, fontsize=7)
    plt.ylabel(r"channel contribution [erg s$^{-1}$ proxy]")
    plt.title("Diagnostic MeV luminosity channel breakdown")
    plt.tight_layout()
    plt.savefig(PLOT_DIR / "mev_luminosity_channel_breakdown.png", dpi=240)
    plt.close()

    finite_nonnegative = all(v >= 0 and np.isfinite(v) for f in base.values() for v in f.values())
    density_increases = bool(l_rho[-1] > l_rho[0])
    temp_increases = bool(l_t[-1] > l_t[0])
    fractions_ok = bool(np.isclose(sum(values), total, rtol=1e-12))
    val_path = VAL_DIR / "mev_luminosity_validation.txt"
    val_path.write_text(
        "\n".join(
            [
                "Diagnostic MeV luminosity proxy validation",
                "status: complete_for_diagnostic_proxy",
                f"finite_nonnegative_luminosities: {finite_nonnegative}",
                f"luminosity_increases_with_rho_normalization: {density_increases}",
                f"luminosity_increases_strongly_with_T: {temp_increases}",
                f"channel_fractions_sum_correctly: {fractions_ok}",
                f"total_L_proxy_erg_s: {total:.8e}",
                "calibration_status: diagnostic_proxy_not_calibrated_absolute_luminosity",
                "",
            ]
        )
    )
    print(f"Saved: {csv_path}")
    print(f"Saved: {md_path}")
    print(f"Saved: {val_path}")


if __name__ == "__main__":
    main()
