#!/usr/bin/env python3
"""Validate diagnostic electron degeneracy approximations for MeV physics."""

from __future__ import annotations

import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
PLOT_DIR = ROOT / "plots" / "mev_physics"
OUT_DIR = ROOT / "output" / "validation"

M_U_G = 1.66053906660e-24
HBAR_C_MEV_CM = 1.973269804e-11
M_E_MEV = 0.510998950


def electron_number_density_cm3(rho, ye):
    return np.maximum(rho, 0.0) * np.clip(ye, 0.0, 1.0) / M_U_G


def electron_fermi_momentum_mev(rho, ye):
    ne = electron_number_density_cm3(rho, ye)
    return HBAR_C_MEV_CM * np.cbrt(3.0 * math.pi * math.pi * ne)


def electron_chemical_potential_mev(rho, ye, t_mev):
    pf = electron_fermi_momentum_mev(rho, ye)
    return np.sqrt(pf * pf + M_E_MEV * M_E_MEV)


def eta_e(rho, ye, t_mev):
    return electron_chemical_potential_mev(rho, ye, t_mev) / np.maximum(t_mev, 1e-300)


def urca_correction(rho, ye, t_mev, flavor):
    eta = eta_e(rho, ye, t_mev)
    shift = 0.25 * np.tanh((eta - 1.0) / 4.0)
    if flavor == "nu_e":
        factor = 1.0 + shift
    elif flavor == "anti_nu_e":
        factor = 1.0 - shift
    else:
        factor = np.ones_like(eta)
    return np.clip(factor, 0.5, 1.5)


def main() -> None:
    PLOT_DIR.mkdir(parents=True, exist_ok=True)
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    rho = np.logspace(6, 13, 240)
    temp = np.logspace(np.log10(0.5), np.log10(30.0), 220)
    rr, tt = np.meshgrid(rho, temp)
    ye = 0.25
    eta = eta_e(rr, ye, tt)

    plt.figure(figsize=(6.5, 4.6))
    im = plt.pcolormesh(rho, temp, np.log10(np.maximum(eta, 1e-30)), shading="auto", cmap="viridis")
    plt.xscale("log")
    plt.yscale("log")
    plt.xlabel(r"$\rho$ [g cm$^{-3}$]")
    plt.ylabel(r"$T$ [MeV]")
    plt.title(r"Diagnostic electron degeneracy, $\eta_e=\mu_e/T$, $Y_e=0.25$")
    plt.colorbar(im, label=r"$\log_{10}\eta_e$")
    plt.tight_layout()
    plt.savefig(PLOT_DIR / "mev_electron_degeneracy_map.png", dpi=240)
    plt.close()

    plt.figure(figsize=(6.5, 4.4))
    for t in [1.0, 3.0, 10.0, 30.0]:
        plt.loglog(rho, eta_e(rho, ye, t), label=f"T={t:g} MeV")
    plt.xlabel(r"$\rho$ [g cm$^{-3}$]")
    plt.ylabel(r"$\eta_e$")
    plt.title("Electron degeneracy increases with density and decreases with temperature")
    plt.legend(frameon=False, fontsize=8)
    plt.tight_layout()
    plt.savefig(PLOT_DIR / "mev_eta_vs_density_temperature.png", dpi=240)
    plt.close()

    eta_rho_low = float(eta_e(1e8, ye, 10.0))
    eta_rho_high = float(eta_e(1e12, ye, 10.0))
    eta_t_low = float(eta_e(1e11, ye, 3.0))
    eta_t_high = float(eta_e(1e11, ye, 30.0))
    corr_nue = urca_correction(rr, ye, tt, "nu_e")
    corr_anue = urca_correction(rr, ye, tt, "anti_nu_e")

    report = OUT_DIR / "mev_degeneracy_validation.txt"
    report.write_text(
        "\n".join(
            [
                "MeV electron degeneracy diagnostic validation",
                "status: complete_for_diagnostic_approximation",
                "eta_definition: eta_e = sqrt(p_F^2 + m_e^2) / T",
                "p_F_definition: hbar*c*(3*pi^2*n_e)^(1/3)",
                "finite_temperature_mu_correction: ignored",
                f"eta_increases_with_rho: {eta_rho_high > eta_rho_low}",
                f"eta_decreases_with_T: {eta_t_low > eta_t_high}",
                f"eta_finite: {bool(np.all(np.isfinite(eta)))}",
                f"correction_finite: {bool(np.all(np.isfinite(corr_nue)) and np.all(np.isfinite(corr_anue)))}",
                f"correction_bounded_0p5_1p5: {bool(np.min(corr_nue) >= 0.5 and np.max(corr_nue) <= 1.5 and np.min(corr_anue) >= 0.5 and np.max(corr_anue) <= 1.5)}",
                f"eta_rho_1e8_T10: {eta_rho_low:.8e}",
                f"eta_rho_1e12_T10: {eta_rho_high:.8e}",
                f"eta_rho_1e11_T3: {eta_t_low:.8e}",
                f"eta_rho_1e11_T30: {eta_t_high:.8e}",
                "calibration_status: diagnostic_not_full_finite_temperature_chemical_potential",
                "",
            ]
        )
    )
    print(f"Saved: {report}")
    print(f"Saved: {PLOT_DIR / 'mev_electron_degeneracy_map.png'}")
    print(f"Saved: {PLOT_DIR / 'mev_eta_vs_density_temperature.png'}")


if __name__ == "__main__":
    main()
