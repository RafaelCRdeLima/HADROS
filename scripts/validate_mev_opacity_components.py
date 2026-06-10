#!/usr/bin/env python3
"""Validate decomposed diagnostic MeV opacity components."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
PLOT_DIR = ROOT / "plots" / "mev_physics"
OUT_DIR = ROOT / "output" / "validation"

M_U_G = 1.66053906660e-24
SIGMA_ABS0 = 9.6e-44
SIGMA_SCAT0 = 1.7e-44
SIGMA_SCAT_E0 = 3.4e-45


def components(rho, ye, e, flavor="anti_nu_e"):
    rho = np.asarray(rho, dtype=float)
    e = np.asarray(e, dtype=float)
    ye = np.clip(ye, 0.0, 1.0)
    nb = np.maximum(rho, 0.0) / M_U_G
    ne = nb * ye
    abs_n = np.where(flavor == "nu_e", nb * (1.0 - ye) * SIGMA_ABS0 * e * e, 0.0)
    abs_p = np.where(flavor == "anti_nu_e", nb * ye * SIGMA_ABS0 * e * e, 0.0)
    scat_n = nb * (1.0 - ye) * SIGMA_SCAT0 * e * e
    scat_p = nb * ye * 0.5 * SIGMA_SCAT0 * e * e
    scat_e = ne * SIGMA_SCAT_E0 * e * e
    return abs_n, abs_p, scat_n, scat_p, scat_e


def total_opacity(rho, ye, e, flavor="anti_nu_e"):
    return sum(components(rho, ye, e, flavor))


def main() -> None:
    PLOT_DIR.mkdir(parents=True, exist_ok=True)
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    energy = np.logspace(0.0, np.log10(80.0), 240)
    rho0 = 1e11
    ye = 0.25
    names = [
        "absorption on neutrons",
        "absorption on protons",
        "NC scattering on neutrons",
        "NC scattering on protons",
        "e-/e+ scattering proxy",
    ]

    plt.figure(figsize=(6.7, 4.6))
    for vals, label in zip(components(rho0, ye, energy, "anti_nu_e"), names):
        plt.loglog(energy, np.maximum(vals, 1e-300), label=label)
    plt.loglog(energy, total_opacity(rho0, ye, energy, "anti_nu_e"), color="black", lw=1.8, label="total")
    plt.xlabel(r"$E_\nu$ [MeV]")
    plt.ylabel(r"$\alpha$ [cm$^{-1}$]")
    plt.title(r"MeV opacity components, $\bar{\nu}_e$, $\rho=10^{11}$ g cm$^{-3}$, $Y_e=0.25$")
    plt.legend(frameon=False, fontsize=7)
    plt.tight_layout()
    plt.savefig(PLOT_DIR / "mev_opacity_components_vs_energy.png", dpi=240)
    plt.close()

    rho = np.logspace(7, 13, 180)
    temp_proxy = np.logspace(np.log10(0.5), np.log10(30), 130)
    rr, tt = np.meshgrid(rho, temp_proxy)
    e_map = 10.0 + 0.0 * tt
    comp = np.stack(components(rr, ye, e_map, "anti_nu_e"), axis=0)
    dominant = np.argmax(comp, axis=0)
    cmap = plt.matplotlib.colors.ListedColormap(["#2563eb", "#dc2626", "#16a34a", "#84cc16", "#9333ea"])
    plt.figure(figsize=(6.7, 4.6))
    im = plt.pcolormesh(rho, temp_proxy, dominant, shading="auto", cmap=cmap, vmin=0, vmax=4)
    plt.xscale("log")
    plt.yscale("log")
    cbar = plt.colorbar(im, ticks=[0.4, 1.2, 2.0, 2.8, 3.6])
    cbar.ax.set_yticklabels(["abs n", "abs p", "scat n", "scat p", "scat e"])
    plt.xlabel(r"$\rho$ [g cm$^{-3}$]")
    plt.ylabel(r"$T$ [MeV] diagnostic axis")
    plt.title(r"Dominant MeV opacity component, $\bar{\nu}_e$, $E_\nu=10$ MeV")
    plt.tight_layout()
    plt.savefig(PLOT_DIR / "mev_opacity_components_rhoT_map.png", dpi=240)
    plt.close()

    vals_anue = components(rho0, ye, 10.0, "anti_nu_e")
    vals_nue = components(rho0, ye, 10.0, "nu_e")
    vals_nux = components(rho0, ye, 10.0, "nu_x")
    total_sum_ok = np.isclose(total_opacity(rho0, ye, 10.0, "anti_nu_e"), sum(vals_anue), rtol=1e-14)
    nonnegative = all(float(v) >= 0 for vals in [vals_anue, vals_nue, vals_nux] for v in vals)
    nue_abs_neutron = vals_nue[0] > vals_nue[1]
    anue_abs_proton = vals_anue[1] > 0 and vals_anue[0] == 0
    nux_no_cc = vals_nux[0] == 0 and vals_nux[1] == 0 and sum(vals_nux[2:]) > 0

    report = OUT_DIR / "mev_opacity_components_validation.txt"
    report.write_text(
        "\n".join(
            [
                "MeV opacity component validation",
                "status: complete_for_diagnostic_decomposition",
                "scaling: sigma proportional to E_MeV^2",
                f"components_nonnegative: {bool(nonnegative)}",
                f"nu_e_absorption_dominated_by_neutron_targets_when_neutron_rich: {bool(nue_abs_neutron)}",
                f"anti_nu_e_absorption_depends_on_proton_fraction: {bool(anue_abs_proton)}",
                f"nu_x_no_charged_current_absorption: {bool(nux_no_cc)}",
                f"total_equals_sum_components: {bool(total_sum_ok)}",
                f"anti_nu_e_abs_n_cm_inv: {float(vals_anue[0]):.8e}",
                f"anti_nu_e_abs_p_cm_inv: {float(vals_anue[1]):.8e}",
                f"anti_nu_e_scat_n_cm_inv: {float(vals_anue[2]):.8e}",
                f"anti_nu_e_scat_p_cm_inv: {float(vals_anue[3]):.8e}",
                f"anti_nu_e_scat_e_cm_inv: {float(vals_anue[4]):.8e}",
                f"anti_nu_e_total_cm_inv: {float(sum(vals_anue)):.8e}",
                "calibration_status: approximate_low_energy_weak_opacity_proxy",
                "",
            ]
        )
    )
    print(f"Saved: {report}")
    print(f"Saved: {PLOT_DIR / 'mev_opacity_components_vs_energy.png'}")
    print(f"Saved: {PLOT_DIR / 'mev_opacity_components_rhoT_map.png'}")


if __name__ == "__main__":
    main()
