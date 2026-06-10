#!/usr/bin/env python3
"""Validate the post-processing MeV neutrino physics module."""

from __future__ import annotations

import math
import os
import shutil
import subprocess
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
PLOT_DIR = ROOT / "plots" / "mev_physics"
OUT_DIR = ROOT / "output" / "validation"
IMG_DIR = ROOT / "output" / "images"
CACHE = ROOT / os.environ.get("SMALL_CACHE_PATH", "output/rays/kerr_geodesics_small.bin")

M_U_G = 1.66053906660e-24
SIGMA_ABS0 = 9.6e-44
SIGMA_SCAT0 = 1.7e-44
SIGMA_SCAT_E0 = 3.4e-45


def thermal_shape(e_mev: np.ndarray | float, t_mev: np.ndarray | float):
    e = np.asarray(e_mev, dtype=float)
    t = np.asarray(t_mev, dtype=float)
    x = np.divide(e, t, out=np.full_like(e + t, np.inf), where=t > 0)
    shape = np.where(x > 120.0, 0.0, e * e / (np.exp(np.minimum(x, 120.0)) + 1.0))
    return shape / np.maximum(t**3, 1.0e-300)


def emissivity_channels(rho, t, ye, e, flavor="anti_nu_e"):
    ye = np.clip(ye, 0.0, 1.0)
    rho10 = rho / 1.0e10
    t10 = t / 10.0
    spec = thermal_shape(e, t)
    if flavor == "nu_x":
        urca = np.zeros_like(np.asarray(t, dtype=float))
    else:
        flavor_factor = ye if flavor == "nu_e" else 1.0 - ye
        urca = 1.0e30 * np.maximum(rho10, 0.0) * np.maximum(t10, 0.0) ** 6 * flavor_factor * spec
    pair_factor = 1.0 if flavor == "nu_x" else 0.7
    brems_factor = 1.0 if flavor == "nu_x" else 0.8
    pair = 3.0e28 * np.maximum(t10, 0.0) ** 9 * pair_factor * spec
    brems = 1.0e27 * np.maximum(rho10, 0.0) ** 2 * np.maximum(t10, 0.0) ** 5 * brems_factor * spec
    return urca, pair, brems


def opacities(rho, ye, e, flavor="anti_nu_e", absorption=True, scattering=True):
    ye = np.clip(ye, 0.0, 1.0)
    nb = rho / M_U_G
    if not absorption or flavor == "nu_x":
        alpha_abs = np.zeros_like(np.asarray(e, dtype=float))
    else:
        target = 1.0 - ye if flavor == "nu_e" else ye
        alpha_abs = nb * SIGMA_ABS0 * np.asarray(e, dtype=float) ** 2 * target
    if not scattering:
        alpha_scat = np.zeros_like(np.asarray(e, dtype=float))
    else:
        nucleon = nb * SIGMA_SCAT0 * np.asarray(e, dtype=float) ** 2 * ((1.0 - ye) + 0.5 * ye)
        electron = nb * ye * SIGMA_SCAT_E0 * np.asarray(e, dtype=float) ** 2
        alpha_scat = nucleon + electron
    return alpha_abs, alpha_scat


def transfer_update(i_in, j, alpha, ds):
    dtau = alpha * ds
    if alpha > 0.0 and dtau > 1.0e-12:
        attenuation = math.exp(-min(dtau, 700.0))
        return i_in * attenuation + (j / alpha) * (1.0 - attenuation)
    return i_in + j * ds


def sci_tag(value: float) -> str:
    text = f"{value:.0e}"
    return text.replace("+", "_").replace("-", "m").replace(".", "p")


def compact_tag(value: float) -> str:
    mantissa, exp = f"{value:.0e}".split("e")
    return f"{mantissa}e{int(exp)}"


def expected_image_path(rho0, enu, mev_enu, mev_norm, cam_theta, profile):
    profile_part = "" if profile == "gaussian" else f"_Profile_{profile}"
    return (
        IMG_DIR
        / (
            "kerr_image_cuda_cache_GBW"
            f"_rho0_torus_{compact_tag(rho0)}"
            f"{profile_part}"
            f"_Enu_{sci_tag(enu)}"
            f"_MeVEnu_{sci_tag(mev_enu)}"
            f"_MeVNorm_{sci_tag(mev_norm)}"
            f"_CamTheta_{cam_theta:.1f}".replace(".", "p")
            + ".dat"
        )
    )


def run_image(mev_model: str, dst: Path):
    enu = float(os.environ.get("SMALL_CACHE_ENU", "1e11"))
    rho0 = float(os.environ.get("SMALL_CACHE_TORUS_RHO0", "1.0"))
    mev_enu = float(os.environ.get("MEV_ENU", os.environ.get("MEV_ENERGY_MEV", "10.0")))
    mev_norm = float(os.environ.get("MEV_NORM", "1.0"))
    cam_theta = float(os.environ.get("CAM_THETA_DEG", "80.0"))
    profile = os.environ.get("DENSITY_PROFILE", "gaussian")

    cmd = [
        "./compute_kerr_image_from_cache",
        str(enu),
        os.environ.get("ASPIN", "0.0001"),
        os.environ.get("MBH_MSUN", "3.0"),
        str(rho0),
        os.environ.get("TORUS_R0_RG", "10.0"),
        os.environ.get("TORUS_SIGMA_RG", "5.0"),
        os.environ.get("TORUS_H_OVER_R", "0.25"),
        os.environ.get("SOURCE_R_RG", "3.5"),
        os.environ.get("SOURCE_SIGMA_RG", "1.0"),
        os.environ.get("SOURCE_THETA_DEG", "15.0"),
        os.environ.get("SOURCE_POWERLAW", "2.0"),
        os.environ.get("SOURCE_EMAX_GEV", "1.0e12"),
        os.environ.get("SOURCE_NORM", "1.0"),
        str(mev_enu),
        str(mev_norm),
        str(cam_theta),
        "data/sigma/sigma_nuN_CC_GBW.dat",
        profile,
        os.environ.get("TORUS_RADIAL_POWER", "2.0"),
        os.environ.get("FUNNEL_DEPLETION", "0.0"),
        os.environ.get("FUNNEL_THETA_DEG", "15.0"),
        os.environ.get("ENVELOPE_RHO0", "0.0"),
        os.environ.get("ENVELOPE_ALPHA", "2.5"),
        os.environ.get("TORUS_R_MIN_RG", "4.0"),
        os.environ.get("TORUS_R_MAX_RG", "60.0"),
        os.environ.get("RHO_FLOOR", "1.0e-99"),
        os.environ.get("CAM_R_OBS_RG", "60.0"),
        os.environ.get("USE_F3", "1"),
        str(CACHE),
        os.environ.get("SOURCE_MODEL", "inner_ring"),
        os.environ.get("SOURCE_FUNNEL_THETA_DEG", "20.0"),
        os.environ.get("SOURCE_DENSITY_Q", "1.0"),
        os.environ.get("SOURCE_RADIAL_S", "2.0"),
        os.environ.get("SOURCE_GRADIENT_DR_RG", "0.1"),
        os.environ.get("SOURCE_GRADIENT_DTHETA_DEG", "1.0"),
        os.environ.get("SOURCE_RHO_REF", str(rho0)),
        os.environ.get("SOURCE_CUTOFF_MIN", "0.0"),
        os.environ.get("SOURCE_CUTOFF_MAX", "1.0e2"),
        os.environ.get("SPECTRAL_MODEL", "monochromatic"),
        os.environ.get("SPECTRAL_GAMMA", "2.0"),
        os.environ.get("SPECTRAL_ECUT_GEV", "1.0e12"),
        os.environ.get("SPECTRAL_E_MIN_GEV", "1.0e5"),
        os.environ.get("SPECTRAL_E_MAX_GEV", "1.0e12"),
        os.environ.get("SPECTRAL_N_BINS", "8"),
        mev_model,
        os.environ.get("MEV_FLAVOR", "anti_nu_e"),
        os.environ.get("MEV_INCLUDE_URCA", "1"),
        os.environ.get("MEV_INCLUDE_PAIR", "1"),
        os.environ.get("MEV_INCLUDE_BREMS", "1"),
        os.environ.get("MEV_INCLUDE_ABSORPTION", "1"),
        os.environ.get("MEV_INCLUDE_SCATTERING", "1"),
        os.environ.get("MEV_THERMAL_PROFILE", "inner_hot_torus"),
        os.environ.get("MEV_YE_PROFILE", "neutron_rich_torus"),
        os.environ.get("MEV_T0_MEV", "6.0"),
        os.environ.get("MEV_T_FLOOR_MEV", "0.1"),
        os.environ.get("MEV_T_POWER", "0.2"),
        os.environ.get("MEV_YE_TORUS", "0.25"),
        os.environ.get("MEV_YE_FUNNEL", "0.55"),
        os.environ.get("MEV_YE_ENVELOPE", "0.45"),
        os.environ.get("MEV_YE_FLOOR", "0.01"),
        os.environ.get("MEV_YE_CEIL", "0.60"),
        os.environ.get("MEV_SPECTRAL_MODE", "monochromatic"),
        os.environ.get("MEV_E_MIN_MEV", "3.0"),
        os.environ.get("MEV_E_MAX_MEV", "50.0"),
        os.environ.get("MEV_N_BINS", "8"),
        os.environ.get("MEV_USE_DEGENERACY_CORRECTION", "0"),
        os.environ.get("MEV_INCLUDE_ABS_N", "1"),
        os.environ.get("MEV_INCLUDE_ABS_P", "1"),
        os.environ.get("MEV_INCLUDE_SCAT_N", "1"),
        os.environ.get("MEV_INCLUDE_SCAT_P", "1"),
        os.environ.get("MEV_INCLUDE_SCAT_E", "1"),
    ]
    subprocess.run(cmd, cwd=ROOT, check=True)
    src = expected_image_path(rho0, enu, mev_enu, mev_norm, cam_theta, profile)
    if not src.exists():
        raise FileNotFoundError(src)
    shutil.copyfile(src, dst)


def load_image_grid(path: Path, column: int):
    data = np.loadtxt(path, comments="#")
    if data.ndim == 1:
        data = data[None, :]
    nx = int(np.max(data[:, 0])) + 1
    ny = int(np.max(data[:, 1])) + 1
    grid = np.zeros((ny, nx))
    for row in data:
        i = int(row[0])
        j = int(row[1])
        grid[ny - 1 - j, i] = row[column]
    return grid, data


def plot_channel_curves():
    temp = np.logspace(np.log10(0.5), np.log10(30.0), 240)
    urca, pair, brems = emissivity_channels(1.0e11, temp, 0.25, 10.0)
    total = urca + pair + brems
    plt.figure(figsize=(6.2, 4.2))
    plt.loglog(temp, urca, label="URCA / beta processes")
    plt.loglog(temp, pair, label="pair annihilation")
    plt.loglog(temp, brems, label="NN bremsstrahlung")
    plt.loglog(temp, total, color="black", linewidth=1.8, label="total")
    plt.xlabel(r"$T$ [MeV]")
    plt.ylabel("local spectral emissivity proxy")
    plt.title(r"MeV emissivity channels, $\rho=10^{11}\,\mathrm{g\,cm^{-3}}$, $E_\nu=10$ MeV")
    plt.legend(frameon=False, fontsize=8)
    plt.tight_layout()
    plt.savefig(PLOT_DIR / "mev_emissivity_channels.png", dpi=240)
    plt.close()


def plot_opacity():
    energies = np.logspace(0.0, np.log10(80.0), 220)
    alpha_abs, alpha_scat = opacities(1.0e11, 0.25, energies)
    plt.figure(figsize=(6.2, 4.2))
    plt.loglog(energies, alpha_abs, label="absorption")
    plt.loglog(energies, alpha_scat, label="scattering")
    plt.loglog(energies, alpha_abs + alpha_scat, color="black", label="total")
    plt.xlabel(r"$E_\nu$ [MeV]")
    plt.ylabel(r"$\alpha$ [cm$^{-1}$]")
    plt.title(r"MeV opacity, $\rho=10^{11}\,\mathrm{g\,cm^{-3}}$, $Y_e=0.25$")
    plt.legend(frameon=False, fontsize=8)
    plt.tight_layout()
    plt.savefig(PLOT_DIR / "mev_opacity_vs_energy.png", dpi=240)
    plt.close()


def plot_dominance_map():
    rho = np.logspace(6.0, 13.0, 160)
    temp = np.logspace(np.log10(0.5), np.log10(30.0), 140)
    rr, tt = np.meshgrid(rho, temp)
    urca, pair, brems = emissivity_channels(rr, tt, 0.25, 10.0)
    dominant = np.argmax(np.stack([urca, pair, brems]), axis=0)
    plt.figure(figsize=(6.4, 4.5))
    cmap = plt.matplotlib.colors.ListedColormap(["#355caa", "#2a9d8f", "#d1495b"])
    plt.pcolormesh(rho, temp, dominant, shading="auto", cmap=cmap, vmin=0, vmax=2)
    plt.xscale("log")
    plt.yscale("log")
    cbar = plt.colorbar(ticks=[0.33, 1.0, 1.67])
    cbar.ax.set_yticklabels(["URCA", "pair", "brems"])
    plt.xlabel(r"$\rho$ [g cm$^{-3}$]")
    plt.ylabel(r"$T$ [MeV]")
    plt.title(r"Dominant MeV emissivity channel, $Y_e=0.25$, $E_\nu=10$ MeV")
    plt.tight_layout()
    plt.savefig(PLOT_DIR / "mev_channel_dominance_map.png", dpi=240)
    plt.close()


def plot_flavor_comparison():
    temp = np.logspace(np.log10(0.5), np.log10(30.0), 220)
    rho = 1.0e11
    ye = 0.25
    e_mev = 10.0
    plt.figure(figsize=(6.4, 4.3))
    for flavor, color in [
        ("nu_e", "#2563eb"),
        ("anti_nu_e", "#dc2626"),
        ("nu_x", "#16a34a"),
    ]:
        urca, pair, brems = emissivity_channels(rho, temp, ye, e_mev, flavor=flavor)
        plt.loglog(temp, urca + pair + brems, color=color, label=f"{flavor} emissivity")
    plt.xlabel(r"$T$ [MeV]")
    plt.ylabel("local spectral emissivity proxy")
    plt.title(r"MeV flavor comparison, $\rho=10^{11}\,\mathrm{g\,cm^{-3}}$, $Y_e=0.25$")
    plt.legend(frameon=False, fontsize=8)
    plt.tight_layout()
    plt.savefig(PLOT_DIR / "mev_flavor_comparison.png", dpi=240)
    plt.close()


def plot_transfer_limits():
    j = 2.0
    alpha = 3.0
    i_in = 0.7
    dtau = np.logspace(-8.0, 3.0, 260)
    i_out = np.array([transfer_update(i_in, j, alpha, x / alpha) for x in dtau])
    thin = i_in + j * dtau / alpha
    source = np.full_like(dtau, j / alpha)
    plt.figure(figsize=(6.4, 4.3))
    plt.loglog(dtau, np.abs(i_out - thin) / np.maximum(np.abs(thin), 1.0e-300), label="thin-limit relative error")
    plt.loglog(dtau, np.abs(i_out - source) / np.maximum(np.abs(source), 1.0e-300), label="source-function relative error")
    plt.axvspan(1.0e-8, 1.0e-4, color="#dbeafe", alpha=0.5, label="optically thin")
    plt.axvspan(30.0, 1.0e3, color="#fee2e2", alpha=0.5, label="optically thick")
    plt.xlabel(r"$\Delta\tau=\alpha ds$")
    plt.ylabel("relative error")
    plt.title("MeV radiative-transfer limits")
    plt.legend(frameon=False, fontsize=8)
    plt.tight_layout()
    plt.savefig(PLOT_DIR / "mev_transfer_limits.png", dpi=240)
    plt.close()


def plot_image_comparison(toy_path: Path, physical_path: Path):
    toy_grid, toy_data = load_image_grid(toy_path, 10)
    phy_grid, phy_data = load_image_grid(physical_path, 10)
    scale = max(float(np.nanmax(toy_grid)), float(np.nanmax(phy_grid)), 1.0e-300)
    fig, axes = plt.subplots(1, 2, figsize=(8.4, 3.8), constrained_layout=True)
    for ax, grid, title in [
        (axes[0], toy_grid, "toy MeV emissivity"),
        (axes[1], phy_grid, "physical MeV emissivity"),
    ]:
        im = ax.imshow(np.log10(np.maximum(grid / scale, 1.0e-12)), origin="lower", cmap="magma")
        ax.set_title(title)
        ax.set_xlabel("camera pixel")
        ax.set_ylabel("camera pixel")
    fig.colorbar(im, ax=axes, label=r"$\log_{10}(I_{\rm MeV}/I_{\max})$")
    plt.savefig(PLOT_DIR / "mev_image_physical_vs_toy.png", dpi=240)
    plt.close()
    return toy_data, phy_data


def main():
    PLOT_DIR.mkdir(parents=True, exist_ok=True)
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    IMG_DIR.mkdir(parents=True, exist_ok=True)

    plot_channel_curves()
    plot_opacity()
    plot_dominance_map()
    plot_flavor_comparison()
    plot_transfer_limits()

    toy_copy = OUT_DIR / "mev_toy_image_validation.dat"
    physical_copy = OUT_DIR / "mev_physical_image_validation.dat"
    run_image("toy", toy_copy)
    run_image("physical", physical_copy)
    toy_data, phy_data = plot_image_comparison(toy_copy, physical_copy)

    rho_test = 1.0e12
    t_test = 10.0
    ye_test = 0.25
    e_test = 10.0
    urca, pair, brems = emissivity_channels(rho_test, t_test, ye_test, e_test)
    alpha_abs, alpha_scat = opacities(rho_test, ye_test, np.array([e_test]))
    alpha = float(alpha_abs[0] + alpha_scat[0])
    j = float(urca + pair + brems)
    ds = 1.0e-8 / max(alpha, 1.0e-300)
    thin_exact = (j / alpha) * (1.0 - math.exp(-alpha * ds)) if alpha > 0 else j * ds
    thin_limit = j * ds
    alpha_zero_ok = math.isclose(thin_exact, thin_limit, rel_tol=1.0e-6, abs_tol=0.0)

    finite_images = np.all(np.isfinite(toy_data)) and np.all(np.isfinite(phy_data))
    positive_emissivity = j > 0.0 and urca > 0.0 and pair > 0.0 and brems > 0.0
    positive_opacity = alpha > 0.0
    physical_diff = not np.allclose(toy_data[:, 10], phy_data[:, 10], rtol=1.0e-6, atol=0.0)
    urca_dominates = urca > pair and urca > brems
    pair_5 = emissivity_channels(1.0e8, 5.0, ye_test, e_test)[1]
    pair_10 = emissivity_channels(1.0e8, 10.0, ye_test, e_test)[1]
    pair_ratio = float(pair_10 / max(pair_5, 1.0e-300))
    pair_steep = pair_ratio > 20.0

    brems_rho1 = emissivity_channels(1.0e10, t_test, ye_test, e_test)[2]
    brems_rho2 = emissivity_channels(2.0e10, t_test, ye_test, e_test)[2]
    brems_ratio = float(brems_rho2 / max(brems_rho1, 1.0e-300))
    brems_rho2_scaling_ok = math.isclose(brems_ratio, 4.0, rel_tol=5.0e-2)

    urca_nue = float(emissivity_channels(rho_test, t_test, ye_test, e_test, "nu_e")[0])
    urca_anue = float(emissivity_channels(rho_test, t_test, ye_test, e_test, "anti_nu_e")[0])
    urca_nux = float(emissivity_channels(rho_test, t_test, ye_test, e_test, "nu_x")[0])
    flavor_ye_behavior = urca_anue > urca_nue and urca_nux == 0.0
    nux_abs, nux_scat = opacities(rho_test, ye_test, np.array([e_test]), flavor="nu_x")
    nux_no_absorption = float(nux_abs[0]) == 0.0 and float(nux_scat[0]) > 0.0

    i_in = 0.7
    ds_thin = 1.0e-9 / max(alpha, 1.0e-300)
    thin_transfer = transfer_update(i_in, j, alpha, ds_thin)
    thin_expected = i_in + j * ds_thin
    optically_thin_limit_ok = math.isclose(thin_transfer, thin_expected, rel_tol=1.0e-6, abs_tol=1.0e-8)
    ds_thick = 100.0 / max(alpha, 1.0e-300)
    thick_transfer = transfer_update(i_in, j, alpha, ds_thick)
    thick_expected = j / alpha
    optically_thick_limit_ok = math.isclose(thick_transfer, thick_expected, rel_tol=1.0e-10, abs_tol=1.0e-10)

    report = OUT_DIR / "mev_physics_validation.txt"
    report.write_text(
        "\n".join(
            [
                "MeV neutrino physics validation",
                "status: complete_for_cpu_postprocessing",
                "MEV_CPU_MODEL: physical",
                "MEV_CUDA_STATUS: legacy_toy_not_equivalent",
                "model_default: physical",
                "flavor_default: anti_nu_e",
                "emissivity_units: local_spectral_emissivity_proxy_per_MeV",
                "opacity_units: cm^-1",
                f"finite_outputs: {finite_images}",
                f"positive_emissivity: {positive_emissivity}",
                f"positive_opacity: {positive_opacity}",
                f"alpha_zero_limit_ok: {alpha_zero_ok}",
                f"optically_thin_limit_ok: {optically_thin_limit_ok}",
                f"optically_thick_source_function_limit_ok: {optically_thick_limit_ok}",
                f"physical_model_differs_from_toy_image: {physical_diff}",
                f"urca_dominates_dense_hot_electron_flavor: {urca_dominates}",
                f"pair_production_steep_temperature_growth: {pair_steep}",
                f"pair_T10_over_T5_ratio: {pair_ratio:.8e}",
                f"brems_rho2_scaling_ok: {brems_rho2_scaling_ok}",
                f"brems_rho2_over_rho1_ratio: {brems_ratio:.8e}",
                f"flavor_Ye_behavior_ok: {flavor_ye_behavior}",
                f"nu_x_no_charged_current_absorption_ok: {nux_no_absorption}",
                f"test_urca_nu_e_emissivity: {urca_nue:.8e}",
                f"test_urca_anti_nu_e_emissivity: {urca_anue:.8e}",
                f"test_urca_nu_x_emissivity: {urca_nux:.8e}",
                f"test_urca_emissivity: {float(urca):.8e}",
                f"test_pair_emissivity: {float(pair):.8e}",
                f"test_brems_emissivity: {float(brems):.8e}",
                f"test_absorption_opacity_cm_inv: {float(alpha_abs[0]):.8e}",
                f"test_scattering_opacity_cm_inv: {float(alpha_scat[0]):.8e}",
                f"test_nu_x_absorption_opacity_cm_inv: {float(nux_abs[0]):.8e}",
                f"test_nu_x_scattering_opacity_cm_inv: {float(nux_scat[0]):.8e}",
                f"toy_image_total_I_mev: {float(np.sum(toy_data[:, 10])):.8e}",
                f"physical_image_total_I_mev: {float(np.sum(phy_data[:, 10])):.8e}",
                "cuda_status: CUDA MeV path remains legacy/toy and is not equivalent to the canonical CPU physical MeV module.",
                "",
            ]
        )
    )
    print(f"Saved: {report}")
    for name in [
        "mev_emissivity_channels.png",
        "mev_opacity_vs_energy.png",
        "mev_channel_dominance_map.png",
        "mev_flavor_comparison.png",
        "mev_transfer_limits.png",
        "mev_image_physical_vs_toy.png",
    ]:
        print(f"Saved: {PLOT_DIR / name}")


if __name__ == "__main__":
    main()
