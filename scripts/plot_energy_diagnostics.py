import glob
import os
import re
import struct
import sys
from pathlib import Path

import matplotlib
import numpy as np

matplotlib.use("Agg")

import matplotlib.pyplot as plt

KIMG_MAGIC = 0x4B494D47
KIMG_HEADER = struct.Struct("<7i4d")
NBINS = 80
PI = np.pi
C_CGS = 2.99792458e10
G_CGS = 6.67430e-8
MSUN_G = 1.98847e33
M_U_G = 1.66053906660e-24

ASPIN = float(os.environ.get("ASPIN", "0.0001"))
MBH_MSUN = float(os.environ.get("MBH_MSUN", "3.0"))
TORUS_R0_RG = float(os.environ.get("TORUS_R0_RG", "10.0"))
TORUS_SIGMA_RG = float(os.environ.get("TORUS_SIGMA_RG", "5.0"))
TORUS_H_OVER_R = float(os.environ.get("TORUS_H_OVER_R", "0.25"))
MEV_ENU = float(os.environ.get("MEV_ENU", "10.0"))
SIGMA_TABLE_PATH = Path(
    os.environ.get("SIGMA_TABLE", "data/sigma/sigma_nuN_CC_GBW.dat")
)
MEV_SIGMA_ABS0_CM2 = 9.6e-44
MEV_SIGMA_SCAT0_CM2 = 1.7e-44


def sci_latex(value, precision=1):
    mantissa, exponent = f"{value:.{precision}e}".split("e")
    exponent = int(exponent)
    return rf"{mantissa}\times 10^{{{exponent}}}"


def sci_latex_compact(value, precision=1):
    mantissa_text, exponent = f"{value:.{precision}e}".split("e")
    mantissa = float(mantissa_text)
    exponent = int(exponent)

    if np.isclose(mantissa, 1.0):
        return rf"10^{{{exponent}}}"

    return rf"{mantissa_text}\times 10^{{{exponent}}}"


def rg_cm(m_bh_msun):
    return G_CGS * m_bh_msun * MSUN_G / (C_CGS * C_CGS)


def rho_folder_tag(value):
    mantissa, exponent = f"{value:.0e}".split("e")
    return f"rho0_torus_{mantissa}e{int(exponent)}"


def rho_value_from_tag(tag):
    match = re.fullmatch(r"rho0_torus_([0-9]+(?:p[0-9]+)?)e(-?[0-9]+)", tag)
    if not match:
        raise ValueError(f"Invalid torus density tag: {tag}")

    mantissa = float(match.group(1).replace("p", "."))
    exponent = int(match.group(2))
    return mantissa * 10.0**exponent


def title_suffix(torus_rho0):
    return (
        rf"$\rho_{{0,\rm torus}}="
        rf"{sci_latex(torus_rho0)}"
        rf"\,\mathrm{{g\,cm^{{-3}}}}$"
    )


def load_kerr_image_binary(path):
    with open(path, "rb") as f:
        header_bytes = f.read(KIMG_HEADER.size)
        if len(header_bytes) != KIMG_HEADER.size:
            raise ValueError(f"Invalid binary image header: {path}")

        (
            magic,
            version,
            backend,
            nx,
            ny,
            ncols,
            reserved,
            enu,
            mev_enu,
            mev_norm,
            cam_theta,
        ) = KIMG_HEADER.unpack(header_bytes)

        if magic != KIMG_MAGIC:
            raise ValueError(f"Invalid binary image magic: {path}")
        if version != 1:
            raise ValueError(f"Unsupported binary image version {version}: {path}")
        if ncols < 8:
            raise ValueError(f"Invalid binary image column count {ncols}: {path}")

        payload = np.fromfile(f, dtype="<f8")

    expected = nx * ny * ncols
    if payload.size != expected:
        raise ValueError(
            f"Invalid binary image payload size in {path}: "
            f"expected {expected}, found {payload.size}"
        )

    return payload.reshape((nx * ny, ncols))


def load_kerr_image(path):
    path = Path(path)
    if path.suffix == ".bin":
        return load_kerr_image_binary(path)

    return np.loadtxt(path)


def load_sigma_table(path):
    data = np.loadtxt(path, comments="#")
    if data.ndim != 2 or data.shape[1] < 3:
        raise ValueError(f"Invalid sigma table: {path}")

    return data[:, 0], data[:, 2]


def sigma_cm2_log_interp(energy_gev, sigma_energy, sigma_cm2):
    if energy_gev < sigma_energy[0] or energy_gev > sigma_energy[-1]:
        return 0.0

    return float(
        np.exp(
            np.interp(
                np.log(energy_gev),
                np.log(sigma_energy),
                np.log(sigma_cm2),
            )
        )
    )


def torus_in(r_rg, theta):
    delta = theta - 0.5 * PI
    return (r_rg > 4.0) and (r_rg < 18.0) and (abs(delta) < 0.45)


def torus_rho(r_rg, theta, rho0):
    if not torus_in(r_rg, theta):
        return 0.0

    delta = theta - 0.5 * PI
    radial = np.exp(-((r_rg - TORUS_R0_RG) / TORUS_SIGMA_RG) ** 2)
    vertical = np.exp(-((delta / TORUS_H_OVER_R) ** 2))
    return rho0 * radial * vertical


def torus_ye(r_rg, theta):
    if not torus_in(r_rg, theta):
        return 0.0

    return 0.2 + 0.1 * np.exp(-((r_rg - TORUS_R0_RG) / 4.0) ** 2)


def kerr_equatorial_radial_dl_dr(r_rg):
    a_bh = np.clip(ASPIN, 0.0, 0.999999)
    delta = r_rg * r_rg - 2.0 * r_rg + a_bh * a_bh

    return np.where(delta > 0.0, r_rg / np.sqrt(delta), np.inf)


def mev_opacity_cm_inv(rho, ye, enu_mev):
    if rho <= 0.0 or enu_mev <= 0.0:
        return 0.0

    nb = rho / M_U_G
    neutron_fraction = np.clip(1.0 - ye, 0.0, 1.0)
    proton_fraction = np.clip(ye, 0.0, 1.0)

    absorption_target = neutron_fraction + 0.25 * proton_fraction
    scattering_composition = 1.0 + 0.5 * np.clip(ye, 0.0, 1.0)

    kappa_abs = (
        nb
        * MEV_SIGMA_ABS0_CM2
        * enu_mev
        * enu_mev
        * absorption_target
    )
    kappa_scat = (
        nb
        * MEV_SIGMA_SCAT0_CM2
        * enu_mev
        * enu_mev
        * scattering_composition
    )

    return kappa_abs + kappa_scat


def find_legacy_files():
    return (
        glob.glob("output/images/kerr_image_cuda_cache_Enu_*.bin")
        + glob.glob("output/images/kerr_image_cuda_Enu_*.bin")
        + glob.glob("output/images/kerr_image_Enu_*.bin")
        + glob.glob("output/images/kerr_image_cuda_cache_Enu_*.dat")
        + glob.glob("output/images/kerr_image_cuda_Enu_*.dat")
        + glob.glob("output/images/kerr_image_Enu_*.dat")
    )


def find_rho_files(rho_tag):
    return (
        glob.glob(f"output/images/kerr_image_cuda_cache_{rho_tag}_Enu_*.bin")
        + glob.glob(f"output/images/kerr_image_cuda_{rho_tag}_Enu_*.bin")
        + glob.glob(f"output/images/kerr_image_{rho_tag}_Enu_*.bin")
        + glob.glob(f"output/images/kerr_image_cuda_cache_{rho_tag}_Enu_*.dat")
        + glob.glob(f"output/images/kerr_image_cuda_{rho_tag}_Enu_*.dat")
        + glob.glob(f"output/images/kerr_image_{rho_tag}_Enu_*.dat")
    )


def discover_rho_tags():
    tags = set()
    pattern = re.compile(r"(rho0_torus_[0-9]+(?:p[0-9]+)?e-?[0-9]+)_Enu_")

    for fname in glob.glob("output/images/kerr_image*_rho0_torus_*_Enu_*"):
        match = pattern.search(Path(fname).name)
        if match:
            tags.add(match.group(1))

    return sorted(tags, key=rho_value_from_tag)


def extract_energy(filename):
    match = re.search(r"Enu_([0-9]+e_[0-9]+)", filename)

    if not match:
        raise ValueError(f"Could not extract energy from {filename}")

    tag = match.group(1).replace("_", "+")
    return float(tag)


def file_priority(filename):
    path = Path(filename)
    is_cuda_cache = path.name.startswith("kerr_image_cuda_cache_")
    is_cuda = path.name.startswith("kerr_image_cuda_")
    is_bin = path.suffix == ".bin"

    if is_cuda_cache and is_bin:
        return 0
    if is_cuda and is_bin:
        return 1
    if is_bin:
        return 2
    if is_cuda_cache:
        return 3
    if is_cuda:
        return 4

    return 5


def select_files_by_energy(candidate_files):
    selected_by_energy = {}

    for fname in candidate_files:
        energy = extract_energy(fname)
        current = selected_by_energy.get(energy)

        if current is None or file_priority(fname) < file_priority(current):
            selected_by_energy[energy] = fname

    return [
        selected_by_energy[energy]
        for energy in sorted(selected_by_energy)
    ]


def make_channels():
    return {
        "UHE": {
            "tau_col": 4,
            "P_col": 5,
            "I_col": 6,
            "flux_total": [],
            "tau_I_weighted": [],
            "P_I_weighted": [],
            "contrast_shadow": [],
            "radial_profiles": [],
        },
        "MeV": {
            "tau_col": 8,
            "P_col": 9,
            "I_col": 10,
            "flux_total": [],
            "tau_I_weighted": [],
            "P_I_weighted": [],
            "contrast_shadow": [],
            "radial_profiles": [],
        },
    }


def compute_channel_diagnostics(data, b, channel):
    tau = data[:, channel["tau_col"]]
    p_surv = data[:, channel["P_col"]]
    i_obs = data[:, channel["I_col"]]

    flux = np.sum(i_obs)

    if flux > 0:
        tau_weighted = np.sum(i_obs * tau) / flux
        p_weighted = np.sum(i_obs * p_surv) / flux
    else:
        tau_weighted = np.nan
        p_weighted = np.nan

    bins = np.linspace(0.0, np.max(b), NBINS + 1)
    b_mid = 0.5 * (bins[:-1] + bins[1:])
    i_radial = np.zeros(NBINS)

    for k in range(NBINS):
        mask = (b >= bins[k]) & (b < bins[k + 1])

        if np.any(mask):
            i_radial[k] = np.mean(i_obs[mask])
        else:
            i_radial[k] = np.nan

    if np.nanmax(i_radial) > 0:
        i_radial_norm = i_radial / np.nanmax(i_radial)
    else:
        i_radial_norm = i_radial

    center_limit = np.percentile(b, 10)
    center_mask = b <= center_limit

    if np.any(center_mask):
        i_center = np.mean(i_obs[center_mask])
    else:
        i_center = np.nan

    i_ring = np.nanmax(i_radial)

    if np.isfinite(i_center) and i_center > 0:
        contrast = i_ring / i_center
    else:
        contrast = np.nan

    return {
        "flux_total": flux,
        "tau_I_weighted": tau_weighted,
        "P_I_weighted": p_weighted,
        "contrast_shadow": contrast,
        "radial_profile": (b_mid, i_radial_norm),
    }


def save_radial_tau_profile(energies, torus_rho0, rho_tag, plot_dir):
    sigma_energy, sigma_cm2 = load_sigma_table(SIGMA_TABLE_PATH)

    a_bh = np.clip(ASPIN, 0.0, 0.999999)
    r_horizon = 1.0 + np.sqrt(1.0 - a_bh * a_bh)
    r_outer = max(20.0, TORUS_R0_RG + 3.0 * TORUS_SIGMA_RG)
    r_start = r_horizon + 1.0e-3
    r_grid = np.linspace(r_start, r_outer, 4000)
    r_mid = 0.5 * (r_grid[:-1] + r_grid[1:])
    dr = np.diff(r_grid)
    dl_cm = kerr_equatorial_radial_dl_dr(r_mid) * dr * rg_cm(MBH_MSUN)

    theta_eq = 0.5 * PI
    rho = np.array([torus_rho(r, theta_eq, torus_rho0) for r in r_mid])
    ye = np.array([torus_ye(r, theta_eq) for r in r_mid])
    nb = rho / M_U_G

    mev_kappa = np.array(
        [
            mev_opacity_cm_inv(local_rho, local_ye, MEV_ENU)
            for local_rho, local_ye in zip(rho, ye)
        ]
    )
    tau_mev = np.concatenate(([0.0], np.cumsum(mev_kappa * dl_cm)))

    plt.figure(figsize=(7, 5))
    for energy in energies:
        sigma = sigma_cm2_log_interp(energy, sigma_energy, sigma_cm2)
        tau_uhe = np.concatenate(([0.0], np.cumsum(nb * sigma * dl_cm)))
        plt.semilogy(
            r_grid - r_horizon,
            np.maximum(tau_uhe, 1.0e-300),
            label=rf"$E_\nu={sci_latex_compact(energy)}\,\mathrm{{GeV}}$",
        )

    plt.axvline(4.0 - r_horizon, color="0.55", linestyle=":", linewidth=1.0)
    plt.axvline(18.0 - r_horizon, color="0.55", linestyle=":", linewidth=1.0)
    plt.xlabel(r"$(r-r_+)/r_g$")
    plt.ylabel(r"Cumulative optical depth $\tau(<r)$")
    plt.title(title_suffix(torus_rho0))
    plt.ylim(bottom=1e-4)
    plt.legend(
        fontsize=7,
        loc="center left",
        bbox_to_anchor="tight",
        borderaxespad=0.0,
    )
    plt.tight_layout(rect=(0.0, 0.0, 0.86, 1.0))

    uhe_output = plot_dir / "diagnostic_tau_radial_from_horizon_UHE.png"
    plt.savefig(uhe_output, dpi=250)
    plt.close()

    plt.figure(figsize=(7, 5))
    plt.semilogy(
        r_grid - r_horizon,
        np.maximum(tau_mev, 1.0e-300),
        color="black",
        linestyle="-",
        linewidth=2.0,
        label=rf"MeV $E_\nu={MEV_ENU:g}\,\mathrm{{MeV}}$",
    )

    plt.axvline(4.0 - r_horizon, color="0.55", linestyle=":", linewidth=1.0)
    plt.axvline(18.0 - r_horizon, color="0.55", linestyle=":", linewidth=1.0)
    plt.xlabel(r"$(r-r_+)/r_g$")
    plt.ylabel(r"Cumulative optical depth $\tau(<r)$")
    plt.title(title_suffix(torus_rho0))
    plt.ylim(bottom=-3)
    plt.legend(fontsize=7)
    plt.tight_layout()

    mev_output = plot_dir / "diagnostic_tau_radial_from_horizon_MeV.png"
    plt.savefig(mev_output, dpi=250)
    plt.close()

    print(f"Saved: {uhe_output}")
    print(f"Saved: {mev_output}")


def save_diagnostics(files, torus_rho0, rho_tag):
    if not files:
        print(f"Skipping {rho_tag}: no image files found")
        return

    files = select_files_by_energy(files)
    plot_dir = Path("plots") / rho_tag / "energy_diagnostics"
    plot_dir.mkdir(parents=True, exist_ok=True)

    energies = []
    channels = make_channels()

    for fname in files:
        energy = extract_energy(fname)
        data = load_kerr_image(fname)

        alpha = data[:, 2]
        beta = data[:, 3]
        b = np.sqrt(alpha**2 + beta**2)

        energies.append(energy)

        for channel in channels.values():
            if data.shape[1] <= channel["I_col"]:
                continue

            diagnostics = compute_channel_diagnostics(data, b, channel)

            channel["flux_total"].append(diagnostics["flux_total"])
            channel["tau_I_weighted"].append(diagnostics["tau_I_weighted"])
            channel["P_I_weighted"].append(diagnostics["P_I_weighted"])
            channel["contrast_shadow"].append(diagnostics["contrast_shadow"])
            channel["radial_profiles"].append(
                (
                    energy,
                    diagnostics["radial_profile"][0],
                    diagnostics["radial_profile"][1],
                )
            )

    order = np.argsort(energies)
    energies = np.array(energies)[order]

    for channel in channels.values():
        channel["flux_total"] = np.array(channel["flux_total"])[order]
        channel["tau_I_weighted"] = np.array(channel["tau_I_weighted"])[order]
        channel["P_I_weighted"] = np.array(channel["P_I_weighted"])[order]
        channel["contrast_shadow"] = np.array(channel["contrast_shadow"])[order]
        channel["radial_profiles"] = [
            channel["radial_profiles"][i]
            for i in order
        ]

    suffix = title_suffix(torus_rho0)
    save_radial_tau_profile(energies, torus_rho0, rho_tag, plot_dir)

    plt.figure(figsize=(7, 5))
    plt.loglog(energies, channels["UHE"]["flux_total"], marker="o", label="UHE")
    plt.loglog(energies, channels["MeV"]["flux_total"], marker="s", label="MeV")
    plt.xlabel(r"$E_\nu\ [\mathrm{GeV}]$")
    plt.ylabel(r"$\sum I_{\rm obs}$")
    plt.title(suffix)
    plt.legend()
    plt.tight_layout()
    plt.savefig(plot_dir / "diagnostic_total_flux_vs_energy.png", dpi=250)

    plt.figure(figsize=(7, 5))
    plt.loglog(
        energies,
        channels["UHE"]["tau_I_weighted"],
        marker="o",
        linestyle="-",
        label=r"UHE $\langle\tau\rangle_I$",
    )
    plt.loglog(
        energies,
        channels["MeV"]["tau_I_weighted"],
        marker="o",
        linestyle="--",
        label=r"MeV $\langle\tau\rangle_I$",
    )
    plt.xlabel(r"$E_\nu\ [\mathrm{GeV}]$")
    plt.ylabel(r"Optical depth")
    plt.title(suffix)
    plt.legend()
    plt.tight_layout()
    plt.savefig(plot_dir / "diagnostic_tau_vs_energy.png", dpi=250)

    plt.figure(figsize=(7, 5))
    plt.semilogx(
        energies,
        channels["UHE"]["P_I_weighted"],
        marker="o",
        linestyle="-",
        label=r"UHE $\langle P_{\rm surv}\rangle_I$",
    )
    plt.semilogx(
        energies,
        channels["MeV"]["P_I_weighted"],
        marker="o",
        linestyle="--",
        label=r"MeV $\langle P_{\rm surv}\rangle_I$",
    )
    plt.xlabel(r"$E_\nu\ [\mathrm{GeV}]$")
    plt.ylabel(r"Survival probability")
    plt.title(suffix)
    plt.legend()
    plt.tight_layout()
    plt.savefig(plot_dir / "diagnostic_survival_vs_energy.png", dpi=250)

    plt.figure(figsize=(7, 5))
    for energy, b_mid, i_radial_norm in channels["UHE"]["radial_profiles"]:
        plt.plot(
            b_mid,
            i_radial_norm,
            linestyle="-",
            label=rf"UHE $E_\nu={energy:.0e}\,\mathrm{{GeV}}$",
        )

    for energy, b_mid, i_radial_norm in channels["MeV"]["radial_profiles"]:
        plt.plot(
            b_mid,
            i_radial_norm,
            linestyle="--",
            label=rf"MeV $E_\nu={energy:.0e}\,\mathrm{{GeV}}$",
        )

    plt.xlabel(r"$b=\sqrt{\alpha^2+\beta^2}$")
    plt.ylabel(r"$\langle I(b)\rangle/I_{\rm max}$")
    plt.title(suffix)
    plt.legend(fontsize=8)
    plt.tight_layout()
    plt.savefig(plot_dir / "diagnostic_radial_profiles.png", dpi=250)

    plt.figure(figsize=(7, 5))
    plt.semilogx(
        energies,
        channels["UHE"]["contrast_shadow"],
        marker="o",
        label="UHE",
    )
    plt.semilogx(
        energies,
        channels["MeV"]["contrast_shadow"],
        marker="s",
        label="MeV",
    )
    plt.xlabel(r"$E_\nu\ [\mathrm{GeV}]$")
    plt.ylabel(r"$I_{\rm ring}/I_{\rm center}$")
    plt.title(suffix)
    plt.legend()
    plt.tight_layout()
    plt.savefig(plot_dir / "diagnostic_shadow_contrast_vs_energy.png", dpi=250)

    plt.close("all")

    print(f"Density: {rho_tag}")
    print("Used files:")
    for fname in files:
        print(f"  {fname}")

    print(f"Saved: {plot_dir / 'diagnostic_total_flux_vs_energy.png'}")
    print(f"Saved: {plot_dir / 'diagnostic_tau_vs_energy.png'}")
    print(f"Saved: {plot_dir / 'diagnostic_survival_vs_energy.png'}")
    print(f"Saved: {plot_dir / 'diagnostic_radial_profiles.png'}")
    print(f"Saved: {plot_dir / 'diagnostic_shadow_contrast_vs_energy.png'}")

    return {
        "rho_tag": rho_tag,
        "torus_rho0": torus_rho0,
        "energies": energies,
        "UHE_tau_I_weighted": channels["UHE"]["tau_I_weighted"],
        "MeV_tau_I_weighted": channels["MeV"]["tau_I_weighted"],
    }


def save_density_comparison(diagnostics_by_density):
    if len(diagnostics_by_density) < 2:
        return

    plot_dir = Path("plots") / "energy_diagnostics_all_densities"
    plot_dir.mkdir(parents=True, exist_ok=True)

    sorted_diagnostics = sorted(
        diagnostics_by_density,
        key=lambda item: item["torus_rho0"],
    )

    plt.figure(figsize=(7, 5))
    for diagnostics in sorted_diagnostics:
        label = (
            rf"$\rho_{{0,\rm torus}}="
            rf"{sci_latex_compact(diagnostics['torus_rho0'])}"
            rf"\,\mathrm{{g\,cm^{{-3}}}}$"
        )

        plt.loglog(
            diagnostics["energies"],
            diagnostics["UHE_tau_I_weighted"],
            marker="o",
            label=label,
        )

    plt.xlabel(r"$E_\nu\ [\mathrm{GeV}]$")
    plt.ylabel(r"$\langle\tau\rangle_I$")
    plt.title("Neutrino UHE type")
    plt.legend(fontsize=8)
    plt.tight_layout()

    uhe_output = plot_dir / "diagnostic_tau_vs_energy_all_torus_densities_UHE.png"
    plt.savefig(uhe_output, dpi=250)
    plt.close()

    plt.figure(figsize=(7, 5))
    for diagnostics in sorted_diagnostics:
        label = (
            rf"$\rho_{{0,\rm torus}}="
            rf"{sci_latex_compact(diagnostics['torus_rho0'])}"
            rf"\,\mathrm{{g\,cm^{{-3}}}}$"
        )

        plt.loglog(
            diagnostics["energies"],
            diagnostics["MeV_tau_I_weighted"],
            marker="o",
            label=label,
        )

    plt.xlabel(r"$E_\nu\ [\mathrm{GeV}]$")
    plt.ylabel(r"$\langle\tau\rangle_I$")
    plt.title("Neutrino MeV type")
    plt.legend(fontsize=8)
    plt.tight_layout()

    mev_output = plot_dir / "diagnostic_tau_vs_energy_all_torus_densities_MeV.png"
    plt.savefig(mev_output, dpi=250)
    plt.close()

    print(f"Saved: {uhe_output}")
    print(f"Saved: {mev_output}")


def main():
    if len(sys.argv) > 1:
        torus_rho0 = float(sys.argv[1])
        rho_tag = rho_folder_tag(torus_rho0)
        files = find_rho_files(rho_tag)

        if not files:
            raise FileNotFoundError(
                f"No files found for {rho_tag}: "
                "output/images/kerr_image*_rho0_torus_*_Enu_*.bin or .dat"
            )

        save_diagnostics(files, torus_rho0, rho_tag)
        return

    rho_tags = discover_rho_tags()

    if not rho_tags:
        torus_rho0 = float(os.environ.get("TORUS_RHO0", "1.0e-2"))
        rho_tag = rho_folder_tag(torus_rho0)
        legacy_files = find_legacy_files()

        if not legacy_files:
            raise FileNotFoundError(
                "No files found: output/images/kerr_image*_Enu_*.bin or .dat"
            )

        print(
            "No rho0_torus tags found in output/images; "
            f"using legacy files as {rho_tag}."
        )
        save_diagnostics(legacy_files, torus_rho0, rho_tag)
        return

    diagnostics_by_density = []

    for rho_tag in rho_tags:
        torus_rho0 = rho_value_from_tag(rho_tag)
        diagnostics = save_diagnostics(find_rho_files(rho_tag), torus_rho0, rho_tag)
        if diagnostics is not None:
            diagnostics_by_density.append(diagnostics)

    save_density_comparison(diagnostics_by_density)


if __name__ == "__main__":
    main()
