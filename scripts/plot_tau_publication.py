"""Generate publication-quality tau vs energy figure for all torus densities.

Run from the BH_Torus_RTX directory:
    python scripts/plot_tau_publication.py

Outputs:
    paper/diagnostic_tau_vs_energy_all_torus_densities_UHE_pub.png
"""

import glob
import os
import re
import struct
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# ---------------------------------------------------------------------------
# Physical / file constants (copied from plot_energy_diagnostics.py)
# ---------------------------------------------------------------------------
KIMG_MAGIC = 0x4B494D47
KIMG_HEADER = struct.Struct("<7i4d")

C_CGS   = 2.99792458e10
G_CGS   = 6.67430e-8
MSUN_G  = 1.98847e33
M_U_G   = 1.66053906660e-24

PAPER_DIR = Path(__file__).resolve().parents[2] / "paper"

# Density to exclude
EXCLUDE_RHO = 1e-2

# ---------------------------------------------------------------------------
# Matplotlib publication style
# ---------------------------------------------------------------------------
plt.rcParams.update({
    "text.usetex": False,
    "font.family": "serif",
    "font.serif": ["DejaVu Serif", "Times New Roman", "Times", "serif"],
    "font.size": 9,
    "axes.labelsize": 10,
    "axes.titlesize": 10,
    "legend.fontsize": 7.5,
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
    "xtick.direction": "in",
    "ytick.direction": "in",
    "xtick.top": True,
    "ytick.right": True,
    "axes.linewidth": 0.8,
    "lines.linewidth": 1.5,
    "figure.dpi": 300,
    "savefig.dpi": 300,
    "savefig.bbox": "tight",
    "savefig.pad_inches": 0.02,
})

# Color palette — 6 levels, visually distinct, colorblind-tolerant
_PALETTE = [
    "#0072B2",   # blue
    "#009E73",   # green
    "#D55E00",   # orange-red
    "#CC79A7",   # pink
    "#56B4E9",   # sky-blue
    "#E69F00",   # amber
]

_LINESTYLES = ["-", "--", "-.", ":", (0,(5,1)), (0,(3,1,1,1))]

# ---------------------------------------------------------------------------
# Binary loader (identical to plot_energy_diagnostics.py)
# ---------------------------------------------------------------------------

def load_kerr_image_binary(path):
    with open(path, "rb") as f:
        header_bytes = f.read(KIMG_HEADER.size)
        if len(header_bytes) != KIMG_HEADER.size:
            raise ValueError(f"Invalid header: {path}")
        (magic, version, backend, nx, ny, ncols, reserved,
         enu, mev_enu, mev_norm, cam_theta) = KIMG_HEADER.unpack(header_bytes)
        if magic != KIMG_MAGIC:
            raise ValueError(f"Invalid magic: {path}")
        if version != 1:
            raise ValueError(f"Unsupported version {version}: {path}")
        payload = np.fromfile(f, dtype="<f8")
    expected = nx * ny * ncols
    if payload.size != expected:
        raise ValueError(
            f"Payload size mismatch in {path}: expected {expected}, got {payload.size}"
        )
    return payload.reshape((nx * ny, ncols))


def load_kerr_image(path):
    path = Path(path)
    if path.suffix == ".bin":
        return load_kerr_image_binary(path)
    return np.loadtxt(path)


# ---------------------------------------------------------------------------
# File discovery helpers
# ---------------------------------------------------------------------------

def rho_value_from_tag(tag):
    match = re.fullmatch(r"rho0_torus_([0-9]+(?:p[0-9]+)?)e(-?[0-9]+)", tag)
    if not match:
        raise ValueError(f"Invalid tag: {tag}")
    mantissa = float(match.group(1).replace("p", "."))
    exponent = int(match.group(2))
    return mantissa * 10.0**exponent


def discover_rho_tags():
    tags = set()
    pattern = re.compile(r"(rho0_torus_[0-9]+(?:p[0-9]+)?e-?[0-9]+)_Enu_")
    for fname in glob.glob("output/images/kerr_image*_rho0_torus_*_Enu_*"):
        m = pattern.search(Path(fname).name)
        if m:
            tags.add(m.group(1))
    return sorted(tags, key=rho_value_from_tag)


def find_rho_files(rho_tag):
    return (
        glob.glob(f"output/images/kerr_image_cuda_cache_{rho_tag}_Enu_*.bin")
        + glob.glob(f"output/images/kerr_image_cuda_{rho_tag}_Enu_*.bin")
        + glob.glob(f"output/images/kerr_image_{rho_tag}_Enu_*.bin")
        + glob.glob(f"output/images/kerr_image_cuda_cache_{rho_tag}_Enu_*.dat")
        + glob.glob(f"output/images/kerr_image_cuda_{rho_tag}_Enu_*.dat")
        + glob.glob(f"output/images/kerr_image_{rho_tag}_Enu_*.dat")
    )


def extract_energy(filename):
    m = re.search(r"Enu_([0-9]+e_[0-9]+)", filename)
    if not m:
        raise ValueError(f"Cannot extract energy from {filename}")
    return float(m.group(1).replace("_", "+"))


def file_priority(filename):
    path = Path(filename)
    is_cc = path.name.startswith("kerr_image_cuda_cache_")
    is_cu = path.name.startswith("kerr_image_cuda_")
    is_bin = path.suffix == ".bin"
    if is_cc and is_bin:  return 0
    if is_cu and is_bin:  return 1
    if is_bin:            return 2
    if is_cc:             return 3
    if is_cu:             return 4
    return 5


def select_files_by_energy(files):
    best = {}
    for f in files:
        e = extract_energy(f)
        if e not in best or file_priority(f) < file_priority(best[e]):
            best[e] = f
    return [best[e] for e in sorted(best)]


# ---------------------------------------------------------------------------
# Diagnostic computation
# ---------------------------------------------------------------------------
UHE_TAU_COL = 4
UHE_I_COL   = 6


def compute_tau_weighted(data):
    tau   = data[:, UHE_TAU_COL]
    i_obs = data[:, UHE_I_COL]
    flux  = np.sum(i_obs)
    if flux > 0:
        return np.sum(i_obs * tau) / flux
    return np.nan


def load_diagnostics(rho_tag):
    files = select_files_by_energy(find_rho_files(rho_tag))
    if not files:
        return None, None
    energies, taus = [], []
    for fname in files:
        data = load_kerr_image(fname)
        if data.shape[1] <= UHE_I_COL:
            continue
        energies.append(extract_energy(fname))
        taus.append(compute_tau_weighted(data))
    order = np.argsort(energies)
    return np.array(energies)[order], np.array(taus)[order]


# ---------------------------------------------------------------------------
# Figure
# ---------------------------------------------------------------------------

def sci_latex_compact(value):
    mantissa_text, exponent = f"{value:.0e}".split("e")
    mantissa = float(mantissa_text)
    exp_int  = int(exponent)
    if np.isclose(mantissa, 1.0):
        return rf"$10^{{{exp_int}}}$"
    return rf"${mantissa_text}\times10^{{{exp_int}}}$"


def make_publication_figure(all_diagnostics):
    # Figure width ~ 84 mm (single-column MNRAS) → 3.31 in
    fig, ax = plt.subplots(figsize=(3.5, 2.9))

    for idx, (rho0, energies, taus) in enumerate(all_diagnostics):
        color = _PALETTE[idx % len(_PALETTE)]
        ls    = _LINESTYLES[idx % len(_LINESTYLES)]
        label = (
            r"$\rho_{0,\rm torus}=" + sci_latex_compact(rho0).strip("$") + r"\,\mathrm{g\,cm^{-3}}$"
        )
        ax.loglog(energies, taus, color=color, linestyle=ls, linewidth=1.5, label=label)

    ax.set_xlabel(r"$E_\nu\ [\mathrm{GeV}]$")
    ax.set_ylabel(r"$\langle\tau\rangle_I$")

    ax.legend(
        loc="lower right",
        framealpha=0.90,
        edgecolor="0.7",
        handlelength=1.5,
        fontsize=6.5,
        borderpad=0.5,
        labelspacing=0.3,
        handletextpad=0.4,
    )

    ax.xaxis.set_tick_params(which="both", direction="in")
    ax.yaxis.set_tick_params(which="both", direction="in")

    fig.tight_layout()
    return fig


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    rho_tags = discover_rho_tags()
    if not rho_tags:
        print("No rho0_torus data found in output/images/. Run from BH_Torus_RTX/.")
        sys.exit(1)

    all_diagnostics = []
    for tag in rho_tags:
        rho0 = rho_value_from_tag(tag)
        if np.isclose(rho0, EXCLUDE_RHO, rtol=0.01):
            print(f"Skipping {tag} (excluded density 1e-2)")
            continue
        energies, taus = load_diagnostics(tag)
        if energies is None or len(energies) == 0:
            print(f"Skipping {tag}: no usable data")
            continue
        all_diagnostics.append((rho0, energies, taus))
        print(f"Loaded {tag}: {len(energies)} energy points")

    if not all_diagnostics:
        print("No data loaded.")
        sys.exit(1)

    fig = make_publication_figure(all_diagnostics)

    out_path = PAPER_DIR / "diagnostic_tau_vs_energy_all_torus_densities_UHE_pub.png"
    fig.savefig(out_path)
    plt.close(fig)
    print(f"Saved: {out_path}")


if __name__ == "__main__":
    main()
