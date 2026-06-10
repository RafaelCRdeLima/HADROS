"""Generate publication-quality tau vs energy figure for all torus densities — IIM model.

Mirror of plot_tau_publication.py but restricted to IIM-tagged output files:
    kerr_image_cuda_cache_IIM_rho0_torus_*

Run from the BH_Torus_RTX directory:
    python scripts/plot_tau_publication_iim.py

Outputs:
    paper/diagnostic_tau_vs_energy_all_torus_densities_IIM_pub.png
"""

import glob
import re
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

PAPER_DIR = Path(__file__).resolve().parents[2] / "paper"

EXCLUDE_RHO = 1e-2   # transparent benchmark density — excluded from plot

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

_PALETTE = [
    "#0072B2",
    "#009E73",
    "#D55E00",
    "#CC79A7",
    "#56B4E9",
    "#E69F00",
]
_LINESTYLES = ["-", "--", "-.", ":", (0,(5,1)), (0,(3,1,1,1))]

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def load_kerr_image(path):
    return np.loadtxt(path)


def rho_value_from_tag(tag):
    m = re.fullmatch(r"rho0_torus_([0-9]+(?:p[0-9]+)?)e(-?[0-9]+)", tag)
    if not m:
        raise ValueError(f"Invalid tag: {tag}")
    mantissa = float(m.group(1).replace("p", "."))
    return mantissa * 10.0 ** int(m.group(2))


def discover_iim_rho_tags():
    tags = set()
    pattern = re.compile(r"IIM_(rho0_torus_[0-9]+(?:p[0-9]+)?e-?[0-9]+)_Enu_")
    for fname in glob.glob("output/images/kerr_image_cuda_cache_IIM_rho0_torus_*_Enu_*.dat"):
        m = pattern.search(Path(fname).name)
        if m:
            tags.add(m.group(1))
    return sorted(tags, key=rho_value_from_tag)


def find_iim_rho_files(rho_tag):
    return glob.glob(
        f"output/images/kerr_image_cuda_cache_IIM_{rho_tag}_Enu_*.dat"
    )


def extract_energy(filename):
    m = re.search(r"Enu_([0-9]+e_[0-9]+)", filename)
    if not m:
        raise ValueError(f"Cannot extract energy from {filename}")
    return float(m.group(1).replace("_", "+"))


def select_files_by_energy(files):
    best = {}
    for f in files:
        e = extract_energy(f)
        if e not in best:
            best[e] = f
    return [best[e] for e in sorted(best)]


def compute_tau_weighted(data):
    tau   = data[:, 4]
    i_obs = data[:, 6]
    flux  = np.sum(i_obs)
    if flux > 0:
        return np.sum(i_obs * tau) / flux
    return np.nan


def load_diagnostics(rho_tag):
    files = select_files_by_energy(find_iim_rho_files(rho_tag))
    if not files:
        return None, None
    energies, taus = [], []
    for fname in files:
        data = load_kerr_image(fname)
        if data.shape[1] <= 6:
            continue
        energies.append(extract_energy(fname))
        taus.append(compute_tau_weighted(data))
    if not energies:
        return None, None
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


def make_figure(all_diagnostics):
    fig, ax = plt.subplots(figsize=(3.5, 2.9))

    for idx, (rho0, energies, taus) in enumerate(all_diagnostics):
        color = _PALETTE[idx % len(_PALETTE)]
        ls    = _LINESTYLES[idx % len(_LINESTYLES)]
        label = (
            r"$\rho_{0,\rm torus}="
            + sci_latex_compact(rho0).strip("$")
            + r"\,\mathrm{g\,cm^{-3}}$"
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
    rho_tags = discover_iim_rho_tags()
    if not rho_tags:
        print("No IIM rho0_torus data found in output/images/. Run from BH_Torus_RTX/.")
        print("Make sure run_iim_tau_scan.sh has been executed first.")
        sys.exit(1)

    all_diagnostics = []
    for tag in rho_tags:
        rho0 = rho_value_from_tag(tag)
        if abs(rho0 - EXCLUDE_RHO) / EXCLUDE_RHO < 0.01:
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

    fig = make_figure(all_diagnostics)
    out = PAPER_DIR / "diagnostic_tau_vs_energy_all_torus_densities_IIM_pub.png"
    fig.savefig(out)
    plt.close(fig)
    print(f"Saved: {out}")


if __name__ == "__main__":
    main()
