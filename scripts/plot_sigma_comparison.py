"""Publication-quality GBW vs IIM sigma_nuN comparison figure.

Run from BH_Torus_RTX/:
    python scripts/plot_sigma_comparison.py
Saves to: ../../paper/sigma_nuN_CC_GBW_vs_IIM.png
"""

from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

PAPER_DIR = Path(__file__).resolve().parents[2] / "paper"

GBW_FILE = Path("data/sigma/sigma_nuN_CC_GBW.dat")
IIM_FILE = Path("data/sigma/sigma_nuN_CC_IIM.dat")

plt.rcParams.update({
    "text.usetex": False,
    "font.family": "serif",
    "font.serif": ["DejaVu Serif", "Times New Roman", "serif"],
    "font.size": 9,
    "axes.labelsize": 10,
    "axes.titlesize": 10,
    "legend.fontsize": 8.5,
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
    "xtick.direction": "in",
    "ytick.direction": "in",
    "xtick.top": True,
    "ytick.right": True,
    "axes.linewidth": 0.8,
    "lines.linewidth": 1.6,
    "figure.dpi": 300,
    "savefig.dpi": 300,
    "savefig.bbox": "tight",
    "savefig.pad_inches": 0.03,
})


def load(path):
    data = np.loadtxt(path, comments="#")
    return data[:, 0], data[:, 2]   # E_GeV, sigma_cm2


def main():
    E_gbw, s_gbw = load(GBW_FILE)
    E_iim, s_iim = load(IIM_FILE)

    fig, ax = plt.subplots(figsize=(3.5, 3.0))

    ax.loglog(E_gbw, s_gbw, color="#0072B2", linestyle="-",
              label="GBW")
    ax.loglog(E_iim, s_iim, color="#D55E00", linestyle="--",
              label="IIM")

    ax.set_xlabel(r"$E_\nu\ [\mathrm{GeV}]$")
    ax.set_ylabel(r"$\sigma_{\nu N}^{CC}\ [\mathrm{cm}^2]$")

    ax.legend(loc="upper left", framealpha=0.95, edgecolor="0.7",
              handlelength=1.8)

    ax.xaxis.set_tick_params(which="both", direction="in")
    ax.yaxis.set_tick_params(which="both", direction="in")

    fig.tight_layout()

    out = PAPER_DIR / "sigma_nuN_CC_GBW_vs_IIM.png"
    fig.savefig(out)
    plt.close(fig)
    print(f"Saved: {out}")


if __name__ == "__main__":
    main()
