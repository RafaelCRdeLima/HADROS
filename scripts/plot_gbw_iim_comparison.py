"""Generate GBW vs IIM comparison figures for the paper.

Creates 4-panel figures (GBW-UHE | IIM-UHE | GBW-MeV | IIM-MeV) for each
scenario defined in SCENARIOS.  The MeV channel does not depend on the DIS
cross-section model, so the two MeV panels are included for completeness.

Run from BH_Torus_RTX/:
    python scripts/plot_gbw_iim_comparison.py
"""

import struct
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.cm import ScalarMappable

PAPER_DIR = Path(__file__).resolve().parents[2] / "paper"

KIMG_MAGIC = 0x4B494D47
KIMG_HEADER = struct.Struct("<7i4d")

# ---------------------------------------------------------------------------
# Publication rcParams
# ---------------------------------------------------------------------------
plt.rcParams.update({
    "text.usetex": False,
    "font.family": "serif",
    "font.serif": ["DejaVu Serif", "Times New Roman", "serif"],
    "font.size": 9,
    "axes.labelsize": 9,
    "axes.titlesize": 9,
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
    "figure.dpi": 300,
    "savefig.dpi": 300,
    "savefig.bbox": "tight",
    "savefig.pad_inches": 0.03,
})

# ---------------------------------------------------------------------------
# Scenarios: (label, rho0, E_uhe_GeV, mev_enu_MeV, mev_norm, cam_theta)
# ---------------------------------------------------------------------------
SCENARIOS = [
    {
        "label": "transparent",
        "rho0": 1e-2,
        "Enu_GeV": 1e12,
        "mev_enu": 10.0,
        "mev_norm": 1.0,
        "cam_theta": 80.0,
        "paper_out": "comparison_gbw_iim_rho1em2_E1e12.png",
    },
    {
        "label": "attenuating",
        "rho0": 1e6,
        "Enu_GeV": 1e3,
        "mev_enu": 10.0,
        "mev_norm": 1.0,
        "cam_theta": 80.0,
        "paper_out": "comparison_gbw_iim_rho1e6_E1e3.png",
    },
]

# ---------------------------------------------------------------------------
# Colour maps
# ---------------------------------------------------------------------------
UHE_CMAP = LinearSegmentedColormap.from_list(
    "uhe_yellow", [(0.0, 0.0, 0.0), (1.0, 1.0, 0.0)]
)
MEV_CMAP = LinearSegmentedColormap.from_list(
    "mev_blue", [(0.0, 0.0, 0.0), (0.0, 0.35, 1.0)]
)

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def sci_tag(value, fmt="%.0e"):
    import re
    s = fmt % value
    s = s.replace("+", "_").replace("-", "m").replace(".", "p")
    return s


def rho_tag(value):
    m, e = ("%.0e" % value).split("e")
    return f"rho0_torus_{m}e{int(e)}"


def make_run_tag(Enu_GeV, mev_enu, mev_norm, cam_theta):
    return (
        f"Enu_{sci_tag(Enu_GeV)}"
        f"_MeVEnu_{sci_tag(mev_enu)}"
        f"_MeVNorm_{sci_tag(mev_norm)}"
        f"_CamTheta_{cam_theta:.1f}".replace(".", "p")
    )


def find_dat(model, rho0, Enu_GeV, mev_enu, mev_norm, cam_theta):
    """Resolve path to the output .dat file for a given model and config."""
    rt = make_run_tag(Enu_GeV, mev_enu, mev_norm, cam_theta)
    rt_cam = rt.replace("_CamTheta_", "_CamTheta_")

    candidates = [
        # IIM/GBW tagged (new naming from modified compute_kerr_image_from_cache)
        Path(f"output/images/kerr_image_cuda_cache_{model}_{rho_tag(rho0)}_{rt}.dat"),
        # GBW legacy naming (original CUDA output, no model tag)
        Path(f"output/images/kerr_image_cuda_cache_{rho_tag(rho0)}_{rt}.dat"),
        Path(f"output/images/kerr_image_cuda_cache_{rho_tag(rho0)}_{rt}.bin"),
    ]
    for p in candidates:
        if p.exists():
            return p
    raise FileNotFoundError(
        f"No image file for model={model}, rho={rho0:.0e}, E={Enu_GeV:.0e}.\n"
        f"Tried:\n" + "\n".join(f"  {p}" for p in candidates)
    )


def load_binary(path):
    with open(path, "rb") as f:
        hdr = f.read(KIMG_HEADER.size)
        if len(hdr) != KIMG_HEADER.size:
            raise ValueError(f"Short header: {path}")
        magic, version, backend, nx, ny, ncols, _, enu, mev_enu, mev_norm, cam = \
            KIMG_HEADER.unpack(hdr)
        if magic != KIMG_MAGIC:
            raise ValueError(f"Bad magic: {path}")
        payload = np.fromfile(f, dtype="<f8")
    expected = nx * ny * ncols
    if payload.size != expected:
        raise ValueError(f"Payload mismatch in {path}")
    return payload.reshape((nx * ny, ncols))


def load_dat(path):
    return np.loadtxt(path, comments="#")


def load_image(path):
    path = Path(path)
    if path.suffix == ".bin":
        return load_binary(path)
    return load_dat(path)


def make_grids(data):
    pi = data[:, 0].astype(int)
    pj = data[:, 1].astype(int)
    alpha = data[:, 2]
    beta = data[:, 3]
    nx = pi.max() + 1
    ny = pj.max() + 1

    I_uhe = data[:, 6]
    I_mev = data[:, 10] if data.shape[1] > 10 else np.zeros_like(I_uhe)

    row = (ny - 1) - pj

    def norm_sqrt(arr):
        mx = arr.max()
        if mx > 0:
            arr = arr / mx
        return np.sqrt(np.clip(arr, 0, 1))

    g_uhe = np.zeros((ny, nx))
    g_mev = np.zeros((ny, nx))
    g_uhe[row, pi] = norm_sqrt(I_uhe)
    g_mev[row, pi] = norm_sqrt(I_mev)

    extent = [alpha.min(), alpha.max(), -beta.max(), -beta.min()]
    return g_uhe, g_mev, extent


# ---------------------------------------------------------------------------
# Figure builder
# ---------------------------------------------------------------------------

def make_comparison_figure(scenario):
    gbw_path = find_dat("GBW", **{k: scenario[k] for k in
                                   ("rho0", "Enu_GeV", "mev_enu", "mev_norm", "cam_theta")})
    iim_path = find_dat("IIM", **{k: scenario[k] for k in
                                   ("rho0", "Enu_GeV", "mev_enu", "mev_norm", "cam_theta")})

    print(f"GBW: {gbw_path}")
    print(f"IIM: {iim_path}")

    gbw = load_image(gbw_path)
    iim = load_image(iim_path)

    g_uhe_gbw, g_mev_gbw, ext = make_grids(gbw)
    g_uhe_iim, g_mev_iim, _   = make_grids(iim)

    rho_str  = f"$\\rho_{{0,\\mathrm{{t}}}}=10^{{{int(round(np.log10(scenario['rho0'])))}}}\\,\\mathrm{{g\\,cm^{{-3}}}}$"
    E_exp    = int(round(np.log10(scenario["Enu_GeV"])))
    E_str    = f"$E_\\nu=10^{{{E_exp}}}\\,\\mathrm{{GeV}}$"
    mev_str  = f"$E_\\nu^{{\\rm MeV}}={scenario['mev_enu']:.0f}\\,\\mathrm{{MeV}}$"

    # 2-row × 2-column layout
    # Row 0: GBW UHE  | IIM UHE
    # Row 1: GBW MeV  | IIM MeV
    fig = plt.figure(figsize=(7.0, 6.4))
    gs  = gridspec.GridSpec(
        2, 2,
        figure=fig,
        hspace=0.35,
        wspace=0.05,
        left=0.09, right=0.88,
        top=0.92, bottom=0.09,
    )

    panels = [
        (0, 0, g_uhe_gbw, UHE_CMAP, f"GBW — UHE ({E_str})",        r"$\sqrt{I/I_{\max}}$"),
        (0, 1, g_uhe_iim, UHE_CMAP, f"IIM — UHE ({E_str})",        r"$\sqrt{I/I_{\max}}$"),
        (1, 0, g_mev_gbw, MEV_CMAP, f"GBW — MeV ({mev_str})",      r"$\sqrt{I/I_{\max}}$"),
        (1, 1, g_mev_iim, MEV_CMAP, f"IIM — MeV ({mev_str})",      r"$\sqrt{I/I_{\max}}$"),
    ]

    axes = []
    for row, col, grid, cmap, title, _ in panels:
        ax = fig.add_subplot(gs[row, col])
        ax.imshow(
            grid,
            origin="lower",
            extent=ext,
            aspect="equal",
            cmap=cmap,
            vmin=0.0,
            vmax=1.0,
        )
        ax.set_title(title, fontsize=8, pad=3)
        ax.set_xlabel(r"$\alpha$", labelpad=1)
        if col == 0:
            ax.set_ylabel(r"$\beta$", labelpad=1)
        else:
            ax.tick_params(labelleft=False)
        axes.append(ax)

    # Shared colorbars on the right
    norm = Normalize(vmin=0.0, vmax=1.0)
    cax_uhe = fig.add_axes([0.90, 0.54, 0.022, 0.36])
    cax_mev = fig.add_axes([0.90, 0.09, 0.022, 0.36])
    fig.colorbar(ScalarMappable(norm=norm, cmap=UHE_CMAP), cax=cax_uhe).set_label(
        r"UHE $\sqrt{I/I_{\max}}$", fontsize=7
    )
    fig.colorbar(ScalarMappable(norm=norm, cmap=MEV_CMAP), cax=cax_mev).set_label(
        r"MeV $\sqrt{I/I_{\max}}$", fontsize=7
    )

    fig.suptitle(
        rho_str,
        fontsize=9, y=0.97
    )

    return fig


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    for scenario in SCENARIOS:
        try:
            fig = make_comparison_figure(scenario)
        except FileNotFoundError as e:
            print(f"SKIP {scenario['label']}: {e}")
            continue

        out = PAPER_DIR / scenario["paper_out"]
        fig.savefig(out)
        plt.close(fig)
        print(f"Saved: {out}")


if __name__ == "__main__":
    main()
