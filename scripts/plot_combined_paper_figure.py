"""Generate the combined RGB-overlay UHE+MeV figure for the paper.

Layout: 3 rows x 2 columns
  - Rows: scenarios (rho0, E_uhe)
  - Columns: GBW (left) | IIM (right)
  - Each panel: single RGB composite image
      Yellow channel (R+G) = UHE neutrinos (sqrt-normalised)
      Blue  channel (B)    = MeV neutrinos (sqrt-normalised)

Scenarios:
  1. rho0=1e-2  g/cm^3,  E_uhe=1e12 GeV  (transparent torus)
  2. rho0=1e6   g/cm^3,  E_uhe=1e6  GeV  (intermediate attenuation)
  3. rho0=1e11  g/cm^3,  E_uhe=1e12 GeV  (opaque torus)

Run from BH_Torus_RTX/:
    python scripts/plot_combined_paper_figure.py
"""

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.cm import ScalarMappable

PAPER_DIR = Path(__file__).resolve().parents[2] / "paper"

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
# Scenarios
# ---------------------------------------------------------------------------
SCENARIOS = [
    {
        "label": "transparent",
        "rho0": 1e-2,
        "Enu_GeV": 1e12,
        "mev_enu": 10.0,
        "mev_norm": 1.0,
        "cam_theta": 80.0,
    },
    {
        "label": "intermediate",
        "rho0": 1e6,
        "Enu_GeV": 1e6,
        "mev_enu": 10.0,
        "mev_norm": 1.0,
        "cam_theta": 80.0,
    },
    {
        "label": "opaque",
        "rho0": 1e11,
        "Enu_GeV": 1e12,
        "mev_enu": 10.0,
        "mev_norm": 1.0,
        "cam_theta": 80.0,
    },
]

OUTPUT_NAME = "comparison_combined_3scenarios.png"

# Colour-map stubs for the colorbars only (actual image is RGB)
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
    s = fmt % value
    return s.replace("+", "_").replace("-", "m").replace(".", "p")


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
    rt = make_run_tag(Enu_GeV, mev_enu, mev_norm, cam_theta)
    candidates = [
        Path(f"output/images/kerr_image_cuda_cache_{model}_{rho_tag(rho0)}_{rt}.dat"),
        # Legacy GBW naming without model tag
        Path(f"output/images/kerr_image_cuda_cache_{rho_tag(rho0)}_{rt}.dat"),
        Path(f"output/images/kerr_image_cuda_cache_{rho_tag(rho0)}_{rt}.bin"),
    ]
    for p in candidates:
        if p.exists():
            return p
    raise FileNotFoundError(
        f"No image file for model={model}, rho={rho0:.0e}, E={Enu_GeV:.0e}.\n"
        "Tried:\n" + "\n".join(f"  {p}" for p in candidates)
    )


def load_dat(path):
    return np.loadtxt(path, comments="#")


def make_rgb_overlay(data):
    """Return (rgb_array, extent) where rgb encodes UHE=yellow, MeV=blue."""
    pi = data[:, 0].astype(int)
    pj = data[:, 1].astype(int)
    alpha = data[:, 2]
    beta  = data[:, 3]
    nx = pi.max() + 1
    ny = pj.max() + 1

    I_uhe = data[:, 6]
    I_mev = data[:, 10] if data.shape[1] > 10 else np.zeros_like(I_uhe)

    row = (ny - 1) - pj

    def norm_sqrt(arr):
        mx = arr.max()
        if mx > 0:
            arr = arr / mx
        return np.sqrt(np.clip(arr, 0.0, 1.0))

    g_uhe = np.zeros((ny, nx))
    g_mev = np.zeros((ny, nx))
    g_uhe[row, pi] = norm_sqrt(I_uhe)
    g_mev[row, pi] = norm_sqrt(I_mev)

    # RGB composite: yellow = UHE (R+G), blue = MeV (B)
    rgb = np.zeros((ny, nx, 3))
    rgb[:, :, 0] = g_uhe   # R
    rgb[:, :, 1] = g_uhe   # G  → R+G = yellow
    rgb[:, :, 2] = g_mev   # B  → blue
    rgb = np.clip(rgb, 0.0, 1.0)

    extent = [alpha.min(), alpha.max(), -beta.max(), -beta.min()]
    return rgb, extent


# ---------------------------------------------------------------------------
# Main figure: 3 rows x 2 columns of RGB overlays
# ---------------------------------------------------------------------------

def make_combined_figure():
    fig = plt.figure(figsize=(6.4, 8.0))
    gs = gridspec.GridSpec(
        3, 2,
        figure=fig,
        hspace=0.28,
        wspace=0.06,
        left=0.09, right=0.88,
        top=0.95, bottom=0.06,
    )

    col_titles = ["GBW", "IIM"]

    for irow, sc in enumerate(SCENARIOS):
        kw = {k: sc[k] for k in ("rho0", "Enu_GeV", "mev_enu", "mev_norm", "cam_theta")}

        rho_exp = int(round(np.log10(sc["rho0"])))
        E_exp   = int(round(np.log10(sc["Enu_GeV"])))

        rho_str = (r"$\rho_{0,\mathrm{t}}=10^{" + str(rho_exp)
                   + r"}\,\mathrm{g\,cm^{-3}}$")
        E_str   = r"$E_\nu^{\rm UHE}=10^{" + str(E_exp) + r"}\,\mathrm{GeV}$"

        for icol, model in enumerate(["GBW", "IIM"]):
            path = find_dat(model, **kw)
            print(f"[{sc['label']}] {model}: {path}")
            data = load_dat(path)
            rgb, extent = make_rgb_overlay(data)

            ax = fig.add_subplot(gs[irow, icol])
            ax.imshow(
                rgb,
                origin="lower",
                extent=extent,
                aspect="equal",
                interpolation="nearest",
            )

            # Column title on top row only
            if irow == 0:
                ax.set_title(col_titles[icol], fontsize=10, pad=4)

            # x-axis label on bottom row only
            if irow == 2:
                ax.set_xlabel(r"$\alpha\;[r_g]$", labelpad=1)
            else:
                ax.tick_params(labelbottom=False)

            # y-axis label on left column only
            if icol == 0:
                ax.set_ylabel(r"$\beta\;[r_g]$", labelpad=1)
            else:
                ax.tick_params(labelleft=False)

            # Row label: top-left corner of the GBW (left) panel only
            if icol == 0:
                ax.text(
                    0.03, 0.97,
                    rho_str + "\n" + E_str,
                    transform=ax.transAxes,
                    va="top", ha="left",
                    fontsize=7,
                    color="white",
                    bbox=dict(boxstyle="round,pad=0.2", fc="black", alpha=0.45,
                              ec="none"),
                )

    # Shared colorbars on the far right
    norm = Normalize(vmin=0.0, vmax=1.0)
    cax_uhe = fig.add_axes([0.91, 0.54, 0.022, 0.39])
    cax_mev = fig.add_axes([0.91, 0.08, 0.022, 0.39])
    fig.colorbar(ScalarMappable(norm=norm, cmap=UHE_CMAP), cax=cax_uhe).set_label(
        r"UHE $\sqrt{I/I_{\max}}$", fontsize=8
    )
    fig.colorbar(ScalarMappable(norm=norm, cmap=MEV_CMAP), cax=cax_mev).set_label(
        r"MeV $\sqrt{I/I_{\max}}$", fontsize=8
    )

    return fig


def main():
    fig = make_combined_figure()
    out = PAPER_DIR / OUTPUT_NAME
    fig.savefig(out)
    plt.close(fig)
    print(f"Saved: {out}")


if __name__ == "__main__":
    main()
