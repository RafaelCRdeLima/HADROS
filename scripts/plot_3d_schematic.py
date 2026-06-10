"""3D schematic: Kerr BH + torus + UHE ring + Schwarzschild null geodesics.

ALL geodesics originate from the SAME camera position (r=52, θ=80°, φ=0).
Full 3D Schwarzschild integration using conserved (E, L, K):
  L = α sinθ_cam   (azimuthal angular momentum)
  K = α² + β²       (Carter constant; α,β = image-plane coords in r_g)
  E = 1              (normalised energy)

Equations of motion:
  dr/dλ = f·pr,          dθ/dλ = pθ/r²,           dφ/dλ = L/(r² sin²θ)
  d(pr)/dλ = −E²/(r²f²) − pr²/r² + pθ²/r³ + L²/(r³sin²θ)
  d(pθ)/dλ = L² cosθ / (r² sin³θ)
where f = 1 − 2M/r.

Run from any directory:
    python scripts/plot_3d_schematic.py [output_path]
"""

import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D          # noqa
from scipy.integrate import solve_ivp
from matplotlib.lines import Line2D

# ── output ────────────────────────────────────────────────────────────────────
PAPER_DIR = Path(__file__).resolve().parents[2] / "paper"
OUT_PATH  = Path(sys.argv[1]) if len(sys.argv) > 1 else \
            PAPER_DIR / "kerr_3d_schematic.png"

# ── geometry (r_g = GM/c² = 1, M = 1) ────────────────────────────────────────
M        = 1.0
R_H      = 2.0 * M        # Schwarzschild horizon
R_SOURCE = 3.5             # UHE source ring
R_TORUS  = 10.0            # torus radial centre
SIG_T    = 4.0             # torus radial half-width (visual)
H_T      = 2.2             # torus half-height
R_CAM    = 52.0            # camera distance
TH_CAM   = np.radians(80.0)  # camera polar angle
PHI_CAM  = 0.0             # camera azimuth

PLOT_LIM = 11.0            # axes half-width for display (r_g)

# ── colours ───────────────────────────────────────────────────────────────────
C_CAP  = "#e05555"         # red   – captured
C_SRC  = "#ffffaa"         # pale-yellow – reaches source ring
C_TOR  = "#f0a040"         # amber – crosses torus (DIS)
C_FREE = "#70c8ff"         # sky-blue  – unattenuated
C_RING = "#ffee33"         # UHE ring
C_TSUF = (0.15, 0.50, 1.0) # torus surface

plt.rcParams.update({
    "font.family": "serif",
    "font.serif":  ["DejaVu Serif", "Times New Roman", "serif"],
    "font.size": 9,
})

# ─────────────────────────────────────────────────────────────────────────────
# 3D Schwarzschild geodesic from camera position
# ─────────────────────────────────────────────────────────────────────────────

def integrate_geodesic(alpha, beta, n=6000):
    """
    Integrate null geodesic starting at the camera.
    alpha, beta : image-plane coords (r_g).
    Returns (x, y, z, r) in Boyer-Lindquist-to-Cartesian conversion.
    """
    E = 1.0
    L = alpha * np.sin(TH_CAM)          # azimuthal angular momentum
    K = alpha**2 + beta**2              # Carter constant

    r0, th0, phi0 = R_CAM, TH_CAM, PHI_CAM
    f0  = 1.0 - 2.0*M/r0

    # Initial radial momentum (inward, null condition)
    pr0 = -np.sqrt(max(E**2/f0**2 - K/(r0**2 * f0), 0.0))

    # Initial polar momentum: pθ² = K − L²/sin²θ at the camera
    pth2 = K - (L/np.sin(th0))**2
    pth0 = -beta                         # sign: β>0 → toward pole → dθ<0 → pθ<0

    def rhs(lam, s):
        r, th, phi, pr, pth = s
        th   = np.clip(th, 1e-4, np.pi - 1e-4)
        sth  = np.sin(th);  cth = np.cos(th)
        f    = 1.0 - 2.0*M/r

        dr   = f * pr
        dth  = pth / r**2
        dphi = L / (r**2 * sth**2)
        dpr  = (-E**2/(r**2*f**2) - pr**2/r**2
                + pth**2/r**3 + L**2/(r**3*sth**2))
        dpth = L**2 * cth / (r**2 * sth**3)
        return [dr, dth, dphi, dpr, dpth]

    # Stop at horizon
    def ev_cap(lam, s):
        return s[0] - R_H * 1.04
    ev_cap.terminal  = True
    ev_cap.direction = -1

    y0  = [r0, th0, phi0, pr0, pth0]
    lam = np.linspace(0, 160, n)
    sol = solve_ivp(rhs, [0, 160], y0, t_eval=lam,
                    events=[ev_cap],
                    max_step=0.15, rtol=1e-8, atol=1e-10)

    r, th, phi = sol.y[0], sol.y[1], sol.y[2]
    x = r * np.sin(th) * np.cos(phi)
    y = r * np.sin(th) * np.sin(phi)
    z = r * np.cos(th)
    return x, y, z, r


def classify(r_arr, z_arr):
    rm = r_arr.min()
    if rm <= R_H * 1.06:
        return "cap"
    if rm <= R_SOURCE * 1.30:
        return "src"
    # Ray crosses torus only if it passes close to the equatorial plane
    # AND within the torus radial extent.
    in_torus = (np.abs(z_arr) < H_T * 1.5) & \
               (r_arr > R_TORUS - SIG_T) & \
               (r_arr < R_TORUS + SIG_T)
    if in_torus.any():
        return "tor"
    return "free"


CMAP = {"cap": C_CAP, "src": C_SRC, "tor": C_TOR, "free": C_FREE}
ALPH = {"cap": 0.95, "src": 0.97, "tor": 0.80, "free": 0.70}
LW   = {"cap": 1.6,  "src": 1.7,  "tor": 1.0,  "free": 0.85}

# ── one ray per fate — (alpha, beta, fate) ───────────────────────────────────
# fate is set explicitly (schematic, not a simulation result)
RAY_PARAMS = [
    ( 3.5,   0.0,  "cap"),   # captured  (b=3.5 < b_c≈5.2, equatorial)
    ( 5.5,   0.0,  "src"),   # reaches source ring (near photon sphere)
    ( 8.5,   5.0,  "tor"),   # crosses equatorial torus (DIS attenuation)
    ( 0.0,  10.0,  "free"),  # polar arc above torus (L=0, stays above equatorial plane)
]

# ─────────────────────────────────────────────────────────────────────────────
# Build figure
# ─────────────────────────────────────────────────────────────────────────────
fig = plt.figure(figsize=(7.0, 5.6), facecolor="black")
ax  = fig.add_subplot(111, projection="3d")
ax.set_facecolor("black")

# ── BH sphere ──────────────────────────────────────────────────────────────
N   = 50
u_b = np.linspace(0, 2*np.pi, N)
v_b = np.linspace(0,   np.pi, N//2)
xs  = R_H * np.outer(np.cos(u_b), np.sin(v_b))
ys  = R_H * np.outer(np.sin(u_b), np.sin(v_b))
zs  = R_H * np.outer(np.ones_like(u_b), np.cos(v_b))
ax.plot_surface(xs, ys, zs, color="black", zorder=6, shade=False, linewidth=0)
# photon-sphere glow ring
ph_eq = np.linspace(0, 2*np.pi, 300)
ax.plot(R_H*np.cos(ph_eq), R_H*np.sin(ph_eq), np.zeros(300),
        color="white", lw=1.2, alpha=0.75, zorder=8)
for inc in [35, 90, 145]:
    ti = np.radians(inc)
    ax.plot(R_H*np.cos(ph_eq),
            R_H*np.sin(ph_eq)*np.cos(ti),
            R_H*np.sin(ph_eq)*np.sin(ti),
            color="white", lw=0.4, alpha=0.25, zorder=5)

# ── Torus surface ──────────────────────────────────────────────────────────
phi_t = np.linspace(0, 2*np.pi, 100)
chi_t = np.linspace(-np.pi, np.pi, 50)
PP, CC = np.meshgrid(phi_t, chi_t)
Rt = R_TORUS + SIG_T*np.cos(CC)
xt = Rt*np.cos(PP);  yt = Rt*np.sin(PP);  zt = H_T*np.sin(CC)
ax.plot_surface(xt, yt, zt, color=C_TSUF, alpha=0.18,
                shade=True, linewidth=0, rstride=2, cstride=4, zorder=2)
for frac in (-1, 0, 1):
    r_rim = R_TORUS + frac*SIG_T*0.65
    ax.plot(r_rim*np.cos(phi_t), r_rim*np.sin(phi_t),
            np.full_like(phi_t, frac*H_T*0.65),
            color=(0.35, 0.65, 1.0), lw=0.5, alpha=0.28, zorder=3)

# ── UHE source ring ────────────────────────────────────────────────────────
phi_s = np.linspace(0, 2*np.pi, 512)
ax.plot(R_SOURCE*np.cos(phi_s), R_SOURCE*np.sin(phi_s), np.zeros(512),
        color=C_RING, lw=2.8, alpha=0.97, zorder=9)

# ── Geodesics ──────────────────────────────────────────────────────────────
print("Integrating geodesics...")
for alpha, beta, fate in RAY_PARAMS:
    try:
        x3, y3, z3, r_arr = integrate_geodesic(alpha, beta)
    except Exception as exc:
        print(f"  Skip (α={alpha}, β={beta}): {exc}")
        continue

    color = CMAP[fate]
    alpha_ = ALPH[fate]
    lw_    = LW[fate]

    # ── plot the in-region portion ──────────────────────────────────────────
    mask = r_arr < PLOT_LIM * 1.12
    if mask.sum() < 5:
        print(f"  (α={alpha:+5.1f}, β={beta:+5.1f})  fate={fate}  SKIPPED (no visible points)")
        continue

    ax.plot(x3[mask], y3[mask], z3[mask],
            color=color, lw=lw_, alpha=alpha_, zorder=7)

    # ── extend each ray to the figure edge in the camera direction ──────────
    # Find the FIRST point inside the visible region (entry from observer side).
    idx_entry = np.where(mask)[0][0]

    # Tangent direction at the entry point: use the first two visible points
    # for better precision; if only one visible point, fall back to cam_dir.
    if mask.sum() >= 2:
        idx2 = np.where(mask)[0][1]
        tang = np.array([x3[idx2]-x3[idx_entry],
                         y3[idx2]-y3[idx_entry],
                         z3[idx2]-z3[idx_entry]])
    else:
        tang = -np.array([np.sin(TH_CAM), 0.0, np.cos(TH_CAM)])

    tang_norm = np.linalg.norm(tang)
    if tang_norm < 1e-12:
        tang = -np.array([np.sin(TH_CAM), 0.0, np.cos(TH_CAM)])
        tang_norm = 1.0
    tang = tang / tang_norm

    # The extension goes BACKWARDS from the entry point (away from BH,
    # toward the observer), i.e., in the -tang direction.
    ext_dir = -tang
    P = np.array([x3[idx_entry], y3[idx_entry], z3[idx_entry]])

    # Find t_max so that P + t*ext_dir hits the box boundary.
    t_max = np.inf
    xlims = (-PLOT_LIM, PLOT_LIM)
    ylims = (-PLOT_LIM, PLOT_LIM)
    zlims = (-8.0,       8.0)
    for i, (lo, hi) in enumerate(zip(
            [xlims[0], ylims[0], zlims[0]],
            [xlims[1], ylims[1], zlims[1]])):
        d = ext_dir[i]
        if abs(d) > 1e-9:
            t_lo = (lo - P[i]) / d
            t_hi = (hi - P[i]) / d
            for t in (t_lo, t_hi):
                if t > 1e-6:
                    t_max = min(t_max, t)

    if np.isfinite(t_max):
        P_edge = P + t_max * ext_dir
        ax.plot([P[0], P_edge[0]],
                [P[1], P_edge[1]],
                [P[2], P_edge[2]],
                color=color, lw=lw_, alpha=alpha_ * 0.55,
                linestyle="--", zorder=6)

    print(f"  (α={alpha:+5.1f}, β={beta:+5.1f})  fate={fate}")

# ── Observer direction arrow ──────────────────────────────────────────────
# Arrow placed near the plot edge, pointing OUTWARD toward the observer.
# Rays enter from this direction so the arrow shows where they come from.
cam_dir = np.array([np.sin(TH_CAM), 0, np.cos(TH_CAM)])  # unit vector to observer
arrow_base = cam_dir * PLOT_LIM * 0.48   # start closer to BH
arrow_vec  = cam_dir * PLOT_LIM * 0.48   # longer arrow reaching near the edge
ax.quiver(*arrow_base, *arrow_vec,
          color="white", lw=2.0, alpha=0.90,
          arrow_length_ratio=0.18, zorder=10)
tip = arrow_base + arrow_vec
ax.text(tip[0]*1.03, tip[1]*1.03, tip[2] + 0.6,
        "Observer", color="white", fontsize=7.5,
        ha="center", va="bottom", zorder=11)

# ── Labels ─────────────────────────────────────────────────────────────────
ax.text(0, 0, -3.8, "Kerr BH",
        color="white", fontsize=8.5, ha="center", va="top", zorder=12)
ax.text(-R_SOURCE*0.5, -R_SOURCE*0.5, 2.5,
        "UHE source\nring", color=C_RING, fontsize=7,
        ha="right", va="bottom", zorder=12)
ax.text(-(R_TORUS + SIG_T + 1.0), -3.0, H_T + 1.0,
        "Dense\ntorus", color=(0.55, 0.80, 1.0), fontsize=7,
        ha="right", va="bottom", zorder=12)

# ── Legend ──────────────────────────────────────────────────────────────────
legend_elements = [
    Line2D([0],[0], color=C_CAP,  lw=1.4, label="Captured by BH"),
    Line2D([0],[0], color=C_SRC,  lw=1.4, label="Reaches source ring"),
    Line2D([0],[0], color=C_TOR,  lw=1.4, label="Crosses torus (DIS)"),
    Line2D([0],[0], color=C_FREE, lw=1.4, label="Free trajectory"),
]
ax.legend(handles=legend_elements, loc="upper left",
          fontsize=6.5, framealpha=0.20, facecolor="black",
          edgecolor="0.4", labelcolor="white", handlelength=1.8,
          borderpad=0.5)

# ── View / axes ─────────────────────────────────────────────────────────────
ax.set_xlim(-PLOT_LIM, PLOT_LIM)
ax.set_ylim(-PLOT_LIM, PLOT_LIM)
ax.set_zlim(-8, 8)
ax.set_box_aspect([1.0, 1.0, 0.45])
ax.set_axis_off()
ax.view_init(elev=28, azim=-50)

plt.tight_layout(pad=0.1)
fig.savefig(OUT_PATH, dpi=300, facecolor="black",
            bbox_inches="tight", pad_inches=0.06)
plt.close(fig)
print(f"Saved: {OUT_PATH}")
