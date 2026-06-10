import struct
import math
import os
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import numpy as np
import matplotlib.pyplot as plt

torus = np.loadtxt("output/profiles/torus_synthetic_grid.dat")


def _env_int(name, default):
    value = os.environ.get(name)
    if value is None:
        return default
    try:
        parsed = int(value)
    except ValueError as exc:
        raise RuntimeError(f"{name} must be an integer, got {value!r}.") from exc
    if parsed <= 0:
        raise RuntimeError(f"{name} must be positive, got {parsed}.")
    return parsed


def _env_float(name, default):
    value = os.environ.get(name)
    if value is None:
        return default
    try:
        parsed = float(value)
    except ValueError as exc:
        raise RuntimeError(f"{name} must be a number, got {value!r}.") from exc
    if parsed <= 0:
        raise RuntimeError(f"{name} must be positive, got {parsed}.")
    return parsed


def load_cached_geodesics(filename, max_rays=5, max_points_per_ray=1200):
    path = Path(filename)
    if not path.exists():
        raise FileNotFoundError(
            f"Missing cached geodesics: {filename}. Run `make kerr-geodesics` first."
        )

    file_header = struct.Struct("=iiiid")
    ray_header = struct.Struct("=iiiii4xdd")
    point_size = struct.Struct("=ddddddd").size

    rays = []

    with path.open("rb") as f:
        raw = f.read(file_header.size)
        if len(raw) != file_header.size:
            raise RuntimeError("Invalid geodesic cache: truncated file header.")

        magic, version, nx, ny, _a_spin = file_header.unpack(raw)
        if magic != 0x4B47454F:
            raise RuntimeError("Invalid geodesic cache: bad magic number.")
        if version != 1:
            raise RuntimeError(f"Unsupported geodesic cache version: {version}.")

        total_rays = nx * ny
        ray_stride = max(1, math.ceil(total_rays / max_rays))

        while True:
            raw = f.read(ray_header.size)
            if not raw:
                break
            if len(raw) != ray_header.size:
                raise RuntimeError("Invalid geodesic cache: truncated ray header.")

            (
                ray_id,
                pixel_i,
                pixel_j,
                captured,
                npoints,
                alpha_rg,
                beta_rg,
            ) = ray_header.unpack(raw)

            if npoints < 0:
                raise RuntimeError("Invalid geodesic cache: negative point count.")

            if ray_id % ray_stride != 0:
                f.seek(npoints * point_size, 1)
                continue

            payload_size = npoints * point_size
            payload = f.read(payload_size)
            if len(payload) != payload_size:
                raise RuntimeError("Invalid geodesic cache: truncated point data.")
            if npoints == 0:
                continue

            point_stride = max(1, math.ceil(npoints / max_points_per_ray))
            points = np.frombuffer(payload, dtype=np.float64).reshape(-1, 7)
            rays.append(points[::point_stride, 2:5].copy())

    if not rays:
        raise RuntimeError("No rays were loaded from the geodesic cache.")

    return rays


rays = load_cached_geodesics(
    "output/rays/kerr_geodesics.bin",
    max_rays=_env_int("KERR_TORUS_MAX_RAYS", 5),
    max_points_per_ray=_env_int("KERR_TORUS_MAX_POINTS_PER_RAY", 1200),
)

r = torus[:, 0]
th = torus[:, 1]
rho = torus[:, 2]

# Keep only dense torus regions
mask = rho > 0.1 * rho.max()

r = r[mask]
th = th[mask]
rho = rho[mask]

# Reduce number of torus points
r = r[::5]
th = th[::5]
rho = rho[::5]

phis = np.linspace(0, 2*np.pi, 96)

X, Y, Z = [], [], []

for ph in phis:
    X.append(r * np.sin(th) * np.cos(ph))
    Y.append(r * np.sin(th) * np.sin(ph))
    Z.append(r * np.cos(th))

X = np.concatenate(X)
Y = np.concatenate(Y)
Z = np.concatenate(Z)

fig = plt.figure(figsize=(8, 7), facecolor="black")
ax = fig.add_subplot(111, projection="3d")
ax.set_facecolor("black")

source_r_rg = _env_float("SOURCE_R_RG", 3.5)
source_phi = np.linspace(0, 2*np.pi, 512)
source_x = source_r_rg * np.cos(source_phi)
source_y = source_r_rg * np.sin(source_phi)
source_z = np.zeros_like(source_phi)

# Lightweight torus
ax.scatter(
    X, Y, Z,
    s=1,
    alpha=0.03
)

for ray in rays:
    ax.plot(
        ray[:, 0],
        ray[:, 1],
        ray[:, 2],
        lw=0.8,
        alpha=0.5
    )

ax.plot(
    source_x,
    source_y,
    source_z,
    color="yellow",
    lw=2.4,
    alpha=0.95
)

ax.set_xlim(-13, 23)
ax.set_ylim(-19, 19)
ax.set_zlim(-15, 15)

ax.set_box_aspect([1.4, 1, 1])
ax.set_axis_off()

plt.tight_layout()

plt.savefig(
    "plots/kerr_rays_on_torus_3d_light.png",
    dpi=200,
    facecolor=fig.get_facecolor(),
    bbox_inches="tight",
    pad_inches=0
)
