import numpy as np
import matplotlib.pyplot as plt

torus = np.loadtxt("output/profiles/torus_synthetic_grid.dat")
rays = np.loadtxt("output/rays/schwarzschild_camera_rays.dat")

r = torus[:, 0]
th = torus[:, 1]
rho = torus[:, 2]

xt = r * np.sin(th)
zt = r * np.cos(th)

ray_id = rays[:, 0].astype(int)
x = rays[:, 5]
z = rays[:, 6]

plt.figure(figsize=(8, 7))

plt.tricontourf(
    xt, zt, np.log10(rho),
    levels=80,
    alpha=0.75
)

plt.colorbar(label=r"$\log_{10}(\rho\,[\mathrm{g\,cm^{-3}}])$")

for rid in np.unique(ray_id):
    mask = ray_id == rid
    plt.plot(x[mask], z[mask], color="black", lw=0.5, alpha=0.45)

horizon = plt.Circle((0, 0), 2.0, color="black", fill=True)
plt.gca().add_patch(horizon)

plt.title("Schwarzschild Camera Rays Through a Synthetic Torus")
plt.xlabel(r"$x/r_g$")
plt.ylabel(r"$z/r_g$")
plt.gca().set_aspect("equal")
plt.xlim(-20, 85)
plt.ylim(-35, 35)
plt.tight_layout()
plt.savefig("plots/camera_rays_on_torus.png", dpi=200)
plt.show()