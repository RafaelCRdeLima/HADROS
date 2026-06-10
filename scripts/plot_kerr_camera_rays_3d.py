import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt("output/rays/kerr_camera_rays.dat")

ray_id = data[:, 0].astype(int)
x = data[:, 5]
y = data[:, 6]
z = data[:, 7]

fig = plt.figure(figsize=(8, 8))
ax = fig.add_subplot(111, projection="3d")

for rid in np.unique(ray_id):
    mask = ray_id == rid
    ax.plot(x[mask], y[mask], z[mask], lw=0.7, alpha=0.55)

ax.set_title("3D Backward Kerr Ray Tracing")
ax.set_xlabel(r"$x/r_g$")
ax.set_ylabel(r"$y/r_g$")
ax.set_zlabel(r"$z/r_g$")

ax.set_box_aspect([1, 1, 1])
plt.tight_layout()
plt.savefig("plots/kerr_camera_rays_3d.png", dpi=200)
plt.show()