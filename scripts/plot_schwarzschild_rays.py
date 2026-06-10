import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt("output/rays/schwarzschild_camera_rays.dat")

ray_id = data[:, 0].astype(int)
x = data[:, 5]
z = data[:, 6]
captured = data[:, 11].astype(int)

plt.figure(figsize=(7, 7))

for rid in np.unique(ray_id):
    mask = ray_id == rid
    lw = 1.0 if captured[mask][0] else 0.6
    plt.plot(x[mask], z[mask], lw=lw)

horizon = plt.Circle((0, 0), 2.0, fill=False, lw=2)
plt.gca().add_patch(horizon)

plt.title("Backward Ray Tracing from a Schwarzschild Camera")
plt.xlabel(r"$x/r_g$")
plt.ylabel(r"$z/r_g$")
plt.gca().set_aspect("equal")
plt.xlim(-5, 85)
plt.ylim(-35, 35)
plt.tight_layout()
plt.savefig("plots/schwarzschild_camera_rays.png", dpi=200)
plt.show()