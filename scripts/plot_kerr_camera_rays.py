import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt("output/rays/kerr_camera_rays.dat")

ray_id = data[:, 0].astype(int)

x = data[:, 5]
z = data[:, 6]

captured = data[:, 11].astype(int)

plt.figure(figsize=(8, 8))

for rid in np.unique(ray_id):

    mask = ray_id == rid

    if captured[mask][0]:
        lw = 1.2
        alpha = 0.9
    else:
        lw = 0.6
        alpha = 0.35

    plt.plot(
        x[mask],
        z[mask],
        lw=lw,
        alpha=alpha
    )

# Kerr horizon (approximate projection)
horizon = plt.Circle(
    (0, 0),
    1.43589,
    color="black",
    fill=True
)

plt.gca().add_patch(horizon)

plt.title("Backward Kerr Ray Tracing")

plt.xlabel(r"$x/r_g$")
plt.ylabel(r"$z/r_g$")

plt.gca().set_aspect("equal")

plt.xlim(-40, 90)
plt.ylim(-50, 50)

plt.tight_layout()

plt.savefig(
    "plots/kerr_camera_rays.png",
    dpi=200
)

plt.show()