import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt("output/tau/tau_camera.dat")

i = data[:, 0]
j = data[:, 1]
alpha = data[:, 2]
beta = data[:, 3]
tau = data[:, 4]
P = data[:, 5]
captured = data[:, 6]

plt.figure(figsize=(7, 6))
plt.tricontourf(alpha, beta, np.log10(tau + 1e-300), levels=80)
plt.colorbar(label=r"$\log_{10}\tau$")
plt.title("DIS Optical Depth Map")
plt.xlabel(r"Camera coordinate $\alpha$")
plt.ylabel(r"Camera coordinate $\beta$")
plt.gca().set_aspect("equal")
plt.tight_layout()
plt.savefig("plots/tau_camera_map.png", dpi=200)
plt.show()

plt.figure(figsize=(7, 6))
plt.tricontourf(alpha, beta, P, levels=80)
plt.colorbar(label=r"$P_{\rm surv}=\exp(-\tau)$")
plt.title("UHE Neutrino Survival Probability Map")
plt.xlabel(r"Camera coordinate $\alpha$")
plt.ylabel(r"Camera coordinate $\beta$")
plt.gca().set_aspect("equal")
plt.tight_layout()
plt.savefig("plots/survival_camera_map.png", dpi=200)
plt.show()

Iobs = P.copy()
Iobs[captured > 0] = 0.0

plt.figure(figsize=(7, 6))

plt.tricontourf(
    alpha,
    beta,
    np.log10(Iobs + 1e-300),
    levels=80
)

plt.colorbar(label=r"$\log_{10}(I_{\rm obs}/I_0)$")

plt.title("UHE Neutrino Image with Black Hole Shadow")

plt.xlabel(r"Camera coordinate $\alpha$")
plt.ylabel(r"Camera coordinate $\beta$")

plt.gca().set_aspect("equal")

plt.tight_layout()

plt.savefig(
    "plots/uhe_neutrino_shadow_log_image.png",
    dpi=200
)

plt.show()