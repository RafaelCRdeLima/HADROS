import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt("output/tau/tau_straight.dat")

b = data[:, 0]
tau = data[:, 1]
P = data[:, 2]

plt.figure(figsize=(7, 5))
plt.semilogy(b, tau)
plt.title("Optical Depth for Straight Rays")
plt.xlabel(r"Impact parameter $b/r_g$")
plt.ylabel(r"$\tau(E_\nu,b)$")
plt.tight_layout()
plt.savefig("plots/tau_straight.png", dpi=200)
plt.show()

plt.figure(figsize=(7, 5))
plt.plot(b, P)
plt.title("Survival Probability for Straight Rays")
plt.xlabel(r"Impact parameter $b/r_g$")
plt.ylabel(r"$P_{\rm surv}=\exp[-\tau]$")
plt.ylim(-0.05, 1.05)
plt.tight_layout()
plt.savefig("plots/survival_straight.png", dpi=200)
plt.show()