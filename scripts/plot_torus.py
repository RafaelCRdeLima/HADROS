import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt("output/profiles/torus_synthetic_grid.dat")

r = data[:, 0]
th = data[:, 1]
rho = data[:, 2]
T = data[:, 3]
Ye = data[:, 4]

x = r * np.sin(th)
z = r * np.cos(th)

plt.figure(figsize=(7, 6))
plt.tricontourf(x, z, np.log10(rho), levels=80)
plt.colorbar(label=r"$\log_{10}(\rho\,[\mathrm{g\,cm^{-3}}])$")
plt.title("Synthetic Torus Density Profile")
plt.xlabel(r"$x/r_g$")
plt.ylabel(r"$z/r_g$")
plt.gca().set_aspect("equal")
plt.tight_layout()
plt.savefig("plots/torus_density.png", dpi=200)
plt.show()

plt.figure(figsize=(7, 6))
plt.tricontourf(x, z, T, levels=80)
plt.colorbar(label=r"$T\,[\mathrm{MeV}]$")
plt.title("Synthetic Torus Temperature Profile")
plt.xlabel(r"$x/r_g$")
plt.ylabel(r"$z/r_g$")
plt.gca().set_aspect("equal")
plt.tight_layout()
plt.savefig("plots/torus_temperature.png", dpi=200)
plt.show()

plt.figure(figsize=(7, 6))
plt.tricontourf(x, z, Ye, levels=80)
plt.colorbar(label=r"$Y_e$")
plt.title("Synthetic Torus Electron Fraction Profile")
plt.xlabel(r"$x/r_g$")
plt.ylabel(r"$z/r_g$")
plt.gca().set_aspect("equal")
plt.tight_layout()
plt.savefig("plots/torus_Ye.png", dpi=200)
plt.show()