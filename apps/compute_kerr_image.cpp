#include "kerr_camera.hpp"
#include "sigma_table.hpp"
#include "torus_profile.hpp"
#include "radiative_transfer.hpp"

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <cmath>

#include <vector>
#include <sstream>
#ifdef _OPENMP
#include <omp.h>
#endif

int main(int argc, char* argv[])
{
    double Enu_obs_GeV = 1.0e9;

    if (argc > 1) {
        Enu_obs_GeV = std::atof(argv[1]);
    }

    if (argc > 2) {
        a_spin = std::atof(argv[2]);
    }

    const double M_bh_msun = 3.0;

KerrCamera camera(
    a_spin,                 // parâmetro de rotação de Kerr a = J/M^2
    80.0,                   // raio do observador [r_g]
    M_PI * 72 / 180.0,      // inclinação do observador theta_obs [rad]
    18.0,                   // campo de visão (FOV) da câmera [graus]
    120,                    // número de pixels no eixo alpha
    120,                    // número de pixels no eixo beta
    120.0,                  // raio máximo da integração [r_g]
    0.001                   // passo máximo do integrador adaptativo RKF45
);

SigmaTable sigma(
    "data/sigma/sigma_nuN_CC_GBW.dat" // tabela de seção de choque DIS neutrino-núcleon
);

TorusProfile torus(
    1.0e-2,                 // densidade máxima do toro rho0 [g/cm^3]
    10.0,                   // raio central do toro r0 [r_g]
    5.0,                    // largura radial do toro sigma_r [r_g]
    0.18                    // H/R (espessura geométrica do toro)
);
    std::ofstream out("output/images/kerr_image.dat");

    out << "# i j alpha beta tau P_surv I_obs captured\n";

    const int nx = camera.nx();
const int ny = camera.ny();

std::vector<std::string> lines(nx * ny);

#pragma omp parallel for collapse(2) schedule(dynamic)
for (int i = 0; i < nx; ++i) {
        for (int j = 0; j < ny; ++j) {

            RayPath ray = camera.trace_pixel(i, j);

            RTResult rt =
                radiative_transfer::integrate_kerr_ray(
                    ray,
                    Enu_obs_GeV,
                    M_bh_msun,
                    torus,
                    sigma
                );

            std::ostringstream ss;

            ss << std::scientific
               << std::setprecision(8)
               << i << " "
               << j << " "
               << ray.alpha_rg << " "
               << ray.beta_rg << " "
               << rt.tau << " "
               << rt.P_surv << " "
               << rt.I_obs << " "
               << ray.captured << "\n";

            lines[i * ny + j] = ss.str();
        }
    }

    for (const auto& line : lines) {
        out << line;
    }

    std::cout << "Saved: output/images/kerr_image.dat\n";
    return 0;
}