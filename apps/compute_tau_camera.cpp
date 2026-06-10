#include "schwarzschild_camera.hpp"
#include "sigma_table.hpp"
#include "torus_profile.hpp"
#include "optical_depth.hpp"

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <cmath>

int main(int argc, char* argv[])
{
    double Enu_inf_GeV = 1.0e9;

    if (argc > 1) {
        Enu_inf_GeV = std::atof(argv[1]);
    }

    SchwarzschildCamera camera(
        80.0,
        M_PI/2.0,
        25.0,
        80,
        80,
        120.0,
        2.0,
        0.02
    );

    SigmaTable sigma("data/sigma/sigma_nuN_CC_GBW.dat");
    TorusProfile torus(
        1.0e3,  // rho0 g/cm^3
        10.0,   // r0
        5.0,    // sigma_r
        0.30    // H/R
    );
    std::ofstream out("output/tau/tau_camera.dat");

    out << "# i j alpha beta tau P_surv captured\n";

    for (int i = 0; i < camera.nx(); ++i) {
        for (int j = 0; j < camera.ny(); ++j) {

            RayPath ray = camera.trace_pixel(i, j);

            double tau =
                optical_depth::tau_along_ray(
                    ray,
                    Enu_inf_GeV,
                    torus,
                    sigma
                );

            double P =
                optical_depth::survival_probability(tau);

            out << std::scientific
                << std::setprecision(8)
                << i << " "
                << j << " "
                << ray.alpha_rg << " "
                << ray.beta_rg << " "
                << tau << " "
                << P << " "
                << ray.captured << "\n";
        }
    }

    std::cout << "Saved: output/tau/tau_camera.dat\n";
    return 0;
}