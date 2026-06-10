#include "schwarzschild_camera.hpp"

#include <fstream>
#include <iostream>
#include <iomanip>
#include <cmath>

int main()
{
    SchwarzschildCamera camera(
        80.0,      // r_obs / r_g
        M_PI/2.0,  // observer at equatorial plane
        25.0,      // field of view in degrees
        15,        // nx
        15,        // ny
        120.0,     // r_max / r_g
        2.0,       // horizon / r_g
        0.02       // step
    );

    std::ofstream out("output/rays/schwarzschild_camera_rays.dat");

    out << "# ray_id pixel_i pixel_j alpha beta "
        << "x_rg z_rg r_rg theta dl_rg redshift captured\n";

    int ray_id = 0;

    for (int i = 0; i < camera.nx(); ++i) {
        for (int j = 0; j < camera.ny(); ++j) {

            RayPath ray = camera.trace_pixel(i, j);

            for (const auto& p : ray.points) {
                out << std::scientific << std::setprecision(8)
                    << ray_id << " "
                    << ray.pixel_i << " "
                    << ray.pixel_j << " "
                    << ray.alpha_rg << " "
                    << ray.beta_rg << " "
                    << p.x_rg << " "
                    << p.z_rg << " "
                    << p.r_rg << " "
                    << p.theta << " "
                    << p.dl_rg << " "
                    << p.redshift_factor << " "
                    << ray.captured << "\n";
            }

            out << "\n";
            ++ray_id;
        }
    }

    std::cout << "Saved: output/rays/schwarzschild_camera_rays.dat\n";

    return 0;
}