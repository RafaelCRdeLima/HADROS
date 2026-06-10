#include "schwarzschild_raytracer.hpp"

#include <fstream>
#include <iostream>
#include <iomanip>
#include <cmath>

int main()
{
    SchwarzschildRayTracer tracer(
        80.0,   // observer distance r_obs / r_g
        120.0,  // outer boundary r_max / r_g
        2.0,    // horizon radius / r_g
        0.02    // integration step
    );

    const int Npix = 9;
    const double fov_rg = 20.0;

    std::ofstream out("output/rays/schwarzschild_rays.dat");

    out << "# ray_id pixel_i pixel_j alpha_rg beta_rg "
        << "x_rg z_rg r_rg theta dl_rg redshift captured\n";

    int ray_id = 0;

    for (int i = 0; i < Npix; ++i) {
        double alpha =
            -0.5 * fov_rg + fov_rg * i / (Npix - 1);

        for (int j = 0; j < Npix; ++j) {
            double beta =
                -0.5 * fov_rg + fov_rg * j / (Npix - 1);

            RayPath ray =
                tracer.trace_ray(alpha, beta, i, j);

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

    std::cout << "Saved: output/rays/schwarzschild_rays.dat\n";

    return 0;
}