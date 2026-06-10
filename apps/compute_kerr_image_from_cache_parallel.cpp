#include "ray.hpp"
#include "sigma_table.hpp"
#include "torus_profile.hpp"
#include "radiative_transfer.hpp"

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cstdlib>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

int main(int argc, char* argv[])
{
    double Enu_obs_GeV = 1.0e9;
    double a_spin = 0.9;

    if (argc > 1) {
        Enu_obs_GeV = std::atof(argv[1]);
    }

    if (argc > 2) {
        a_spin = std::atof(argv[2]);
    }

    const double M_bh_msun = 3.0;

    SigmaTable sigma("data/sigma/sigma_nuN_CC_GBW.dat");

    TorusProfile torus(
        1.0e-2,
        10.0,
        5.0,
        0.18
    );

    std::ifstream in("output/rays/kerr_geodesics.dat");

    if (!in) {
        std::cerr << "Could not open output/rays/kerr_geodesics.dat\n";
        return 1;
    }

    std::string line;

    int current_ray_id = -1;

    RayPath ray;
    ray.pixel_i = -1;
    ray.pixel_j = -1;
    ray.a_bh = a_spin;

    std::vector<RayPath> rays;

    auto store_ray = [&]() {
        if (ray.pixel_i < 0 || ray.pixel_j < 0) {
            return;
        }

        rays.push_back(ray);
    };

    while (std::getline(in, line)) {

        if (line.empty() || line[0] == '#') {
            continue;
        }

        std::istringstream iss(line);

        int ray_id;
        int i, j;
        double alpha, beta;
        double x, y, z, r, theta, dl, redshift;
        int captured;

        iss >> ray_id
            >> i >> j
            >> alpha >> beta
            >> x >> y >> z
            >> r >> theta >> dl >> redshift
            >> captured;

        if (!iss) {
            continue;
        }

        if (current_ray_id == -1) {
            current_ray_id = ray_id;

            ray = RayPath{};
            ray.a_bh = a_spin;
            ray.pixel_i = i;
            ray.pixel_j = j;
            ray.alpha_rg = alpha;
            ray.beta_rg = beta;
            ray.captured = static_cast<bool>(captured);
        }

        if (ray_id != current_ray_id) {
            store_ray();

            current_ray_id = ray_id;

            ray = RayPath{};
            ray.a_bh = a_spin;
            ray.pixel_i = i;
            ray.pixel_j = j;
            ray.alpha_rg = alpha;
            ray.beta_rg = beta;
            ray.captured = static_cast<bool>(captured);
        }

        PathPoint p;
        p.x_rg = x;
        p.y_rg = y;
        p.z_rg = z;
        p.r_rg = r;
        p.theta = theta;
        p.dl_rg = dl;
        p.redshift_factor = redshift;

        ray.points.push_back(p);
    }

    store_ray();

    std::vector<RTResult> results(rays.size());

    for (std::size_t k = 0; k < results.size(); ++k) {
        results[k].tau = 0.0;
        results[k].P_surv = 1.0;
        results[k].I_obs = 0.0;
    }

#ifdef _OPENMP
    std::cout << "OpenMP threads: " << omp_get_max_threads() << "\n";
#endif

#pragma omp parallel for schedule(dynamic)
    for (int k = 0; k < static_cast<int>(rays.size()); ++k) {
        if (!rays[k].points.empty()) {
            results[k] = radiative_transfer::integrate_kerr_ray(
                rays[k],
                Enu_obs_GeV,
                M_bh_msun,
                torus,
                sigma
            );
        }
    }

    std::ofstream out("output/images/kerr_image.dat");

    if (!out) {
        std::cerr << "Could not open output/images/kerr_image.dat\n";
        return 1;
    }

    out << "# i j alpha beta tau P_surv I_obs captured\n";

    for (std::size_t k = 0; k < rays.size(); ++k) {
        const RayPath& ray_out = rays[k];
        const RTResult& rt = results[k];

        out << std::scientific
            << std::setprecision(8)
            << ray_out.pixel_i << " "
            << ray_out.pixel_j << " "
            << ray_out.alpha_rg << " "
            << ray_out.beta_rg << " "
            << rt.tau << " "
            << rt.P_surv << " "
            << rt.I_obs << " "
            << ray_out.captured << "\n";
    }

    std::cout << "Loaded rays: " << rays.size() << "\n";
    std::cout << "Saved: output/images/kerr_image.dat\n";

    return 0;
}
