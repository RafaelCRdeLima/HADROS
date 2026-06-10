#include "constants.hpp"
#include "radiative_transfer.hpp"
#include "ray.hpp"
#include "sigma_table.hpp"
#include "torus_profile.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>

namespace {

struct BinaryRayHeader {
    std::int32_t ray_id;
    std::int32_t pixel_i;
    std::int32_t pixel_j;
    std::int32_t captured;
    std::int32_t npoints;
    double alpha_rg;
    double beta_rg;
};

struct Metrics {
    double tau_sum = 0.0;
    double tau_max = 0.0;
    double psurv_sum = 0.0;
    double intensity_sum = 0.0;
    std::size_t rays = 0;
    std::size_t rays_reaching_torus = 0;
};

double legacy_gaussian_rho(
    double r_rg,
    double theta,
    double rho0,
    double r0,
    double sigma_r,
    double h_over_r
)
{
    constexpr double pi = 3.141592653589793238462643383279502884;
    const double delta = theta - 0.5 * pi;

    if (r_rg <= 4.0 || r_rg >= 18.0 || std::abs(delta) >= 0.45) {
        return 0.0;
    }

    const double radial =
        std::exp(-std::pow((r_rg - r0) / sigma_r, 2.0));

    const double vertical =
        std::exp(-std::pow(delta / h_over_r, 2.0));

    return rho0 * radial * vertical;
}

double zamo_g_from_point(const PathPoint& p)
{
    return 1.0 / std::max(p.redshift_factor, 1.0e-300);
}

Metrics accumulate_metrics(
    const std::string& cache_path,
    bool use_new_profile,
    double Enu_obs_GeV,
    double M_bh_msun,
    double rho0,
    double r0,
    double sigma_r,
    double h_over_r,
    double rho_floor,
    const SigmaTable& sigma,
    const UHECircularSourceParams& source
)
{
    std::ifstream in(cache_path, std::ios::binary);
    if (!in) {
        throw std::runtime_error("Could not open cache: " + cache_path);
    }

    std::int32_t magic = 0;
    std::int32_t version = 0;
    std::int32_t nx = 0;
    std::int32_t ny = 0;
    double cache_spin = 0.0;

    in.read(reinterpret_cast<char*>(&magic), sizeof(magic));
    in.read(reinterpret_cast<char*>(&version), sizeof(version));
    in.read(reinterpret_cast<char*>(&nx), sizeof(nx));
    in.read(reinterpret_cast<char*>(&ny), sizeof(ny));
    in.read(reinterpret_cast<char*>(&cache_spin), sizeof(cache_spin));

    if (magic != 0x4B47454F || version != 1) {
        throw std::runtime_error("Invalid geodesic cache: " + cache_path);
    }

    TorusProfile new_torus(
        rho0,
        r0,
        sigma_r,
        h_over_r,
        "gaussian",
        2.0,
        0.0,
        15.0 * constants::pi / 180.0,
        0.0,
        2.5,
        4.0,
        60.0,
        rho_floor
    );

    const double rg_cm = constants::rg_cm(M_bh_msun);
    Metrics metrics;

    while (true) {
        BinaryRayHeader header{};
        in.read(reinterpret_cast<char*>(&header), sizeof(header));
        if (!in) {
            break;
        }

        if (header.npoints < 0) {
            throw std::runtime_error("Negative point count in cache.");
        }

        double tau = 0.0;
        double intensity = 0.0;
        bool reaches_torus = false;

        for (int k = 0; k < header.npoints; ++k) {
            PathPoint p{};
            in.read(reinterpret_cast<char*>(&p), sizeof(p));
            if (!in) {
                throw std::runtime_error("Unexpected EOF in cache.");
            }

            const double raw_legacy =
                legacy_gaussian_rho(
                    p.r_rg,
                    p.theta,
                    rho0,
                    r0,
                    sigma_r,
                    h_over_r
                );

            const double rho =
                use_new_profile
                ? new_torus.rho(p.r_rg, p.theta)
                : raw_legacy;

            const double raw_new =
                new_torus.raw_rho(p.r_rg, p.theta);

            if ((use_new_profile ? raw_new : raw_legacy) > 0.0) {
                reaches_torus = true;
            }

            const double g = zamo_g_from_point(p);
            const double Enu_local_GeV =
                Enu_obs_GeV / std::max(g, 1.0e-300);

            const double j =
                radiative_transfer::emissivity_collapsar_ring(
                    p.r_rg,
                    p.theta,
                    Enu_local_GeV,
                    source
                );

            const double dl_cm = p.dl_rg * rg_cm;

            intensity +=
                std::pow(g, 3.0)
                * j
                * std::exp(-tau)
                * dl_cm;

            if (rho <= 0.0 ||
                Enu_local_GeV < sigma.Emin() ||
                Enu_local_GeV > sigma.Emax()) {
                continue;
            }

            tau +=
                (rho / constants::m_u_g)
                * sigma.sigma_cm2(Enu_local_GeV)
                * dl_cm;
        }

        metrics.rays += 1;
        metrics.tau_sum += tau;
        metrics.tau_max = std::max(metrics.tau_max, tau);
        metrics.psurv_sum += std::exp(-tau);
        metrics.intensity_sum += intensity;
        if (reaches_torus) {
            metrics.rays_reaching_torus += 1;
        }
    }

    return metrics;
}

double relative_difference(double value, double reference)
{
    if (reference == 0.0) {
        return value == 0.0 ? 0.0 : std::numeric_limits<double>::infinity();
    }

    return (value - reference) / reference;
}

void write_metric(
    std::ofstream& out,
    const std::string& name,
    double legacy,
    double current
)
{
    out << name << " "
        << std::setprecision(16)
        << legacy << " "
        << current << " "
        << relative_difference(current, legacy) << "\n";
}

} // namespace

int main(int argc, char* argv[])
{
    std::string cache_path = "output/rays/kerr_geodesics_small.bin";
    std::string report_path =
        "plots/validation_density_backgrounds/gaussian_backward_compatibility.txt";

    if (argc > 1) {
        cache_path = argv[1];
    }
    if (argc > 2) {
        report_path = argv[2];
    }

    const double Enu_obs_GeV = 1.0e5;
    const double M_bh_msun = 3.0;
    const double rho0 = 1.0e-2;
    const double r0 = 10.0;
    const double sigma_r = 5.0;
    const double h_over_r = 0.25;
    const double rho_floor = 1.0e-99;

    SigmaTable sigma("data/sigma/sigma_nuN_CC_GBW.dat");

    UHECircularSourceParams source;
    source.r_center_rg = 3.5;
    source.sigma_r_rg = 1.0;
    source.theta_width_rad = 15.0 * constants::pi / 180.0;
    source.powerlaw = 2.0;
    source.emax_GeV = 1.0e12;
    source.norm = 1.0;

    const Metrics legacy =
        accumulate_metrics(
            cache_path,
            false,
            Enu_obs_GeV,
            M_bh_msun,
            rho0,
            r0,
            sigma_r,
            h_over_r,
            rho_floor,
            sigma,
            source
        );

    const Metrics current =
        accumulate_metrics(
            cache_path,
            true,
            Enu_obs_GeV,
            M_bh_msun,
            rho0,
            r0,
            sigma_r,
            h_over_r,
            rho_floor,
            sigma,
            source
        );

    std::ofstream out(report_path);
    if (!out) {
        std::cerr << "Could not open report: " << report_path << "\n";
        return 1;
    }

    out << "# Gaussian backward-compatibility validation\n";
    out << "# reference legacy_gaussian_formula\n";
    out << "# current TorusProfile_gaussian\n";
    out << "# cache_path " << cache_path << "\n";
    out << "# metric legacy current relative_difference\n";

    write_metric(
        out,
        "mean_optical_depth",
        legacy.tau_sum / static_cast<double>(legacy.rays),
        current.tau_sum / static_cast<double>(current.rays)
    );

    write_metric(out, "maximum_optical_depth", legacy.tau_max, current.tau_max);

    write_metric(
        out,
        "mean_survival_probability",
        legacy.psurv_sum / static_cast<double>(legacy.rays),
        current.psurv_sum / static_cast<double>(current.rays)
    );

    write_metric(
        out,
        "total_image_intensity",
        legacy.intensity_sum,
        current.intensity_sum
    );

    out << "rays_reaching_torus "
        << legacy.rays_reaching_torus << " "
        << current.rays_reaching_torus << " "
        << relative_difference(
            static_cast<double>(current.rays_reaching_torus),
            static_cast<double>(legacy.rays_reaching_torus)
        ) << "\n";

    std::cout << "Saved: " << report_path << "\n";
    return 0;
}
