#include "radiative_transfer.hpp"
#include "torus_profile.hpp"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace {

struct Args {
    int nx = 64;
    int ny = 64;
    int nz = 64;
    double box_rg = 80.0;
    double Enu_GeV = 1.0e5;
    double rho0 = 1.0e-2;
    double r0 = 10.0;
    double sigma_r = 5.0;
    double H_over_R = 0.25;
    std::string density_profile = "gaussian";
    double radial_power = 2.0;
    double funnel_depletion = 0.0;
    double funnel_theta_deg = 15.0;
    double envelope_rho0 = 0.0;
    double envelope_alpha = 2.5;
    double r_min = 4.0;
    double r_max = 60.0;
    double rho_floor = 1.0e-99;
    UHESourceParams source;
    std::string output_path = "output/paraview/bh_torus_fields.vtk";
};

double arg_double(char* argv[], int argc, int index, double fallback)
{
    return argc > index ? std::atof(argv[index]) : fallback;
}

int arg_int(char* argv[], int argc, int index, int fallback)
{
    return argc > index ? std::atoi(argv[index]) : fallback;
}

std::string arg_string(char* argv[], int argc, int index, const std::string& fallback)
{
    return argc > index ? std::string(argv[index]) : fallback;
}

double finite_or_zero(double value)
{
    return std::isfinite(value) ? value : 0.0;
}

void write_scalar(
    std::ofstream& out,
    const std::string& name,
    const std::vector<double>& values
)
{
    out << "SCALARS " << name << " double 1\n";
    out << "LOOKUP_TABLE default\n";

    for (double value : values) {
        out << std::scientific << std::setprecision(10)
            << finite_or_zero(value) << "\n";
    }
}

} // namespace

int main(int argc, char* argv[])
{
    Args args;
    args.nx = arg_int(argv, argc, 1, args.nx);
    args.ny = arg_int(argv, argc, 2, args.ny);
    args.nz = arg_int(argv, argc, 3, args.nz);
    args.box_rg = arg_double(argv, argc, 4, args.box_rg);
    args.Enu_GeV = arg_double(argv, argc, 5, args.Enu_GeV);
    args.rho0 = arg_double(argv, argc, 6, args.rho0);
    args.r0 = arg_double(argv, argc, 7, args.r0);
    args.sigma_r = arg_double(argv, argc, 8, args.sigma_r);
    args.H_over_R = arg_double(argv, argc, 9, args.H_over_R);
    args.density_profile = arg_string(argv, argc, 10, args.density_profile);
    args.radial_power = arg_double(argv, argc, 11, args.radial_power);
    args.funnel_depletion = arg_double(argv, argc, 12, args.funnel_depletion);
    args.funnel_theta_deg = arg_double(argv, argc, 13, args.funnel_theta_deg);
    args.envelope_rho0 = arg_double(argv, argc, 14, args.envelope_rho0);
    args.envelope_alpha = arg_double(argv, argc, 15, args.envelope_alpha);
    args.r_min = arg_double(argv, argc, 16, args.r_min);
    args.r_max = arg_double(argv, argc, 17, args.r_max);
    args.rho_floor = arg_double(argv, argc, 18, args.rho_floor);
    args.source.r_center_rg = arg_double(argv, argc, 19, args.source.r_center_rg);
    args.source.sigma_r_rg = arg_double(argv, argc, 20, args.source.sigma_r_rg);
    args.source.theta_width_rad =
        arg_double(argv, argc, 21, args.source.theta_width_rad * 180.0 / 3.14159265358979323846)
        * 3.14159265358979323846 / 180.0;
    args.source.powerlaw = arg_double(argv, argc, 22, args.source.powerlaw);
    args.source.emax_GeV = arg_double(argv, argc, 23, args.source.emax_GeV);
    args.source.norm = arg_double(argv, argc, 24, args.source.norm);
    args.source.model_name = arg_string(argv, argc, 25, args.source.model_name);
    args.source.model = radiative_transfer::parse_uhe_source_model(args.source.model_name);
    args.source.funnel_theta_rad =
        arg_double(argv, argc, 26, args.source.funnel_theta_rad * 180.0 / 3.14159265358979323846)
        * 3.14159265358979323846 / 180.0;
    args.source.density_power_q = arg_double(argv, argc, 27, args.source.density_power_q);
    args.source.radial_power_s = arg_double(argv, argc, 28, args.source.radial_power_s);
    args.source.gradient_dr_rg = arg_double(argv, argc, 29, args.source.gradient_dr_rg);
    args.source.gradient_dtheta_rad =
        arg_double(argv, argc, 30, args.source.gradient_dtheta_rad * 180.0 / 3.14159265358979323846)
        * 3.14159265358979323846 / 180.0;
    args.source.rho_ref_gcm3 = arg_double(argv, argc, 31, args.source.rho_ref_gcm3);
    args.source.cutoff_min = arg_double(argv, argc, 32, args.source.cutoff_min);
    args.source.cutoff_max = arg_double(argv, argc, 33, args.source.cutoff_max);
    args.output_path = arg_string(argv, argc, 34, args.output_path);

    if (args.nx <= 1 || args.ny <= 1 || args.nz <= 1 || args.box_rg <= 0.0) {
        std::cerr << "Invalid ParaView grid parameters.\n";
        return 1;
    }

    TorusProfile torus(
        args.rho0,
        args.r0,
        args.sigma_r,
        args.H_over_R,
        args.density_profile,
        args.radial_power,
        args.funnel_depletion,
        args.funnel_theta_deg * 3.14159265358979323846 / 180.0,
        args.envelope_rho0,
        args.envelope_alpha,
        args.r_min,
        args.r_max,
        args.rho_floor
    );

    const std::size_t npoints =
        static_cast<std::size_t>(args.nx)
        * static_cast<std::size_t>(args.ny)
        * static_cast<std::size_t>(args.nz);

    std::vector<double> density(npoints);
    std::vector<double> log_density(npoints);
    std::vector<double> emissivity(npoints);
    std::vector<double> log_emissivity(npoints);
    std::vector<double> r_rg(npoints);
    std::vector<double> theta_values(npoints);
    std::vector<double> phi_values(npoints);
    std::vector<double> normalized_source(npoints);

    const double dx = 2.0 * args.box_rg / static_cast<double>(args.nx - 1);
    const double dy = 2.0 * args.box_rg / static_cast<double>(args.ny - 1);
    const double dz = 2.0 * args.box_rg / static_cast<double>(args.nz - 1);
    const double rho_log_floor = std::max(args.rho_floor, 1.0e-300);
    const double emissivity_floor = 1.0e-300;

    double max_source = 0.0;

#pragma omp parallel for reduction(max:max_source) schedule(static)
    for (int k = 0; k < args.nz; ++k) {
        for (int j = 0; j < args.ny; ++j) {
            for (int i = 0; i < args.nx; ++i) {
                const std::size_t idx =
                    static_cast<std::size_t>(i)
                    + static_cast<std::size_t>(args.nx)
                    * (static_cast<std::size_t>(j)
                    + static_cast<std::size_t>(args.ny) * static_cast<std::size_t>(k));

                const double x = -args.box_rg + dx * static_cast<double>(i);
                const double y = -args.box_rg + dy * static_cast<double>(j);
                const double z = -args.box_rg + dz * static_cast<double>(k);
                const double r = std::sqrt(x * x + y * y + z * z);
                const double theta =
                    r > 0.0
                    ? std::acos(std::clamp(z / r, -1.0, 1.0))
                    : 0.0;
                const double phi = r > 0.0 ? std::atan2(y, x) : 0.0;

                const double rho = finite_or_zero(torus.rho(r, theta));
                const double source_weight =
                    finite_or_zero(
                        radiative_transfer::uhe_source_spatial_weight(
                            r,
                            theta,
                            torus,
                            args.source
                        )
                    );
                const double j_uhe =
                    finite_or_zero(
                        radiative_transfer::emissivity_uhe(
                            r,
                            theta,
                            args.Enu_GeV,
                            torus,
                            args.source
                        )
                    );

                density[idx] = rho;
                log_density[idx] = std::log10(std::max(rho, rho_log_floor));
                emissivity[idx] = j_uhe;
                log_emissivity[idx] = std::log10(std::max(j_uhe, emissivity_floor));
                r_rg[idx] = r;
                theta_values[idx] = theta;
                phi_values[idx] = phi;
                normalized_source[idx] = source_weight;
                max_source = std::max(max_source, source_weight);
            }
        }
    }

    if (max_source > 0.0) {
        for (double& value : normalized_source) {
            value = finite_or_zero(value / max_source);
        }
    }

    std::ofstream out(args.output_path);
    if (!out) {
        std::cerr << "Could not open " << args.output_path << "\n";
        return 1;
    }

    out << "# vtk DataFile Version 3.0\n";
    out << "BH_Torus_RTX axisymmetric semi-analytic model sampled on Cartesian grid\n";
    out << "ASCII\n";
    out << "DATASET STRUCTURED_POINTS\n";
    out << "DIMENSIONS " << args.nx << " " << args.ny << " " << args.nz << "\n";
    out << std::scientific << std::setprecision(10);
    out << "ORIGIN " << -args.box_rg << " " << -args.box_rg << " " << -args.box_rg << "\n";
    out << "SPACING " << dx << " " << dy << " " << dz << "\n";
    out << "POINT_DATA " << npoints << "\n";

    write_scalar(out, "density_gcm3", density);
    write_scalar(out, "log10_density", log_density);
    write_scalar(out, "uhe_emissivity", emissivity);
    write_scalar(out, "log10_uhe_emissivity", log_emissivity);
    write_scalar(out, "r_rg", r_rg);
    write_scalar(out, "theta", theta_values);
    write_scalar(out, "phi", phi_values);
    write_scalar(out, "normalized_source", normalized_source);

    std::cout << "Saved: " << args.output_path << "\n";
    std::cout << "Grid points: " << npoints << "\n";
    return 0;
}
