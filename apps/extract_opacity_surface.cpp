#include "constants.hpp"
#include "sigma_table.hpp"
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

namespace {

struct Args {
    double Enu_GeV = 1.0e5;
    double M_bh_msun = 3.0;
    std::string sigma_path = "data/sigma/sigma_nuN_CC_GBW.dat";
    std::string profile = "gaussian";
    double rho0 = 1.0e-2;
    double r0 = 10.0;
    double sigma_r = 5.0;
    double H_over_R = 0.25;
    double radial_power = 2.0;
    double funnel_depletion = 0.0;
    double funnel_theta_deg = 15.0;
    double envelope_rho0 = 0.0;
    double envelope_alpha = 2.5;
    double R_min = 4.0;
    double R_max = 60.0;
    double rho_floor = 1.0e-99;
    double tau_surface = 1.0;
    int ntheta = 181;
    int nr = 1200;
    std::string output_path = "output/opacity_surfaces/tau_surface_tau1.dat";
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

double classify_tau(double tau)
{
    if (tau < 0.1) {
        return 0.0;
    }
    if (tau < 1.0) {
        return 1.0;
    }
    if (tau < 3.0) {
        return 2.0;
    }
    return 3.0;
}

}

int main(int argc, char* argv[])
{
    Args args;
    args.Enu_GeV = arg_double(argv, argc, 1, args.Enu_GeV);
    args.M_bh_msun = arg_double(argv, argc, 2, args.M_bh_msun);
    args.sigma_path = arg_string(argv, argc, 3, args.sigma_path);
    args.profile = arg_string(argv, argc, 4, args.profile);
    args.rho0 = arg_double(argv, argc, 5, args.rho0);
    args.r0 = arg_double(argv, argc, 6, args.r0);
    args.sigma_r = arg_double(argv, argc, 7, args.sigma_r);
    args.H_over_R = arg_double(argv, argc, 8, args.H_over_R);
    args.radial_power = arg_double(argv, argc, 9, args.radial_power);
    args.funnel_depletion = arg_double(argv, argc, 10, args.funnel_depletion);
    args.funnel_theta_deg = arg_double(argv, argc, 11, args.funnel_theta_deg);
    args.envelope_rho0 = arg_double(argv, argc, 12, args.envelope_rho0);
    args.envelope_alpha = arg_double(argv, argc, 13, args.envelope_alpha);
    args.R_min = arg_double(argv, argc, 14, args.R_min);
    args.R_max = arg_double(argv, argc, 15, args.R_max);
    args.rho_floor = arg_double(argv, argc, 16, args.rho_floor);
    args.tau_surface = arg_double(argv, argc, 17, args.tau_surface);
    args.ntheta = arg_int(argv, argc, 18, args.ntheta);
    args.nr = arg_int(argv, argc, 19, args.nr);
    args.output_path = arg_string(argv, argc, 20, args.output_path);

    if (args.ntheta < 2 || args.nr < 4 || args.R_max <= args.R_min) {
        std::cerr << "Invalid opacity-surface grid parameters.\n";
        return 1;
    }

    SigmaTable sigma(args.sigma_path);
    if (args.Enu_GeV < sigma.Emin() || args.Enu_GeV > sigma.Emax()) {
        std::cerr << "Energy outside sigma table range.\n";
        return 1;
    }

    TorusProfile torus(
        args.rho0,
        args.r0,
        args.sigma_r,
        args.H_over_R,
        args.profile,
        args.radial_power,
        args.funnel_depletion,
        args.funnel_theta_deg * constants::pi / 180.0,
        args.envelope_rho0,
        args.envelope_alpha,
        args.R_min,
        args.R_max,
        args.rho_floor
    );

    std::ofstream out(args.output_path);
    if (!out) {
        std::cerr << "Could not open " << args.output_path << "\n";
        return 1;
    }

    const double sigma_cm2 = sigma.sigma_cm2(args.Enu_GeV);
    const double rg_cm = constants::rg_cm(args.M_bh_msun);
    const double dr = (args.R_max - args.R_min) / static_cast<double>(args.nr - 1);

    out << "# opacity_surface axisymmetric_radial\n"
        << "# Enu_GeV " << args.Enu_GeV << "\n"
        << "# tau_surface " << args.tau_surface << "\n"
        << "# profile_type " << args.profile << "\n"
        << "# rho0 " << args.rho0 << "\n"
        << "# r0 " << args.r0 << "\n"
        << "# sigma_r " << args.sigma_r << "\n"
        << "# H_over_R " << args.H_over_R << "\n"
        << "# radial_power " << args.radial_power << "\n"
        << "# funnel_depletion " << args.funnel_depletion << "\n"
        << "# funnel_theta_deg " << args.funnel_theta_deg << "\n"
        << "# envelope_rho0 " << args.envelope_rho0 << "\n"
        << "# envelope_alpha " << args.envelope_alpha << "\n"
        << "# R_min " << args.R_min << "\n"
        << "# R_max " << args.R_max << "\n"
        << "# rho_floor " << args.rho_floor << "\n"
        << "# sigma_path " << args.sigma_path << "\n"
        << "# sigma_cm2 " << sigma_cm2 << "\n"
        << "# M_bh_msun " << args.M_bh_msun << "\n"
        << "# angular_sampling_uniform_theta " << args.ntheta << "\n"
        << "# radial_sampling_uniform_rg " << args.nr << "\n"
        << "# interpolation linear_in_tau_between_radial_samples\n"
        << "# no_crossing_r_tau_rg -1\n"
        << "# future_interface r_tau(theta,phi)\n"
        << "# columns theta_rad theta_deg r_tau_rg tau_total classification crossing_found\n";

    std::vector<double> rho(static_cast<std::size_t>(args.nr), 0.0);
    std::vector<double> tau(static_cast<std::size_t>(args.nr), 0.0);

    for (int it = 0; it < args.ntheta; ++it) {
        const double theta =
            constants::pi * static_cast<double>(it) / static_cast<double>(args.ntheta - 1);

        for (int ir = 0; ir < args.nr; ++ir) {
            const double r = args.R_min + dr * static_cast<double>(ir);
            rho[static_cast<std::size_t>(ir)] = torus.rho(r, theta);
        }

        tau[static_cast<std::size_t>(args.nr - 1)] = 0.0;
        for (int ir = args.nr - 2; ir >= 0; --ir) {
            const double rho_mid =
                0.5 * (rho[static_cast<std::size_t>(ir)] + rho[static_cast<std::size_t>(ir + 1)]);
            const double nb = rho_mid / constants::m_u_g;
            const double dtau = nb * sigma_cm2 * dr * rg_cm;
            tau[static_cast<std::size_t>(ir)] =
                tau[static_cast<std::size_t>(ir + 1)] + dtau;
        }

        const double tau_total = tau.front();
        double r_tau = -1.0;
        int crossing_found = 0;

        if (tau_total >= args.tau_surface) {
            for (int ir = 0; ir < args.nr - 1; ++ir) {
                const double t0 = tau[static_cast<std::size_t>(ir)];
                const double t1 = tau[static_cast<std::size_t>(ir + 1)];
                if (t0 >= args.tau_surface && t1 <= args.tau_surface) {
                    const double r0 = args.R_min + dr * static_cast<double>(ir);
                    const double frac =
                        (args.tau_surface - t0) / std::max(t1 - t0, -1.0e-300);
                    r_tau = r0 + std::clamp(frac, 0.0, 1.0) * dr;
                    crossing_found = 1;
                    break;
                }
            }
        }

        out << std::scientific << std::setprecision(10)
            << theta << " "
            << theta * 180.0 / constants::pi << " "
            << r_tau << " "
            << tau_total << " "
            << classify_tau(tau_total) << " "
            << crossing_found << "\n";
    }

    std::cout << "Saved: " << args.output_path << "\n";
    return 0;
}
