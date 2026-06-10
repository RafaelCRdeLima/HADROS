#include "torus_profile.hpp"

#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>
#include <cstdlib>
#include <string>
#include <algorithm>

int main(int argc, char* argv[])
{
    double rho0 = 1.0e10;
    double r0_rg = 10.0;
    double sigma_rg = 5.0;
    double h_over_r = 0.01;
    std::string density_profile = "gaussian";
    double radial_power = 2.0;
    double funnel_depletion = 0.0;
    double funnel_theta_deg = 15.0;
    double envelope_rho0 = 0.0;
    double envelope_alpha = 2.5;
    double r_min_rg = 4.0;
    double r_max_rg = 60.0;
    double rho_floor = 1.0e-99;
    int Nr = 300;
    int Nth = 200;

    if (argc > 1) {
        rho0 = std::atof(argv[1]);
    }
    if (argc > 2) {
        r0_rg = std::atof(argv[2]);
    }
    if (argc > 3) {
        sigma_rg = std::atof(argv[3]);
    }
    if (argc > 4) {
        h_over_r = std::atof(argv[4]);
    }
    if (argc > 5) {
        density_profile = argv[5];
    }
    if (argc > 6) {
        radial_power = std::atof(argv[6]);
    }
    if (argc > 7) {
        funnel_depletion = std::atof(argv[7]);
    }
    if (argc > 8) {
        funnel_theta_deg = std::atof(argv[8]);
    }
    if (argc > 9) {
        envelope_rho0 = std::atof(argv[9]);
    }
    if (argc > 10) {
        envelope_alpha = std::atof(argv[10]);
    }
    if (argc > 11) {
        r_min_rg = std::atof(argv[11]);
    }
    if (argc > 12) {
        r_max_rg = std::atof(argv[12]);
    }
    if (argc > 13) {
        rho_floor = std::atof(argv[13]);
    }
    if (argc > 14) {
        Nr = std::atoi(argv[14]);
    }
    if (argc > 15) {
        Nth = std::atoi(argv[15]);
    }

    TorusProfile torus(
        rho0,
        r0_rg,
        sigma_rg,
        h_over_r,
        density_profile,
        radial_power,
        funnel_depletion,
        funnel_theta_deg * M_PI / 180.0,
        envelope_rho0,
        envelope_alpha,
        r_min_rg,
        r_max_rg,
        rho_floor
    );

    if (Nr <= 1) {
        Nr = 2;
    }
    if (Nth <= 1) {
        Nth = 2;
    }

    const double rmin = 1.0;
    const double rmax = std::max(30.0, r_max_rg);

    const double thmin = 0.0;
    const double thmax = M_PI;

    std::ofstream out("output/profiles/torus_synthetic_grid.dat");

    out << "# profile_type " << density_profile << "\n"
        << "# rho0 " << rho0 << "\n"
        << "# rho0_gcm3 " << rho0 << "\n"
        << "# r0 " << r0_rg << "\n"
        << "# r0_rg " << r0_rg << "\n"
        << "# sigma_r " << sigma_rg << "\n"
        << "# sigma_r_rg " << sigma_rg << "\n"
        << "# H_over_R " << h_over_r << "\n"
        << "# radial_power " << radial_power << "\n"
        << "# funnel_depletion " << funnel_depletion << "\n"
        << "# funnel_theta_deg " << funnel_theta_deg << "\n"
        << "# rho_floor " << rho_floor << "\n"
        << "# envelope_rho0_gcm3 " << envelope_rho0 << "\n"
        << "# envelope_rho0 " << envelope_rho0 << "\n"
        << "# envelope_alpha " << envelope_alpha << "\n"
        << "# R_min " << r_min_rg << "\n"
        << "# R_min_rg " << r_min_rg << "\n"
        << "# R_max " << r_max_rg << "\n"
        << "# R_max_rg " << r_max_rg << "\n"
        << "# rho_floor_gcm3 " << rho_floor << "\n";
    out << "# r_rg theta rho_gcm3 T_MeV Ye\n";

    for (int i = 0; i < Nr; ++i) {
        double r = rmin + (rmax - rmin) * i / (Nr - 1);

        for (int j = 0; j < Nth; ++j) {
            double th = thmin + (thmax - thmin) * j / (Nth - 1);

            out << std::scientific << std::setprecision(8)
                << r << " "
                << th << " "
                << torus.rho(r, th) << " "
                << torus.temperature_MeV(r, th) << " "
                << torus.Ye(r, th) << "\n";
        }

        out << "\n";
    }

    std::cout << "Saved: output/profiles/torus_synthetic_grid.dat\n";

    return 0;
}
