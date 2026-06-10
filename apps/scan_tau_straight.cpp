#include "sigma_table.hpp"
#include "torus_profile.hpp"
#include "optical_depth.hpp"

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>

int main(int argc, char* argv[])
{
    double Enu_GeV = 1.0e9;

    if (argc > 1) {
        Enu_GeV = std::atof(argv[1]);
    }

    // ============================================
    // Black hole mass
    // ============================================

    const double M_bh_msun = 3.0;

    // ============================================
    // Straight ray integration parameters
    // ============================================

    const double xmax_rg = 40.0;
    const double dx_rg   = 0.02;

    // ============================================
    // Impact parameter scan
    // ============================================

    const double bmin = 0.0;
    const double bmax = 25.0;

    const int Nb = 300;

    // ============================================
    // Objects
    // ============================================

    SigmaTable sigma(
        "data/sigma/sigma_nuN_CC_GBW.dat"
    );

    TorusProfile torus;

    // ============================================
    // Output
    // ============================================

    std::ofstream out(
        "output/tau/tau_straight.dat"
    );

    out << "# b_rg tau P_surv\n";

    // ============================================
    // Scan
    // ============================================

    for (int i = 0; i < Nb; ++i) {

        double b =
            bmin + (bmax - bmin) * i / (Nb - 1);

        double tau =
            optical_depth::tau_straight_ray(
                b,
                Enu_GeV,
                xmax_rg,
                dx_rg,
                M_bh_msun,
                torus,
                sigma
            );

        double P =
            optical_depth::survival_probability(tau);

        out << std::scientific
            << std::setprecision(8)
            << b   << " "
            << tau << " "
            << P   << "\n";
    }

    std::cout << "Saved: output/tau/tau_straight.dat\n";

    return 0;
}