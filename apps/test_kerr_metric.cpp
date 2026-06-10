#include "kerr_metric.hpp"

#include <iostream>
#include <iomanip>
#include <cmath>

int main()
{
    const double a = 0.9;

    KerrMetric kerr(a);

    const double r  = 5.0;
    const double th = M_PI / 2.0;

    double g[4][4];
    double ginv[4][4];

    kerr.metric(r, th, g);
    kerr.inverse_metric(r, th, ginv);

    std::cout << std::scientific
              << std::setprecision(8);

    std::cout << "a = " << a << "\n";

    std::cout << "Horizon radius = "
              << kerr.horizon_radius()
              << "\n\n";

    std::cout << "Metric components at r="
              << r
              << ", theta=pi/2\n\n";

    std::cout << "g_tt   = " << g[0][0] << "\n";
    std::cout << "g_tphi = " << g[0][3] << "\n";
    std::cout << "g_rr   = " << g[1][1] << "\n";
    std::cout << "g_thth = " << g[2][2] << "\n";
    std::cout << "g_phph = " << g[3][3] << "\n\n";

    std::cout << "Inverse metric:\n\n";

    std::cout << "gtt    = " << ginv[0][0] << "\n";
    std::cout << "gtphi  = " << ginv[0][3] << "\n";
    std::cout << "grr    = " << ginv[1][1] << "\n";
    std::cout << "gthth  = " << ginv[2][2] << "\n";
    std::cout << "gphph  = " << ginv[3][3] << "\n\n";

    std::cout << "Lapse = "
              << kerr.lapse(r, th)
              << "\n";

    std::cout << "Frame dragging omega = "
              << kerr.omega_frame_drag(r, th)
              << "\n";

    return 0;
}