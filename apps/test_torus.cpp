#include "torus_profile.hpp"
#include <iostream>
#include <cmath>

int main()
{
    TorusProfile torus;

    double r = 10.0;
    double th = M_PI / 2.0;

    std::cout << "rho = " << torus.rho(r, th) << " g/cm^3\n";
    std::cout << "T   = " << torus.temperature_MeV(r, th) << " MeV\n";
    std::cout << "Ye  = " << torus.Ye(r, th) << "\n";

    return 0;
}