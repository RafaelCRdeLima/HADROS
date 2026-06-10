#include "sigma_table.hpp"
#include <iostream>

int main()
{
    SigmaTable sig("data/sigma/sigma_nuN_CC_GBW.dat");

    double E = 1.0e9;

    std::cout << "N = " << sig.size() << "\n";
    std::cout << "E range = " << sig.Emin() << " -- " << sig.Emax() << " GeV\n";
    std::cout << "sigma(" << E << " GeV) = "
              << sig.sigma_cm2(E) << " cm^2\n";

    return 0;
}