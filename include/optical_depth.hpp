#ifndef OPTICAL_DEPTH_HPP
#define OPTICAL_DEPTH_HPP

#include "sigma_table.hpp"
#include "torus_profile.hpp"
#include "ray.hpp"


namespace optical_depth {

    double tau_straight_ray(
        double impact_parameter_rg,
        double Enu_GeV,
        double xmax_rg,
        double dx_rg,
        double M_bh_msun,
        const TorusProfile& torus,
        const SigmaTable& sigma
    );

    double tau_along_ray(
        const RayPath& ray,
        double Enu_inf_GeV,
        const TorusProfile& torus,
        const SigmaTable& sigma
    );

    double survival_probability(double tau);

}

#endif