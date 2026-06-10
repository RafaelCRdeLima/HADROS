#ifndef RAY_HPP
#define RAY_HPP

#include <vector>

struct PathPoint {
    double r_rg;
    double theta;
    double x_rg;
    double y_rg;
    double z_rg;
    double dl_rg;
    double redshift_factor;
};

struct RayPath {
    int pixel_i = 0;
    int pixel_j = 0;

    double alpha_rg = 0.0;
    double beta_rg  = 0.0;
    double impact_parameter_rg = 0.0;

    double a_bh = 0.95;

    bool captured = false;

    std::vector<PathPoint> points;
};

#endif