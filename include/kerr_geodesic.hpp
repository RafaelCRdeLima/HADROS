#ifndef KERR_GEODESIC_HPP
#define KERR_GEODESIC_HPP

#include "kerr_metric.hpp"
#include "geodesic_state.hpp"

class KerrGeodesic {
public:
    explicit KerrGeodesic(
        KerrMetric metric,
        double h = 0.02,
        double tolerance = 1.0e-6
    );

    double hamiltonian(
        const GeodesicState& y
    ) const;

    GeodesicState rhs(
        const GeodesicState& y
    ) const;

    void step_rk4(
    GeodesicState& y
    ) const;

    void step_adaptive(
        GeodesicState& y
    ) const;

private:
    KerrMetric metric_;
    double h_;
    double tolerance_;

    double dg_inv_dr(
        int mu,
        int nu,
        double r,
        double th
    ) const;

    double dg_inv_dtheta(
        int mu,
        int nu,
        double r,
        double th
    ) const;
};

#endif