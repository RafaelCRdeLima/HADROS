#include <cuda_runtime.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

constexpr double PI = 3.141592653589793238462643383279502884;
constexpr double C_CGS = 2.99792458e10;
constexpr double G_CGS = 6.67430e-8;
constexpr double MSUN_G = 1.98847e33;
constexpr double M_U_G = 1.66053906660e-24;

#define CUDA_CHECK(call)                                                       \
    do {                                                                       \
        cudaError_t err__ = (call);                                            \
        if (err__ != cudaSuccess) {                                            \
            std::ostringstream oss__;                                          \
            oss__ << "CUDA error at " << __FILE__ << ":" << __LINE__          \
                  << ": " << cudaGetErrorString(err__);                       \
            throw std::runtime_error(oss__.str());                             \
        }                                                                      \
    } while (0)

struct HostSigmaTable {
    std::vector<double> energy_GeV;
    std::vector<double> sigma_cm2;
};

struct Params {
    double Enu_obs_GeV;
    double M_bh_msun;
    double a_spin;
    double r_obs_rg;
    double theta_obs_rad;
    double fov_rad;
    double r_max_rg;
    double h;
    double tolerance;
    int nx;
    int ny;
    int max_steps;
    double torus_rho0;
    double torus_r0_rg;
    double torus_sigma_rg;
    double torus_h_over_r;
    int density_profile_kind;
    double torus_radial_power;
    double funnel_depletion;
    double funnel_theta_rad;
    double envelope_rho0;
    double envelope_alpha;
    double torus_r_min_rg;
    double torus_r_max_rg;
    double rho_floor;
    double source_r_center_rg;
    double source_sigma_r_rg;
    double source_theta_width_rad;
    double source_powerlaw;
    double source_emax_GeV;
    double source_norm;
    int source_model_kind;
    double source_funnel_theta_rad;
    double source_density_power_q;
    double source_radial_power_s;
    double source_rho_ref_gcm3;
    double source_cutoff_min;
    double source_cutoff_max;
    double source_gradient_dr_rg;
    double source_gradient_dtheta_rad;
    double mev_Enu_obs_MeV;
    double mev_norm;
    double mev_sigma_abs0_cm2;
    double mev_sigma_scat0_cm2;
    double mev_neutrinosphere_tau;
    int sigma_n;
};

enum SourceModelKind {
    SOURCE_INNER_RING = 0,
    SOURCE_FUNNEL_WALL = 1,
    SOURCE_JET_BASE = 2,
    SOURCE_SHOCK_LAYER = 3,
    SOURCE_DENSITY_WEIGHTED = 4
};

struct State {
    double t;
    double r;
    double theta;
    double phi;
    double pt;
    double pr;
    double ptheta;
    double pphi;
};

int parse_density_profile_kind(const std::string& profile_name)
{
    if (profile_name == "gaussian" ||
        profile_name == "torus" ||
        profile_name == "gaussian_torus") {
        return 0;
    }

    if (profile_name == "powerlaw" ||
        profile_name == "powerlaw_disk") {
        return 1;
    }

    if (profile_name == "funnel" ||
        profile_name == "gaussian_funnel") {
        return 2;
    }

    if (profile_name == "powerlaw_funnel") {
        return 3;
    }

    if (profile_name == "gaussian_envelope") {
        return 4;
    }

    if (profile_name == "powerlaw_envelope") {
        return 5;
    }

    if (profile_name == "powerlaw_funnel_envelope") {
        return 6;
    }

    throw std::runtime_error(
        "Unknown density profile '" + profile_name +
        "'. Use gaussian, powerlaw, gaussian_funnel, powerlaw_funnel, "
        "gaussian_envelope, powerlaw_envelope, or powerlaw_funnel_envelope."
    );
}

int parse_source_model_kind(const std::string& model_name)
{
    if (model_name == "inner_ring" ||
        model_name == "INNER_RING" ||
        model_name == "ring") {
        return SOURCE_INNER_RING;
    }

    if (model_name == "funnel_wall" ||
        model_name == "FUNNEL_WALL") {
        return SOURCE_FUNNEL_WALL;
    }

    if (model_name == "jet_base" ||
        model_name == "JET_BASE") {
        return SOURCE_JET_BASE;
    }

    if (model_name == "shock_layer" ||
        model_name == "SHOCK_LAYER") {
        return SOURCE_SHOCK_LAYER;
    }

    if (model_name == "density_weighted" ||
        model_name == "DENSITY_WEIGHTED") {
        return SOURCE_DENSITY_WEIGHTED;
    }

    throw std::runtime_error(
        "Unknown UHE source model '" + model_name +
        "'. Use inner_ring, funnel_wall, jet_base, shock_layer, or density_weighted."
    );
}

const char* source_model_name_from_kind(int kind)
{
    switch (kind) {
        case SOURCE_FUNNEL_WALL:
            return "funnel_wall";
        case SOURCE_JET_BASE:
            return "jet_base";
        case SOURCE_SHOCK_LAYER:
            return "shock_layer";
        case SOURCE_DENSITY_WEIGHTED:
            return "density_weighted";
        case SOURCE_INNER_RING:
        default:
            return "inner_ring";
    }
}

struct RayAccum {
    int pixel_i;
    int pixel_j;
    int captured;
    int finished;
    int steps_done;
    double alpha_rg;
    double beta_rg;
    double prev_r;
    double prev_theta;
    double prev_phi;
    double energy_obs;
    double tau_uhe;
    double I_uhe;
    double tau_mev;
    double I_mev;
    double r_neutrinosphere_rg;
    double leakage_sum;
    double leakage_weight;
    State y;
};

struct OutputPixel {
    int pixel_i;
    int pixel_j;
    int captured;
    double alpha_rg;
    double beta_rg;
    double tau_uhe;
    double P_surv_uhe;
    double I_uhe;
    double tau_mev;
    double P_surv_mev;
    double I_mev;
    double r_neutrinosphere_rg;
    double leakage_mev;
};

struct CachePoint {
    double r_rg;
    double theta;
    double x_rg;
    double y_rg;
    double z_rg;
    double dl_rg;
    double redshift_factor;
};

struct BinaryRayHeader {
    std::int32_t ray_id;
    std::int32_t pixel_i;
    std::int32_t pixel_j;
    std::int32_t captured;
    std::int32_t npoints;
    double alpha_rg;
    double beta_rg;
};

struct CacheRayDevice {
    int pixel_i;
    int pixel_j;
    int captured;
    int npoints;
    int offset;
    double alpha_rg;
    double beta_rg;
};

void write_i32(std::ofstream& out, std::int32_t value)
{
    out.write(reinterpret_cast<const char*>(&value), sizeof(value));
}

void write_f64(std::ofstream& out, double value)
{
    out.write(reinterpret_cast<const char*>(&value), sizeof(value));
}

void write_image_binary(
    const std::string& filename,
    const std::vector<OutputPixel>& image,
    const Params& params,
    double cam_theta_deg
)
{
    std::ofstream out(filename, std::ios::binary);
    if (!out) {
        throw std::runtime_error("Could not open binary output file: " + filename);
    }

    const std::int32_t magic = 0x4B494D47; // "KIMG"
    const std::int32_t version = 1;
    const std::int32_t backend_cuda = 2;
    const std::int32_t ncols = 13;
    const std::int32_t reserved = 0;

    write_i32(out, magic);
    write_i32(out, version);
    write_i32(out, backend_cuda);
    write_i32(out, static_cast<std::int32_t>(params.nx));
    write_i32(out, static_cast<std::int32_t>(params.ny));
    write_i32(out, ncols);
    write_i32(out, reserved);

    write_f64(out, params.Enu_obs_GeV);
    write_f64(out, params.mev_Enu_obs_MeV);
    write_f64(out, params.mev_norm);
    write_f64(out, cam_theta_deg);

    for (const OutputPixel& p : image) {
        const double row[ncols] = {
            static_cast<double>(p.pixel_i),
            static_cast<double>(p.pixel_j),
            p.alpha_rg,
            p.beta_rg,
            p.tau_uhe,
            p.P_surv_uhe,
            p.I_uhe,
            static_cast<double>(p.captured),
            p.tau_mev,
            p.P_surv_mev,
            p.I_mev,
            p.r_neutrinosphere_rg,
            p.leakage_mev
        };

        out.write(
            reinterpret_cast<const char*>(row),
            sizeof(row)
        );
    }
}

HostSigmaTable load_sigma_table(const std::string& filename)
{
    std::ifstream in(filename);
    if (!in) {
        throw std::runtime_error("Could not open sigma table: " + filename);
    }

    HostSigmaTable table;
    std::string line;

    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#') {
            continue;
        }

        std::istringstream iss(line);
        double E = 0.0;
        double sigma_GeV2 = 0.0;
        double sigma_cm2 = 0.0;

        if (!(iss >> E >> sigma_GeV2 >> sigma_cm2)) {
            continue;
        }

        if (E > 0.0 && sigma_GeV2 > 0.0 && sigma_cm2 > 0.0) {
            table.energy_GeV.push_back(E);
            table.sigma_cm2.push_back(sigma_cm2);
        }
    }

    if (table.energy_GeV.size() < 2) {
        throw std::runtime_error("Sigma table must contain at least two valid data lines.");
    }

    for (std::size_t i = 1; i < table.energy_GeV.size(); ++i) {
        if (table.energy_GeV[i] <= table.energy_GeV[i - 1]) {
            throw std::runtime_error("Energy grid in sigma table must be strictly increasing.");
        }
    }

    return table;
}

std::string make_tag(const char* fmt, double value)
{
    char buffer[64];
    std::snprintf(buffer, sizeof(buffer), fmt, value);
    std::string tag(buffer);
    std::replace(tag.begin(), tag.end(), '+', '_');
    std::replace(tag.begin(), tag.end(), '-', 'm');
    std::replace(tag.begin(), tag.end(), '.', 'p');
    return tag;
}

std::string make_compact_scientific_tag(double value)
{
    char buffer[64];
    std::snprintf(buffer, sizeof(buffer), "%.0e", value);

    std::string tag(buffer);
    const std::size_t epos = tag.find('e');

    if (epos == std::string::npos) {
        return tag;
    }

    const std::string mantissa = tag.substr(0, epos);
    const int exponent = std::stoi(tag.substr(epos + 1));

    return mantissa + "e" + std::to_string(exponent);
}

__host__ __device__ double rg_cm(double M_msun)
{
    return G_CGS * M_msun * MSUN_G / (C_CGS * C_CGS);
}

__device__ double clamp_device(double x, double lo, double hi)
{
    return fmin(fmax(x, lo), hi);
}

__device__ double kerr_sigma(double r, double th, double a)
{
    const double c = cos(th);
    return r * r + a * a * c * c;
}

__device__ double kerr_delta(double r, double a)
{
    return r * r - 2.0 * r + a * a;
}

__device__ double kerr_big_a(double r, double th, double a)
{
    const double s = sin(th);
    const double s2 = s * s;
    const double rr_aa = r * r + a * a;
    return rr_aa * rr_aa - a * a * kerr_delta(r, a) * s2;
}

__device__ double kerr_horizon(double a)
{
    return 1.0 + sqrt(fmax(1.0 - a * a, 0.0));
}

__device__ void kerr_metric(double r, double th, double a, double g[4][4])
{
    const double sig = kerr_sigma(r, th, a);
    const double del = kerr_delta(r, a);
    const double s = sin(th);
    const double s2 = s * s;

    for (int mu = 0; mu < 4; ++mu) {
        for (int nu = 0; nu < 4; ++nu) {
            g[mu][nu] = 0.0;
        }
    }

    g[0][0] = -(1.0 - 2.0 * r / sig);
    g[0][3] = -2.0 * a * r * s2 / sig;
    g[3][0] = g[0][3];
    g[1][1] = sig / del;
    g[2][2] = sig;
    g[3][3] = (r * r + a * a + 2.0 * a * a * r * s2 / sig) * s2;
}

__device__ void kerr_inverse_metric(double r, double th, double a, double ginv[4][4])
{
    const double sig = kerr_sigma(r, th, a);
    const double del = kerr_delta(r, a);
    const double s = sin(th);
    const double s2 = fmax(s * s, 1.0e-300);
    const double bigA = kerr_big_a(r, th, a);

    for (int mu = 0; mu < 4; ++mu) {
        for (int nu = 0; nu < 4; ++nu) {
            ginv[mu][nu] = 0.0;
        }
    }

    ginv[0][0] = -bigA / (sig * del);
    ginv[0][3] = -2.0 * a * r / (sig * del);
    ginv[3][0] = ginv[0][3];
    ginv[1][1] = del / sig;
    ginv[2][2] = 1.0 / sig;
    ginv[3][3] = (del - a * a * s2) / (sig * del * s2);
}

__device__ double kerr_lapse(double r, double th, double a)
{
    const double sig = kerr_sigma(r, th, a);
    const double del = kerr_delta(r, a);
    const double bigA = kerr_big_a(r, th, a);
    return sqrt(fmax(sig * del / bigA, 1.0e-300));
}

__device__ double kerr_omega(double r, double th, double a)
{
    const double s = sin(th);
    const double s2 = s * s;
    return 2.0 * a * r / kerr_big_a(r, th, a);
}

__device__ double zamo_energy(const Params& params, double r, double theta, double p_t, double p_phi)
{
    const double alpha = kerr_lapse(r, theta, params.a_spin);
    const double omega = kerr_omega(r, theta, params.a_spin);
    return -(p_t + omega * p_phi) / alpha;
}

__device__ double wrapped_delta_phi(double phi_new, double phi_old)
{
    double dphi = phi_new - phi_old;

    while (dphi > PI) {
        dphi -= 2.0 * PI;
    }

    while (dphi < -PI) {
        dphi += 2.0 * PI;
    }

    return dphi;
}

__device__ double zamo_spatial_interval_rg(const Params& params, const State& current, const RayAccum& ray)
{
    const double r_mid = 0.5 * (current.r + ray.prev_r);
    const double theta_mid = 0.5 * (current.theta + ray.prev_theta);
    double g[4][4];
    kerr_metric(r_mid, theta_mid, params.a_spin, g);

    const double dr = current.r - ray.prev_r;
    const double dtheta = current.theta - ray.prev_theta;
    const double dphi = wrapped_delta_phi(current.phi, ray.prev_phi);
    const double dl2 =
        g[1][1] * dr * dr
        + g[2][2] * dtheta * dtheta
        + g[3][3] * dphi * dphi;

    return sqrt(fmax(dl2, 0.0));
}

__device__ double dg_inv_dr(int mu, int nu, double r, double th, double a)
{
    const double eps = 1.0e-5 * fmax(1.0, fabs(r));
    double gp[4][4];
    double gm[4][4];
    kerr_inverse_metric(r + eps, th, a, gp);
    kerr_inverse_metric(r - eps, th, a, gm);
    return (gp[mu][nu] - gm[mu][nu]) / (2.0 * eps);
}

__device__ double dg_inv_dtheta(int mu, int nu, double r, double th, double a)
{
    const double eps = 1.0e-5;
    double gp[4][4];
    double gm[4][4];
    kerr_inverse_metric(r, th + eps, a, gp);
    kerr_inverse_metric(r, th - eps, a, gm);
    return (gp[mu][nu] - gm[mu][nu]) / (2.0 * eps);
}

__device__ State geodesic_rhs(const Params& params, const State& y)
{
    double ginv[4][4];
    kerr_inverse_metric(y.r, y.theta, params.a_spin, ginv);

    const double p[4] = {y.pt, y.pr, y.ptheta, y.pphi};
    State dydl{};

    for (int nu = 0; nu < 4; ++nu) {
        dydl.t += ginv[0][nu] * p[nu];
        dydl.r += ginv[1][nu] * p[nu];
        dydl.theta += ginv[2][nu] * p[nu];
        dydl.phi += ginv[3][nu] * p[nu];
    }

    dydl.pt = 0.0;
    dydl.pphi = 0.0;

    for (int mu = 0; mu < 4; ++mu) {
        for (int nu = 0; nu < 4; ++nu) {
            dydl.pr -= 0.5 * dg_inv_dr(mu, nu, y.r, y.theta, params.a_spin) * p[mu] * p[nu];
            dydl.ptheta -= 0.5 * dg_inv_dtheta(mu, nu, y.r, y.theta, params.a_spin) * p[mu] * p[nu];
        }
    }

    return dydl;
}

__device__ State add_scaled(const State& y, const State& k, double h)
{
    State out;
    out.t = y.t + h * k.t;
    out.r = y.r + h * k.r;
    out.theta = y.theta + h * k.theta;
    out.phi = y.phi + h * k.phi;
    out.pt = y.pt + h * k.pt;
    out.pr = y.pr + h * k.pr;
    out.ptheta = y.ptheta + h * k.ptheta;
    out.pphi = y.pphi + h * k.pphi;
    return out;
}

__device__ State add_scaled5(
    const State& y,
    const State& k1,
    double a1,
    const State& k2,
    double a2,
    const State& k3,
    double a3,
    const State& k4,
    double a4,
    const State& k5,
    double a5
)
{
    State out;
    out.t = y.t + a1 * k1.t + a2 * k2.t + a3 * k3.t + a4 * k4.t + a5 * k5.t;
    out.r = y.r + a1 * k1.r + a2 * k2.r + a3 * k3.r + a4 * k4.r + a5 * k5.r;
    out.theta = y.theta + a1 * k1.theta + a2 * k2.theta + a3 * k3.theta + a4 * k4.theta + a5 * k5.theta;
    out.phi = y.phi + a1 * k1.phi + a2 * k2.phi + a3 * k3.phi + a4 * k4.phi + a5 * k5.phi;
    out.pt = y.pt + a1 * k1.pt + a2 * k2.pt + a3 * k3.pt + a4 * k4.pt + a5 * k5.pt;
    out.pr = y.pr + a1 * k1.pr + a2 * k2.pr + a3 * k3.pr + a4 * k4.pr + a5 * k5.pr;
    out.ptheta = y.ptheta + a1 * k1.ptheta + a2 * k2.ptheta + a3 * k3.ptheta + a4 * k4.ptheta + a5 * k5.ptheta;
    out.pphi = y.pphi + a1 * k1.pphi + a2 * k2.pphi + a3 * k3.pphi + a4 * k4.pphi + a5 * k5.pphi;
    return out;
}

__device__ double error_norm(const State& a, const State& b)
{
    double err = 0.0;
    err = fmax(err, fabs(a.r - b.r));
    err = fmax(err, fabs(a.theta - b.theta));
    err = fmax(err, fabs(a.phi - b.phi));
    err = fmax(err, fabs(a.pr - b.pr));
    err = fmax(err, fabs(a.ptheta - b.ptheta));
    err = fmax(err, fabs(a.pphi - b.pphi));
    return err;
}

__device__ void rk4_step(const Params& params, State& y)
{
    const double h = params.h;
    const State k1 = geodesic_rhs(params, y);
    const State k2 = geodesic_rhs(params, add_scaled(y, k1, 0.5 * h));
    const State k3 = geodesic_rhs(params, add_scaled(y, k2, 0.5 * h));
    const State k4 = geodesic_rhs(params, add_scaled(y, k3, h));

    y.t += h * (k1.t + 2.0 * k2.t + 2.0 * k3.t + k4.t) / 6.0;
    y.r += h * (k1.r + 2.0 * k2.r + 2.0 * k3.r + k4.r) / 6.0;
    y.theta += h * (k1.theta + 2.0 * k2.theta + 2.0 * k3.theta + k4.theta) / 6.0;
    y.phi += h * (k1.phi + 2.0 * k2.phi + 2.0 * k3.phi + k4.phi) / 6.0;
    y.pt += h * (k1.pt + 2.0 * k2.pt + 2.0 * k3.pt + k4.pt) / 6.0;
    y.pr += h * (k1.pr + 2.0 * k2.pr + 2.0 * k3.pr + k4.pr) / 6.0;
    y.ptheta += h * (k1.ptheta + 2.0 * k2.ptheta + 2.0 * k3.ptheta + k4.ptheta) / 6.0;
    y.pphi += h * (k1.pphi + 2.0 * k2.pphi + 2.0 * k3.pphi + k4.pphi) / 6.0;
}

__device__ void rkf45_step_adaptive(const Params& params, State& y)
{
    double h = params.h;
    const double h_min = params.h * 1.0e-5;

    for (int attempt = 0; attempt < 50; ++attempt) {
        const State zero{};
        const State k1 = geodesic_rhs(params, y);
        const State k2 = geodesic_rhs(params, add_scaled5(
            y, k1, h * 1.0 / 4.0, zero, 0.0, zero, 0.0, zero, 0.0, zero, 0.0));
        const State k3 = geodesic_rhs(params, add_scaled5(
            y, k1, h * 3.0 / 32.0, k2, h * 9.0 / 32.0, zero, 0.0, zero, 0.0, zero, 0.0));
        const State k4 = geodesic_rhs(params, add_scaled5(
            y, k1, h * 1932.0 / 2197.0, k2, h * -7200.0 / 2197.0,
            k3, h * 7296.0 / 2197.0, zero, 0.0, zero, 0.0));
        const State k5 = geodesic_rhs(params, add_scaled5(
            y, k1, h * 439.0 / 216.0, k2, h * -8.0,
            k3, h * 3680.0 / 513.0, k4, h * -845.0 / 4104.0, zero, 0.0));
        const State k6 = geodesic_rhs(params, add_scaled5(
            y, k1, h * -8.0 / 27.0, k2, h * 2.0,
            k3, h * -3544.0 / 2565.0, k4, h * 1859.0 / 4104.0,
            k5, h * -11.0 / 40.0));

        State y4;
        y4.t = y.t + h * (25.0 / 216.0 * k1.t + 1408.0 / 2565.0 * k3.t + 2197.0 / 4104.0 * k4.t - 1.0 / 5.0 * k5.t);
        y4.r = y.r + h * (25.0 / 216.0 * k1.r + 1408.0 / 2565.0 * k3.r + 2197.0 / 4104.0 * k4.r - 1.0 / 5.0 * k5.r);
        y4.theta = y.theta + h * (25.0 / 216.0 * k1.theta + 1408.0 / 2565.0 * k3.theta + 2197.0 / 4104.0 * k4.theta - 1.0 / 5.0 * k5.theta);
        y4.phi = y.phi + h * (25.0 / 216.0 * k1.phi + 1408.0 / 2565.0 * k3.phi + 2197.0 / 4104.0 * k4.phi - 1.0 / 5.0 * k5.phi);
        y4.pt = y.pt + h * (25.0 / 216.0 * k1.pt + 1408.0 / 2565.0 * k3.pt + 2197.0 / 4104.0 * k4.pt - 1.0 / 5.0 * k5.pt);
        y4.pr = y.pr + h * (25.0 / 216.0 * k1.pr + 1408.0 / 2565.0 * k3.pr + 2197.0 / 4104.0 * k4.pr - 1.0 / 5.0 * k5.pr);
        y4.ptheta = y.ptheta + h * (25.0 / 216.0 * k1.ptheta + 1408.0 / 2565.0 * k3.ptheta + 2197.0 / 4104.0 * k4.ptheta - 1.0 / 5.0 * k5.ptheta);
        y4.pphi = y.pphi + h * (25.0 / 216.0 * k1.pphi + 1408.0 / 2565.0 * k3.pphi + 2197.0 / 4104.0 * k4.pphi - 1.0 / 5.0 * k5.pphi);

        State y5;
        y5.t = y.t + h * (16.0 / 135.0 * k1.t + 6656.0 / 12825.0 * k3.t + 28561.0 / 56430.0 * k4.t - 9.0 / 50.0 * k5.t + 2.0 / 55.0 * k6.t);
        y5.r = y.r + h * (16.0 / 135.0 * k1.r + 6656.0 / 12825.0 * k3.r + 28561.0 / 56430.0 * k4.r - 9.0 / 50.0 * k5.r + 2.0 / 55.0 * k6.r);
        y5.theta = y.theta + h * (16.0 / 135.0 * k1.theta + 6656.0 / 12825.0 * k3.theta + 28561.0 / 56430.0 * k4.theta - 9.0 / 50.0 * k5.theta + 2.0 / 55.0 * k6.theta);
        y5.phi = y.phi + h * (16.0 / 135.0 * k1.phi + 6656.0 / 12825.0 * k3.phi + 28561.0 / 56430.0 * k4.phi - 9.0 / 50.0 * k5.phi + 2.0 / 55.0 * k6.phi);
        y5.pt = y.pt + h * (16.0 / 135.0 * k1.pt + 6656.0 / 12825.0 * k3.pt + 28561.0 / 56430.0 * k4.pt - 9.0 / 50.0 * k5.pt + 2.0 / 55.0 * k6.pt);
        y5.pr = y.pr + h * (16.0 / 135.0 * k1.pr + 6656.0 / 12825.0 * k3.pr + 28561.0 / 56430.0 * k4.pr - 9.0 / 50.0 * k5.pr + 2.0 / 55.0 * k6.pr);
        y5.ptheta = y.ptheta + h * (16.0 / 135.0 * k1.ptheta + 6656.0 / 12825.0 * k3.ptheta + 28561.0 / 56430.0 * k4.ptheta - 9.0 / 50.0 * k5.ptheta + 2.0 / 55.0 * k6.ptheta);
        y5.pphi = y.pphi + h * (16.0 / 135.0 * k1.pphi + 6656.0 / 12825.0 * k3.pphi + 28561.0 / 56430.0 * k4.pphi - 9.0 / 50.0 * k5.pphi + 2.0 / 55.0 * k6.pphi);

        const double err = error_norm(y4, y5);

        if (err < params.tolerance || h <= h_min) {
            y = y5;
            return;
        }

        const double safety = 0.8;
        const double factor =
            safety * pow(params.tolerance / fmax(err, 1.0e-30), 0.25);
        h *= clamp_device(factor, 0.1, 0.5);

        if (h < h_min) {
            h = h_min;
        }
    }

    rk4_step(params, y);
}

__device__ double torus_in(const Params& params, double r, double theta)
{
    (void)theta;
    return (r >= params.torus_r_min_rg && r <= params.torus_r_max_rg)
        ? 1.0
        : 0.0;
}

__device__ bool profile_uses_powerlaw(const Params& params)
{
    return
        params.density_profile_kind == 1 ||
        params.density_profile_kind == 3 ||
        params.density_profile_kind == 5 ||
        params.density_profile_kind == 6;
}

__device__ bool profile_uses_funnel(const Params& params)
{
    return
        params.density_profile_kind == 2 ||
        params.density_profile_kind == 3 ||
        params.density_profile_kind == 6;
}

__device__ bool profile_uses_envelope(const Params& params)
{
    return
        params.density_profile_kind == 4 ||
        params.density_profile_kind == 5 ||
        params.density_profile_kind == 6;
}

__device__ double gaussian_disk_shape(const Params& params, double r, double theta)
{
    const double delta = theta - 0.5 * PI;

    if (r <= 4.0 || r >= 18.0 || fabs(delta) >= 0.45) {
        return 0.0;
    }

    const double radial =
        exp(-pow((r - params.torus_r0_rg) / params.torus_sigma_rg, 2.0));

    const double vertical =
        exp(-pow(delta / params.torus_h_over_r, 2.0));

    return radial * vertical;
}

__device__ double powerlaw_disk_shape(const Params& params, double r, double theta)
{
    if (torus_in(params, r, theta) == 0.0 ||
        r <= 0.0 ||
        params.torus_r0_rg <= 0.0) {
        return 0.0;
    }

    const double vertical =
        exp(-pow(cos(theta) / params.torus_h_over_r, 2.0));

    const double inner_taper =
        1.0 - exp(-pow((r - params.torus_r_min_rg) / params.torus_sigma_rg, 2.0));

    const double outer_taper =
        exp(-pow(r / params.torus_r_max_rg, 4.0));

    const double radial =
        pow(r / params.torus_r0_rg, -params.torus_radial_power);

    return radial * vertical * inner_taper * outer_taper;
}

__device__ double funnel_factor(const Params& params, double theta)
{
    if (!profile_uses_funnel(params)) {
        return 1.0;
    }

    const double theta_f = fmax(params.funnel_theta_rad, 1.0e-6);

    const double north =
        exp(-pow(theta / theta_f, 2.0));

    const double south =
        exp(-pow((PI - theta) / theta_f, 2.0));

    return clamp_device(
        1.0 - params.funnel_depletion * (north + south),
        0.0,
        1.0
    );
}

__device__ double envelope_rho(const Params& params, double r)
{
    if (!profile_uses_envelope(params) ||
        params.envelope_rho0 <= 0.0 ||
        r < params.torus_r0_rg ||
        r > params.torus_r_max_rg) {
        return 0.0;
    }

    return params.envelope_rho0 * pow(r / params.torus_r0_rg, -params.envelope_alpha);
}

__device__ double torus_rho(const Params& params, double r, double theta)
{
    const double shape = profile_uses_powerlaw(params)
        ? powerlaw_disk_shape(params, r, theta)
        : gaussian_disk_shape(params, r, theta);

    const double disk =
        params.torus_rho0 * shape * funnel_factor(params, theta);

    return fmax(disk + envelope_rho(params, r), params.rho_floor);
}

__device__ double torus_temperature_MeV(const Params& params, double r, double theta)
{
    const double shape =
        clamp_device(torus_rho(params, r, theta) / fmax(params.torus_rho0, 1.0e-300), 0.0, 1.0);

    const double T = 6.0 * pow(shape, 0.2);
    return fmax(T, 1.0e-10);
}

__device__ double torus_Ye(const Params& params, double r, double theta)
{
    if (torus_rho(params, r, theta) <= 0.0) {
        return 0.0;
    }

    return 0.2 + 0.1 * exp(-pow((r - params.torus_r0_rg) / 4.0, 2.0));
}

__device__ double sigma_interp_cm2(
    double E_GeV,
    const double* sigma_E,
    const double* sigma_cm2,
    int n
)
{
    int lo = 0;
    int hi = n - 1;

    while (hi - lo > 1) {
        const int mid = (lo + hi) / 2;
        if (sigma_E[mid] <= E_GeV) {
            lo = mid;
        } else {
            hi = mid;
        }
    }

    const double lx = log(E_GeV);
    const double lx1 = log(sigma_E[lo]);
    const double lx2 = log(sigma_E[hi]);
    const double ly1 = log(sigma_cm2[lo]);
    const double ly2 = log(sigma_cm2[hi]);
    const double t = (lx - lx1) / (lx2 - lx1);
    return exp(ly1 + t * (ly2 - ly1));
}

__device__ double fermi_dirac_shape(double Enu_MeV, double T_MeV)
{
    if (Enu_MeV <= 0.0 || T_MeV <= 0.0) {
        return 0.0;
    }

    const double x = Enu_MeV / T_MeV;
    if (x > 120.0) {
        return 0.0;
    }

    return Enu_MeV * Enu_MeV / (exp(x) + 1.0);
}

__device__ double emissivity_collapsar_ring(
    const Params& params,
    double r,
    double theta,
    double Enu_local_GeV
)
{
    if (Enu_local_GeV <= 0.0) {
        return 0.0;
    }

    const double delta_r = (r - params.source_r_center_rg) / params.source_sigma_r_rg;
    const double delta_theta = (theta - 0.5 * PI) / params.source_theta_width_rad;
    const double radial_profile = exp(-delta_r * delta_r);
    const double vertical_profile = exp(-delta_theta * delta_theta);
    const double spectral_profile =
        pow(fmax(Enu_local_GeV, 1.0e-300), -params.source_powerlaw)
        * exp(-Enu_local_GeV / params.source_emax_GeV);

    return params.source_norm * radial_profile * vertical_profile * spectral_profile;
}

__device__ double spectral_uhe(const Params& params, double Enu_local_GeV)
{
    if (Enu_local_GeV <= 0.0 || params.source_emax_GeV <= 0.0) {
        return 0.0;
    }

    return
        pow(fmax(Enu_local_GeV, 1.0e-300), -params.source_powerlaw)
        * exp(-Enu_local_GeV / params.source_emax_GeV);
}

__device__ double bipolar_gaussian_source(
    double theta,
    double theta0,
    double width
)
{
    if (width <= 0.0) {
        return 0.0;
    }

    const double north = theta - theta0;
    const double south = theta - (PI - theta0);

    return
        exp(-(north / width) * (north / width))
        + exp(-(south / width) * (south / width));
}

__device__ double uhe_source_spatial_weight(
    const Params& params,
    double r,
    double theta
)
{
    if (r <= 0.0 ||
        params.source_sigma_r_rg <= 0.0 ||
        params.source_theta_width_rad <= 0.0) {
        return 0.0;
    }

    if (params.source_model_kind == SOURCE_INNER_RING) {
        const double delta_r =
            (r - params.source_r_center_rg) / params.source_sigma_r_rg;
        const double delta_theta =
            (theta - 0.5 * PI) / params.source_theta_width_rad;

        return exp(-delta_r * delta_r) * exp(-delta_theta * delta_theta);
    }

    if (params.source_model_kind == SOURCE_FUNNEL_WALL) {
        const double delta_r =
            (r - params.source_r_center_rg) / params.source_sigma_r_rg;

        return exp(-delta_r * delta_r)
            * bipolar_gaussian_source(
                theta,
                params.source_funnel_theta_rad,
                params.source_theta_width_rad
            );
    }

    if (params.source_model_kind == SOURCE_JET_BASE) {
        const double radial =
            exp(-pow(r / fmax(params.source_r_center_rg, 1.0e-300), 2.0));

        const double polar =
            bipolar_gaussian_source(theta, 0.0, params.source_theta_width_rad);

        return radial * polar;
    }

    if (params.source_model_kind == SOURCE_DENSITY_WEIGHTED) {
        const double rho = torus_rho(params, r, theta);
        const double rho_norm =
            fmax(
                params.source_rho_ref_gcm3 > 0.0
                ? params.source_rho_ref_gcm3
                : torus_rho(params, params.source_r_center_rg, 0.5 * PI),
                1.0e-300
            );

        const double density_weight =
            pow(fmax(rho / rho_norm, 0.0), params.source_density_power_q);

        const double radial_weight =
            pow(
                fmax(r / fmax(params.source_r_center_rg, 1.0e-300), 1.0e-300),
                -params.source_radial_power_s
            );

        const double cutoff_min = fmax(params.source_cutoff_min, 0.0);
        const double cutoff_max = fmax(params.source_cutoff_max, cutoff_min);

        return clamp_device(density_weight * radial_weight, cutoff_min, cutoff_max);
    }

    if (params.source_model_kind == SOURCE_SHOCK_LAYER) {
        const double dr = fmax(params.source_gradient_dr_rg, 1.0e-4);
        const double dtheta = fmax(params.source_gradient_dtheta_rad, 1.0e-5);
        const double r_minus = fmax(r - dr, 1.0e-6);
        const double theta_minus = clamp_device(theta - dtheta, 0.0, PI);
        const double theta_plus = clamp_device(theta + dtheta, 0.0, PI);

        const double rho_r_plus = torus_rho(params, r + dr, theta);
        const double rho_r_minus = torus_rho(params, r_minus, theta);
        const double rho_t_plus = torus_rho(params, r, theta_plus);
        const double rho_t_minus = torus_rho(params, r, theta_minus);
        const double rho_norm =
            fmax(torus_rho(params, params.source_r_center_rg, 0.5 * PI), 1.0e-300);

        const double d_r = (rho_r_plus - rho_r_minus) / (2.0 * dr);
        const double d_theta =
            (rho_t_plus - rho_t_minus) / (2.0 * dtheta * fmax(r, 1.0e-6));

        const double gradient = sqrt(d_r * d_r + d_theta * d_theta) / rho_norm;
        const double radial_taper =
            exp(-pow(r / fmax(params.source_r_center_rg, 1.0e-300), 2.0));

        const double cutoff_min = fmax(params.source_cutoff_min, 0.0);
        const double cutoff_max = fmax(params.source_cutoff_max, cutoff_min);

        return clamp_device(gradient * radial_taper, cutoff_min, cutoff_max);
    }

    return 0.0;
}

__device__ double emissivity_uhe(
    const Params& params,
    double r,
    double theta,
    double Enu_local_GeV
)
{
    if (params.source_model_kind == SOURCE_INNER_RING) {
        return emissivity_collapsar_ring(params, r, theta, Enu_local_GeV);
    }

    return params.source_norm
        * uhe_source_spatial_weight(params, r, theta)
        * spectral_uhe(params, Enu_local_GeV);
}

__device__ double emissivity_mev_thermal(
    const Params& params,
    double rho,
    double T_MeV,
    double Ye,
    double Enu_local_MeV
)
{
    if (rho <= 0.0 || T_MeV <= 0.0 || Enu_local_MeV <= 0.0) {
        return 0.0;
    }

    const double nb = rho / M_U_G;
    const double neutron_fraction = clamp_device(1.0 - Ye, 0.0, 1.0);
    const double proton_fraction = clamp_device(Ye, 0.0, 1.0);
    const double charged_current_weight = neutron_fraction + 0.5 * proton_fraction;
    const double thermal_spectrum = fermi_dirac_shape(Enu_local_MeV, T_MeV);
    const double capture_like = nb * charged_current_weight * pow(T_MeV, 6.0);
    const double pair_like = 0.05 * nb * pow(T_MeV, 9.0);

    return params.mev_norm * (capture_like + pair_like) * thermal_spectrum;
}

__device__ double opacity_mev_absorption_cm_inv(
    const Params& params,
    double rho,
    double Ye,
    double Enu_local_MeV
)
{
    if (rho <= 0.0 || Enu_local_MeV <= 0.0) {
        return 0.0;
    }

    const double nb = rho / M_U_G;
    const double neutron_fraction = clamp_device(1.0 - Ye, 0.0, 1.0);
    const double proton_fraction = clamp_device(Ye, 0.0, 1.0);
    const double target_fraction = neutron_fraction + 0.25 * proton_fraction;

    return nb * params.mev_sigma_abs0_cm2 * Enu_local_MeV * Enu_local_MeV * target_fraction;
}

__device__ double opacity_mev_scattering_cm_inv(
    const Params& params,
    double rho,
    double Ye,
    double Enu_local_MeV
)
{
    if (rho <= 0.0 || Enu_local_MeV <= 0.0) {
        return 0.0;
    }

    const double nb = rho / M_U_G;
    const double composition_factor = 1.0 + 0.5 * clamp_device(Ye, 0.0, 1.0);
    return nb * params.mev_sigma_scat0_cm2 * Enu_local_MeV * Enu_local_MeV * composition_factor;
}

__device__ void accumulate_radiative_transfer(
    const Params& params,
    const double* sigma_E,
    const double* sigma_cm2,
    RayAccum& ray
)
{
    const State& y = ray.y;
    const double r_h = kerr_horizon(params.a_spin);

    if (y.r < r_h + 1.0e-3) {
        return;
    }

    const double dl_rg = zamo_spatial_interval_rg(params, y, ray);
    if (dl_rg <= 0.0) {
        return;
    }

    const double local_over_obs =
        zamo_energy(params, y.r, y.theta, y.pt, y.pphi)
        / fmax(ray.energy_obs, 1.0e-300);

    const double g_obs_over_local = 1.0 / fmax(local_over_obs, 1.0e-300);
    const double Enu_local_GeV = params.Enu_obs_GeV * local_over_obs;
    const double Enu_local_MeV = params.mev_Enu_obs_MeV * local_over_obs;
    const double dl_cm = dl_rg * rg_cm(params.M_bh_msun);
    const double rho = torus_rho(params, y.r, y.theta);
    const double T_MeV = torus_temperature_MeV(params, y.r, y.theta);
    const double Ye = torus_Ye(params, y.r, y.theta);

    const double kappa_mev =
        opacity_mev_absorption_cm_inv(params, rho, Ye, Enu_local_MeV)
        + opacity_mev_scattering_cm_inv(params, rho, Ye, Enu_local_MeV);

    const double leakage = 1.0 / (1.0 + ray.tau_mev * ray.tau_mev);
    const double j_mev =
        emissivity_mev_thermal(params, rho, T_MeV, Ye, Enu_local_MeV);

    ray.I_mev +=
        pow(g_obs_over_local, 3.0)
        * j_mev
        * leakage
        * exp(-ray.tau_mev)
        * dl_cm;

    ray.leakage_sum += leakage * j_mev * dl_cm;
    ray.leakage_weight += j_mev * dl_cm;
    ray.tau_mev += kappa_mev * dl_cm;

    if (ray.r_neutrinosphere_rg < 0.0 &&
        ray.tau_mev >= params.mev_neutrinosphere_tau) {
        ray.r_neutrinosphere_rg = y.r;
    }

    const double j_uhe =
        emissivity_uhe(params, y.r, y.theta, Enu_local_GeV);

    ray.I_uhe +=
        pow(g_obs_over_local, 3.0)
        * j_uhe
        * exp(-ray.tau_uhe)
        * dl_cm;

    if (rho <= 0.0 ||
        Enu_local_GeV < sigma_E[0] ||
        Enu_local_GeV > sigma_E[params.sigma_n - 1]) {
        return;
    }

    const double sigma = sigma_interp_cm2(
        Enu_local_GeV,
        sigma_E,
        sigma_cm2,
        params.sigma_n
    );

    const double nb = rho / M_U_G;
    ray.tau_uhe += nb * sigma * dl_cm;
}

__device__ void accumulate_cache_point(
    const Params& params,
    const double* sigma_E,
    const double* sigma_cm2,
    const CachePoint& p,
    OutputPixel& out
)
{
    const double r_h = kerr_horizon(params.a_spin);

    if (p.r_rg < r_h + 1.0e-3 || p.dl_rg <= 0.0) {
        return;
    }

    const double local_over_obs =
        fmax(p.redshift_factor, 1.0e-300);

    const double g_obs_over_local =
        1.0 / local_over_obs;

    const double Enu_local_GeV =
        params.Enu_obs_GeV * local_over_obs;

    const double Enu_local_MeV =
        params.mev_Enu_obs_MeV * local_over_obs;

    const double dl_cm =
        p.dl_rg * rg_cm(params.M_bh_msun);

    const double rho =
        torus_rho(params, p.r_rg, p.theta);

    const double T_MeV =
        torus_temperature_MeV(params, p.r_rg, p.theta);

    const double Ye =
        torus_Ye(params, p.r_rg, p.theta);

    const double kappa_mev =
        opacity_mev_absorption_cm_inv(params, rho, Ye, Enu_local_MeV)
        + opacity_mev_scattering_cm_inv(params, rho, Ye, Enu_local_MeV);

    const double leakage =
        1.0 / (1.0 + out.tau_mev * out.tau_mev);

    const double j_mev =
        emissivity_mev_thermal(params, rho, T_MeV, Ye, Enu_local_MeV);

    out.I_mev +=
        pow(g_obs_over_local, 3.0)
        * j_mev
        * leakage
        * exp(-out.tau_mev)
        * dl_cm;

    out.leakage_mev += leakage * j_mev * dl_cm;
    out.P_surv_mev += j_mev * dl_cm;
    out.tau_mev += kappa_mev * dl_cm;

    if (out.r_neutrinosphere_rg < 0.0 &&
        out.tau_mev >= params.mev_neutrinosphere_tau) {
        out.r_neutrinosphere_rg = p.r_rg;
    }

    const double j_uhe =
        emissivity_uhe(params, p.r_rg, p.theta, Enu_local_GeV);

    out.I_uhe +=
        pow(g_obs_over_local, 3.0)
        * j_uhe
        * exp(-out.tau_uhe)
        * dl_cm;

    if (rho <= 0.0 ||
        Enu_local_GeV < sigma_E[0] ||
        Enu_local_GeV > sigma_E[params.sigma_n - 1]) {
        return;
    }

    const double sigma =
        sigma_interp_cm2(
            Enu_local_GeV,
            sigma_E,
            sigma_cm2,
            params.sigma_n
        );

    const double nb = rho / M_U_G;
    out.tau_uhe += nb * sigma * dl_cm;
}

__global__ void radiative_from_cache_kernel(
    Params params,
    const double* sigma_E,
    const double* sigma_cm2,
    const CacheRayDevice* rays,
    const CachePoint* points,
    OutputPixel* outputs,
    int n_rays
)
{
    const int id = blockIdx.x * blockDim.x + threadIdx.x;
    if (id >= n_rays) {
        return;
    }

    const CacheRayDevice ray = rays[id];
    OutputPixel out{};
    out.pixel_i = ray.pixel_i;
    out.pixel_j = ray.pixel_j;
    out.captured = ray.captured;
    out.alpha_rg = ray.alpha_rg;
    out.beta_rg = ray.beta_rg;
    out.P_surv_mev = 0.0;
    out.r_neutrinosphere_rg = -1.0;

    for (int k = 0; k < ray.npoints; ++k) {
        accumulate_cache_point(
            params,
            sigma_E,
            sigma_cm2,
            points[ray.offset + k],
            out
        );
    }

    const double leakage_weight = out.P_surv_mev;
    out.P_surv_uhe = exp(-out.tau_uhe);
    out.P_surv_mev = exp(-out.tau_mev);
    out.leakage_mev =
        leakage_weight > 0.0 ? out.leakage_mev / leakage_weight : 1.0;

    outputs[id] = out;
}

__global__ void init_rays_kernel(
    Params params,
    int tile_i0,
    int tile_j0,
    int tile_nx,
    int tile_ny,
    RayAccum* rays,
    int* active_count
)
{
    const int local_id = blockIdx.x * blockDim.x + threadIdx.x;
    const int n = tile_nx * tile_ny;
    if (local_id >= n) {
        return;
    }

    const int local_i = local_id % tile_nx;
    const int local_j = local_id / tile_nx;
    const int i = tile_i0 + local_i;
    const int j = tile_j0 + local_j;

    RayAccum ray{};
    ray.pixel_i = i;
    ray.pixel_j = j;
    ray.alpha_rg = (2.0 * (i + 0.5) / params.nx - 1.0) * tan(0.5 * params.fov_rad);
    ray.beta_rg = (2.0 * (j + 0.5) / params.ny - 1.0) * tan(0.5 * params.fov_rad);
    ray.captured = 0;
    ray.finished = 0;
    ray.steps_done = 0;
    ray.r_neutrinosphere_rg = -1.0;

    const double u = ray.alpha_rg;
    const double v = ray.beta_rg;
    const double norm = sqrt(1.0 + u * u + v * v);
    const double n_r = -1.0 / norm;
    const double n_theta = v / norm;
    const double n_phi = u / norm;

    double g[4][4];
    kerr_metric(params.r_obs_rg, params.theta_obs_rad, params.a_spin, g);

    const double alpha = kerr_lapse(params.r_obs_rg, params.theta_obs_rad, params.a_spin);
    const double omega = kerr_omega(params.r_obs_rg, params.theta_obs_rad, params.a_spin);
    const double p_contra[4] = {
        1.0 / alpha,
        n_r / sqrt(g[1][1]),
        n_theta / sqrt(g[2][2]),
        n_phi / sqrt(g[3][3]) + omega / alpha
    };

    double p_cov[4] = {0.0, 0.0, 0.0, 0.0};
    for (int mu = 0; mu < 4; ++mu) {
        for (int nu = 0; nu < 4; ++nu) {
            p_cov[mu] += g[mu][nu] * p_contra[nu];
        }
    }

    ray.y.t = 0.0;
    ray.y.r = params.r_obs_rg;
    ray.y.theta = params.theta_obs_rad;
    ray.y.phi = 0.0;
    ray.y.pt = p_cov[0];
    ray.y.pr = p_cov[1];
    ray.y.ptheta = p_cov[2];
    ray.y.pphi = p_cov[3];
    ray.prev_r = ray.y.r;
    ray.prev_theta = ray.y.theta;
    ray.prev_phi = ray.y.phi;
    ray.energy_obs =
        zamo_energy(params, params.r_obs_rg, params.theta_obs_rad, ray.y.pt, ray.y.pphi);

    rays[local_id] = ray;
    atomicAdd(active_count, 1);
}

__global__ void init_linear_rays_kernel(
    Params params,
    int ray_start,
    int n_rays,
    RayAccum* rays
)
{
    const int local_id = blockIdx.x * blockDim.x + threadIdx.x;
    if (local_id >= n_rays) {
        return;
    }

    const int ray_id = ray_start + local_id;
    const int i = ray_id / params.ny;
    const int j = ray_id - i * params.ny;

    RayAccum ray{};
    ray.pixel_i = i;
    ray.pixel_j = j;
    ray.alpha_rg = (2.0 * (i + 0.5) / params.nx - 1.0) * tan(0.5 * params.fov_rad);
    ray.beta_rg = (2.0 * (j + 0.5) / params.ny - 1.0) * tan(0.5 * params.fov_rad);
    ray.r_neutrinosphere_rg = -1.0;

    const double u = ray.alpha_rg;
    const double v = ray.beta_rg;
    const double norm = sqrt(1.0 + u * u + v * v);
    const double n_r = -1.0 / norm;
    const double n_theta = v / norm;
    const double n_phi = u / norm;

    double g[4][4];
    kerr_metric(params.r_obs_rg, params.theta_obs_rad, params.a_spin, g);

    const double alpha = kerr_lapse(params.r_obs_rg, params.theta_obs_rad, params.a_spin);
    const double omega = kerr_omega(params.r_obs_rg, params.theta_obs_rad, params.a_spin);
    const double p_contra[4] = {
        1.0 / alpha,
        n_r / sqrt(g[1][1]),
        n_theta / sqrt(g[2][2]),
        n_phi / sqrt(g[3][3]) + omega / alpha
    };

    double p_cov[4] = {0.0, 0.0, 0.0, 0.0};
    for (int mu = 0; mu < 4; ++mu) {
        for (int nu = 0; nu < 4; ++nu) {
            p_cov[mu] += g[mu][nu] * p_contra[nu];
        }
    }

    ray.y.t = 0.0;
    ray.y.r = params.r_obs_rg;
    ray.y.theta = params.theta_obs_rad;
    ray.y.phi = 0.0;
    ray.y.pt = p_cov[0];
    ray.y.pr = p_cov[1];
    ray.y.ptheta = p_cov[2];
    ray.y.pphi = p_cov[3];
    ray.prev_r = ray.y.r;
    ray.prev_theta = ray.y.theta;
    ray.prev_phi = ray.y.phi;
    ray.energy_obs =
        zamo_energy(params, params.r_obs_rg, params.theta_obs_rad, ray.y.pt, ray.y.pphi);

    rays[local_id] = ray;
}

__global__ void trace_cache_kernel(
    Params params,
    RayAccum* rays,
    CachePoint* points,
    int* counts,
    int n_rays,
    int max_points
)
{
    const int id = blockIdx.x * blockDim.x + threadIdx.x;
    if (id >= n_rays) {
        return;
    }

    RayAccum ray = rays[id];
    const double r_stop = kerr_horizon(params.a_spin) + 1.0e-3;
    int npoints = 0;

    for (int step = 0; step < params.max_steps && npoints < max_points; ++step) {
        if (ray.y.r <= r_stop) {
            ray.captured = 1;
            ray.finished = 1;
            break;
        }

        if (ray.y.r >= params.r_max_rg && step > 10) {
            ray.finished = 1;
            break;
        }

        CachePoint p{};
        p.r_rg = ray.y.r;
        p.theta = ray.y.theta;
        p.x_rg = ray.y.r * sin(ray.y.theta) * cos(ray.y.phi);
        p.y_rg = ray.y.r * sin(ray.y.theta) * sin(ray.y.phi);
        p.z_rg = ray.y.r * cos(ray.y.theta);
        p.dl_rg = (step == 0) ? 0.0 : zamo_spatial_interval_rg(params, ray.y, ray);
        p.redshift_factor =
            zamo_energy(params, ray.y.r, ray.y.theta, ray.y.pt, ray.y.pphi)
            / fmax(ray.energy_obs, 1.0e-300);

        points[id * max_points + npoints] = p;
        npoints += 1;

        ray.prev_r = ray.y.r;
        ray.prev_theta = ray.y.theta;
        ray.prev_phi = ray.y.phi;
        rkf45_step_adaptive(params, ray.y);
        ray.steps_done += 1;
    }

    if (!ray.finished && npoints >= max_points) {
        ray.finished = 1;
    }

    counts[id] = npoints;
    rays[id] = ray;
}

__global__ void advance_rays_kernel(
    Params params,
    const double* sigma_E,
    const double* sigma_cm2,
    RayAccum* rays,
    int n_rays,
    int steps_per_launch,
    int* active_count
)
{
    const int id = blockIdx.x * blockDim.x + threadIdx.x;
    if (id >= n_rays) {
        return;
    }

    RayAccum ray = rays[id];
    if (ray.finished) {
        return;
    }

    const double r_h = kerr_horizon(params.a_spin);
    const double r_stop = r_h + 1.0e-3;

    for (int step = 0; step < steps_per_launch; ++step) {
        if (ray.y.r <= r_stop) {
            ray.captured = 1;
            ray.finished = 1;
            break;
        }

        if (ray.y.r >= params.r_max_rg && ray.steps_done > 10) {
            ray.finished = 1;
            break;
        }

        if (ray.steps_done >= params.max_steps) {
            ray.finished = 1;
            break;
        }

        if (ray.steps_done > 0) {
            accumulate_radiative_transfer(params, sigma_E, sigma_cm2, ray);
        }

        ray.prev_r = ray.y.r;
        ray.prev_theta = ray.y.theta;
        ray.prev_phi = ray.y.phi;
        rkf45_step_adaptive(params, ray.y);
        ray.steps_done += 1;
    }

    if (!ray.finished) {
        atomicAdd(active_count, 1);
    }

    rays[id] = ray;
}

void write_tile_outputs(
    const std::vector<RayAccum>& rays,
    std::vector<OutputPixel>& image,
    int nx
)
{
    for (const RayAccum& ray : rays) {
        OutputPixel out{};
        out.pixel_i = ray.pixel_i;
        out.pixel_j = ray.pixel_j;
        out.captured = ray.captured;
        out.alpha_rg = ray.alpha_rg;
        out.beta_rg = ray.beta_rg;
        out.tau_uhe = ray.tau_uhe;
        out.P_surv_uhe = std::exp(-ray.tau_uhe);
        out.I_uhe = ray.I_uhe;
        out.tau_mev = ray.tau_mev;
        out.P_surv_mev = std::exp(-ray.tau_mev);
        out.I_mev = ray.I_mev;
        out.r_neutrinosphere_rg = ray.r_neutrinosphere_rg;
        out.leakage_mev =
            ray.leakage_weight > 0.0
                ? ray.leakage_sum / ray.leakage_weight
                : 1.0;

        image[static_cast<std::size_t>(out.pixel_j * nx + out.pixel_i)] = out;
    }
}

int run_geodesic_cache_gpu(int argc, char* argv[])
{
    Params params{};
    params.a_spin = 0.9;
    params.r_obs_rg = 60.0;
    params.theta_obs_rad = 80.0 * PI / 180.0;
    params.fov_rad = 25.0 * PI / 180.0;
    params.nx = 100;
    params.ny = 100;
    params.r_max_rg = 120.0;
    params.h = 0.001;
    params.tolerance = 1.0e-8;
    params.max_steps = 200000;

    int cache_rays = 32;
    int cache_max_points = 50000;

    if (argc > 2) params.a_spin = std::atof(argv[2]);
    if (argc > 3) params.r_obs_rg = std::atof(argv[3]);
    if (argc > 4) params.theta_obs_rad = std::atof(argv[4]) * PI / 180.0;
    if (argc > 5) params.fov_rad = std::atof(argv[5]) * PI / 180.0;
    if (argc > 6) params.nx = std::atoi(argv[6]);
    if (argc > 7) params.ny = std::atoi(argv[7]);
    if (argc > 8) params.r_max_rg = std::atof(argv[8]);
    if (argc > 9) params.h = std::atof(argv[9]);
    if (argc > 10) cache_rays = std::atoi(argv[10]);
    if (argc > 11) cache_max_points = std::atoi(argv[11]);
    if (argc > 12) params.tolerance = std::atof(argv[12]);

    if (params.nx <= 0 || params.ny <= 0 || cache_rays <= 0 ||
        cache_max_points <= 0) {
        std::cerr
            << "Usage: " << argv[0]
            << " --geodesic-cache ASPIN CAM_R_OBS CAM_THETA_DEG CAM_FOV_DEG"
            << " CAM_NX CAM_NY CAM_R_MAX CAM_STEP GPU_CACHE_RAYS"
            << " GPU_CACHE_MAX_POINTS GPU_TOL\n";
        return 1;
    }

    int device_count = 0;
    CUDA_CHECK(cudaGetDeviceCount(&device_count));
    if (device_count <= 0) {
        throw std::runtime_error("No CUDA device found.");
    }

    cudaDeviceProp prop{};
    CUDA_CHECK(cudaGetDeviceProperties(&prop, 0));

    std::cout << "CUDA device 0: " << prop.name << "\n";
    std::cout << "Geodesic cache grid: " << params.nx << " x " << params.ny << "\n";
    std::cout << "Batch rays: " << cache_rays
              << ", max points/ray: " << cache_max_points
              << ", RKF45 tol: " << params.tolerance << "\n";

    std::ofstream out(
        "output/rays/kerr_geodesics_cuda.bin",
        std::ios::binary
    );

    if (!out) {
        throw std::runtime_error("Could not open output/rays/kerr_geodesics_cuda.bin");
    }

    const std::int32_t magic = 0x4B47454F; // "KGEO"
    const std::int32_t version = 1;
    const std::int32_t nx_i32 = static_cast<std::int32_t>(params.nx);
    const std::int32_t ny_i32 = static_cast<std::int32_t>(params.ny);

    out.write(reinterpret_cast<const char*>(&magic), sizeof(magic));
    out.write(reinterpret_cast<const char*>(&version), sizeof(version));
    out.write(reinterpret_cast<const char*>(&nx_i32), sizeof(nx_i32));
    out.write(reinterpret_cast<const char*>(&ny_i32), sizeof(ny_i32));
    out.write(reinterpret_cast<const char*>(&params.a_spin), sizeof(params.a_spin));

    RayAccum* d_rays = nullptr;
    CachePoint* d_points = nullptr;
    int* d_counts = nullptr;

    CUDA_CHECK(cudaMalloc(&d_rays, cache_rays * sizeof(RayAccum)));
    CUDA_CHECK(cudaMalloc(
        &d_points,
        static_cast<std::size_t>(cache_rays)
        * static_cast<std::size_t>(cache_max_points)
        * sizeof(CachePoint)
    ));
    CUDA_CHECK(cudaMalloc(&d_counts, cache_rays * sizeof(int)));

    std::vector<RayAccum> rays(static_cast<std::size_t>(cache_rays));
    std::vector<int> counts(static_cast<std::size_t>(cache_rays));
    std::vector<CachePoint> points(
        static_cast<std::size_t>(cache_rays)
        * static_cast<std::size_t>(cache_max_points)
    );

    const int threads = 128;
    const int total_rays = params.nx * params.ny;
    std::uint64_t total_points = 0;
    std::uint64_t truncated_rays = 0;
    const auto t_start = std::chrono::steady_clock::now();

    for (int ray_start = 0; ray_start < total_rays; ray_start += cache_rays) {
        const int n_batch = std::min(cache_rays, total_rays - ray_start);
        const int blocks = (n_batch + threads - 1) / threads;

        init_linear_rays_kernel<<<blocks, threads>>>(
            params,
            ray_start,
            n_batch,
            d_rays
        );
        CUDA_CHECK(cudaGetLastError());
        CUDA_CHECK(cudaDeviceSynchronize());

        trace_cache_kernel<<<blocks, threads>>>(
            params,
            d_rays,
            d_points,
            d_counts,
            n_batch,
            cache_max_points
        );
        CUDA_CHECK(cudaGetLastError());
        CUDA_CHECK(cudaDeviceSynchronize());

        CUDA_CHECK(cudaMemcpy(
            rays.data(),
            d_rays,
            n_batch * sizeof(RayAccum),
            cudaMemcpyDeviceToHost
        ));
        CUDA_CHECK(cudaMemcpy(
            counts.data(),
            d_counts,
            n_batch * sizeof(int),
            cudaMemcpyDeviceToHost
        ));
        CUDA_CHECK(cudaMemcpy(
            points.data(),
            d_points,
            static_cast<std::size_t>(n_batch)
            * static_cast<std::size_t>(cache_max_points)
            * sizeof(CachePoint),
            cudaMemcpyDeviceToHost
        ));

        for (int k = 0; k < n_batch; ++k) {
            const int ray_id = ray_start + k;
            BinaryRayHeader header{};
            header.ray_id = static_cast<std::int32_t>(ray_id);
            header.pixel_i = static_cast<std::int32_t>(rays[k].pixel_i);
            header.pixel_j = static_cast<std::int32_t>(rays[k].pixel_j);
            header.captured = rays[k].captured ? 1 : 0;
            header.npoints = static_cast<std::int32_t>(counts[k]);
            header.alpha_rg = rays[k].alpha_rg;
            header.beta_rg = rays[k].beta_rg;

            out.write(reinterpret_cast<const char*>(&header), sizeof(header));
            out.write(
                reinterpret_cast<const char*>(&points[k * cache_max_points]),
                static_cast<std::streamsize>(counts[k] * sizeof(CachePoint))
            );
            total_points += static_cast<std::uint64_t>(counts[k]);

            if (counts[k] >= cache_max_points && !rays[k].captured) {
                truncated_rays += 1;
                std::cerr
                    << "Warning: ray " << ray_id
                    << " reached GPU_CACHE_MAX_POINTS=" << cache_max_points
                    << " and may be truncated.\n";
            }
        }

        std::cout << "Cached rays "
                  << ray_start << "-"
                  << (ray_start + n_batch - 1) << "\n";
    }

    const auto t_end = std::chrono::steady_clock::now();
    const double elapsed_s =
        std::chrono::duration<double>(t_end - t_start).count();

    CUDA_CHECK(cudaFree(d_rays));
    CUDA_CHECK(cudaFree(d_points));
    CUDA_CHECK(cudaFree(d_counts));

    std::cout << "Saved: output/rays/kerr_geodesics_cuda.bin\n";
    std::cout << "Cached points: " << total_points << "\n";
    std::cout << "Geodesic GPU timing: " << elapsed_s << " s, "
              << (static_cast<double>(total_rays) / std::max(elapsed_s, 1.0e-300))
              << " rays/s\n";

    if (truncated_rays > 0) {
        std::cerr
            << "ERROR: " << truncated_rays
            << " rays reached GPU_CACHE_MAX_POINTS. "
            << "Do not use this cache for science. Increase GPU_CACHE_MAX_POINTS "
            << "or CAM_STEP, then regenerate output/rays/kerr_geodesics_cuda.bin.\n";

        return 2;
    }

    return 0;
}

int run_image_from_cache_gpu(int argc, char* argv[])
{
    Params params{};
    params.Enu_obs_GeV = 1.0e9;
    params.a_spin = 0.9;
    params.M_bh_msun = 3.0;
    params.r_obs_rg = 60.0;
    params.theta_obs_rad = 80.0 * PI / 180.0;
    params.torus_rho0 = 1.0e-2;
    params.torus_r0_rg = 10.0;
    params.torus_sigma_rg = 5.0;
    params.torus_h_over_r = 0.25;
    params.density_profile_kind = 0;
    params.torus_radial_power = 2.0;
    params.funnel_depletion = 0.0;
    params.funnel_theta_rad = 15.0 * PI / 180.0;
    params.envelope_rho0 = 0.0;
    params.envelope_alpha = 2.5;
    params.torus_r_min_rg = 4.0;
    params.torus_r_max_rg = 60.0;
    params.rho_floor = 1.0e-99;
    params.source_r_center_rg = 3.5;
    params.source_sigma_r_rg = 1.0;
    params.source_theta_width_rad = 15.0 * PI / 180.0;
    params.source_powerlaw = 2.0;
    params.source_emax_GeV = 1.0e12;
    params.source_norm = 1.0;
    params.source_model_kind = SOURCE_INNER_RING;
    params.source_funnel_theta_rad = 20.0 * PI / 180.0;
    params.source_density_power_q = 1.0;
    params.source_radial_power_s = 2.0;
    params.source_rho_ref_gcm3 = -1.0;
    params.source_cutoff_min = 0.0;
    params.source_cutoff_max = 1.0e2;
    params.source_gradient_dr_rg = 0.1;
    params.source_gradient_dtheta_rad = 1.0 * PI / 180.0;
    params.mev_Enu_obs_MeV = 10.0;
    params.mev_norm = 1.0;
    params.mev_sigma_abs0_cm2 = 9.6e-44;
    params.mev_sigma_scat0_cm2 = 1.7e-44;
    params.mev_neutrinosphere_tau = 2.0 / 3.0;

    double cam_theta_deg = 80.0;
    int chunk_rays = 256;
    std::string density_profile = "gaussian";
    std::string source_model = "inner_ring";
    double funnel_theta_deg = 15.0;

    if (argc > 2) params.Enu_obs_GeV = std::atof(argv[2]);
    if (argc > 3) params.a_spin = std::atof(argv[3]);
    if (argc > 4) params.M_bh_msun = std::atof(argv[4]);
    if (argc > 5) params.torus_rho0 = std::atof(argv[5]);
    if (argc > 6) params.torus_r0_rg = std::atof(argv[6]);
    if (argc > 7) params.torus_sigma_rg = std::atof(argv[7]);
    if (argc > 8) params.torus_h_over_r = std::atof(argv[8]);
    if (argc > 9) params.source_r_center_rg = std::atof(argv[9]);
    if (argc > 10) params.source_sigma_r_rg = std::atof(argv[10]);
    if (argc > 11) params.source_theta_width_rad = std::atof(argv[11]) * PI / 180.0;
    if (argc > 12) params.source_powerlaw = std::atof(argv[12]);
    if (argc > 13) params.source_emax_GeV = std::atof(argv[13]);
    if (argc > 14) params.source_norm = std::atof(argv[14]);
    if (argc > 15) params.mev_Enu_obs_MeV = std::atof(argv[15]);
    if (argc > 16) params.mev_norm = std::atof(argv[16]);
    if (argc > 17) cam_theta_deg = std::atof(argv[17]);
    if (argc > 18) chunk_rays = std::atoi(argv[18]);
    if (argc > 19) {
        density_profile = argv[19];
        params.density_profile_kind =
            parse_density_profile_kind(density_profile);
    }
    if (argc > 20) params.torus_radial_power = std::atof(argv[20]);
    if (argc > 21) params.funnel_depletion = std::atof(argv[21]);
    if (argc > 22) {
        funnel_theta_deg = std::atof(argv[22]);
        params.funnel_theta_rad = funnel_theta_deg * PI / 180.0;
    }
    if (argc > 23) params.envelope_rho0 = std::atof(argv[23]);
    if (argc > 24) params.envelope_alpha = std::atof(argv[24]);
    if (argc > 25) params.torus_r_min_rg = std::atof(argv[25]);
    if (argc > 26) params.torus_r_max_rg = std::atof(argv[26]);
    if (argc > 27) params.rho_floor = std::atof(argv[27]);
    if (argc > 28) {
        source_model = argv[28];
        params.source_model_kind = parse_source_model_kind(source_model);
    } else {
        source_model = source_model_name_from_kind(params.source_model_kind);
    }
    if (argc > 29) params.source_funnel_theta_rad = std::atof(argv[29]) * PI / 180.0;
    if (argc > 30) params.source_density_power_q = std::atof(argv[30]);
    if (argc > 31) params.source_radial_power_s = std::atof(argv[31]);
    if (argc > 32) params.source_gradient_dr_rg = std::atof(argv[32]);
    if (argc > 33) params.source_gradient_dtheta_rad = std::atof(argv[33]) * PI / 180.0;
    if (argc > 34) params.source_rho_ref_gcm3 = std::atof(argv[34]);
    if (argc > 35) params.source_cutoff_min = std::atof(argv[35]);
    if (argc > 36) params.source_cutoff_max = std::atof(argv[36]);

    if (chunk_rays <= 0) {
        std::cerr
            << "Usage: " << argv[0]
            << " --image-from-cache ENU ASPIN MBH_MSUN TORUS_RHO0 TORUS_R0"
            << " TORUS_SIGMA TORUS_H_OVER_R SOURCE_R SOURCE_SIGMA SOURCE_THETA_DEG"
            << " SOURCE_POWERLAW SOURCE_EMAX SOURCE_NORM MEV_ENU MEV_NORM"
            << " CAM_THETA_DEG GPU_CACHE_RAYS [PROFILE RADIAL_POWER"
            << " FUNNEL_DEPLETION FUNNEL_THETA_DEG ENVELOPE_RHO0"
            << " ENVELOPE_ALPHA R_MIN R_MAX RHO_FLOOR"
            << " SOURCE_MODEL SOURCE_FUNNEL_THETA_DEG SOURCE_DENSITY_Q"
            << " SOURCE_RADIAL_S SOURCE_GRADIENT_DR SOURCE_GRADIENT_DTHETA_DEG]\n";
        return 1;
    }

    const HostSigmaTable sigma =
        load_sigma_table("data/sigma/sigma_nuN_CC_GBW.dat");
    params.sigma_n = static_cast<int>(sigma.energy_GeV.size());

    std::ifstream in(
        "output/rays/kerr_geodesics_cuda.bin",
        std::ios::binary
    );

    if (!in) {
        throw std::runtime_error("Could not open output/rays/kerr_geodesics_cuda.bin");
    }

    std::int32_t magic = 0;
    std::int32_t version = 0;
    std::int32_t nx = 0;
    std::int32_t ny = 0;
    double cache_spin = 0.0;

    in.read(reinterpret_cast<char*>(&magic), sizeof(magic));
    in.read(reinterpret_cast<char*>(&version), sizeof(version));
    in.read(reinterpret_cast<char*>(&nx), sizeof(nx));
    in.read(reinterpret_cast<char*>(&ny), sizeof(ny));
    in.read(reinterpret_cast<char*>(&cache_spin), sizeof(cache_spin));

    if (magic != 0x4B47454F || version != 1) {
        throw std::runtime_error("Invalid output/rays/kerr_geodesics_cuda.bin");
    }

    params.nx = nx;
    params.ny = ny;

    if (std::abs(cache_spin - params.a_spin) > 1.0e-12) {
        std::cerr
            << "Warning: cache spin a=" << cache_spin
            << " differs from requested ASPIN=" << params.a_spin << "\n";
    }

    double* d_sigma_E = nullptr;
    double* d_sigma_cm2 = nullptr;
    CacheRayDevice* d_rays = nullptr;
    CachePoint* d_points = nullptr;
    OutputPixel* d_outputs = nullptr;

    CUDA_CHECK(cudaMalloc(&d_sigma_E, sigma.energy_GeV.size() * sizeof(double)));
    CUDA_CHECK(cudaMalloc(&d_sigma_cm2, sigma.sigma_cm2.size() * sizeof(double)));
    CUDA_CHECK(cudaMemcpy(
        d_sigma_E,
        sigma.energy_GeV.data(),
        sigma.energy_GeV.size() * sizeof(double),
        cudaMemcpyHostToDevice
    ));
    CUDA_CHECK(cudaMemcpy(
        d_sigma_cm2,
        sigma.sigma_cm2.data(),
        sigma.sigma_cm2.size() * sizeof(double),
        cudaMemcpyHostToDevice
    ));
    CUDA_CHECK(cudaMalloc(&d_rays, chunk_rays * sizeof(CacheRayDevice)));
    CUDA_CHECK(cudaMalloc(&d_outputs, chunk_rays * sizeof(OutputPixel)));

    std::vector<OutputPixel> image(static_cast<std::size_t>(nx * ny));
    std::vector<CacheRayDevice> rays;
    std::vector<CachePoint> points;
    std::vector<OutputPixel> outputs(static_cast<std::size_t>(chunk_rays));

    const int threads = 128;
    std::size_t total_loaded = 0;
    std::size_t total_points = 0;
    const auto t_start = std::chrono::steady_clock::now();

    auto process_chunk = [&]() {
        if (rays.empty()) {
            return;
        }

        CachePoint* current_points = nullptr;
        CUDA_CHECK(cudaMalloc(&current_points, points.size() * sizeof(CachePoint)));
        d_points = current_points;

        CUDA_CHECK(cudaMemcpy(
            d_rays,
            rays.data(),
            rays.size() * sizeof(CacheRayDevice),
            cudaMemcpyHostToDevice
        ));
        CUDA_CHECK(cudaMemcpy(
            d_points,
            points.data(),
            points.size() * sizeof(CachePoint),
            cudaMemcpyHostToDevice
        ));

        const int n_rays = static_cast<int>(rays.size());
        const int blocks = (n_rays + threads - 1) / threads;

        radiative_from_cache_kernel<<<blocks, threads>>>(
            params,
            d_sigma_E,
            d_sigma_cm2,
            d_rays,
            d_points,
            d_outputs,
            n_rays
        );
        CUDA_CHECK(cudaGetLastError());
        CUDA_CHECK(cudaDeviceSynchronize());

        CUDA_CHECK(cudaMemcpy(
            outputs.data(),
            d_outputs,
            rays.size() * sizeof(OutputPixel),
            cudaMemcpyDeviceToHost
        ));

        for (std::size_t k = 0; k < rays.size(); ++k) {
            const OutputPixel& out = outputs[k];
            image[static_cast<std::size_t>(out.pixel_j * nx + out.pixel_i)] = out;
        }

        CUDA_CHECK(cudaFree(d_points));
        d_points = nullptr;
        rays.clear();
        points.clear();
    };

    while (true) {
        BinaryRayHeader header{};
        in.read(reinterpret_cast<char*>(&header), sizeof(header));
        if (!in) {
            break;
        }

        if (header.npoints < 0) {
            throw std::runtime_error("Invalid negative point count in CUDA geodesic cache");
        }

        CacheRayDevice ray{};
        ray.pixel_i = header.pixel_i;
        ray.pixel_j = header.pixel_j;
        ray.captured = header.captured;
        ray.npoints = header.npoints;
        ray.offset = static_cast<int>(points.size());
        ray.alpha_rg = header.alpha_rg;
        ray.beta_rg = header.beta_rg;

        points.resize(points.size() + static_cast<std::size_t>(header.npoints));

        if (header.npoints > 0) {
            in.read(
                reinterpret_cast<char*>(&points[ray.offset]),
                static_cast<std::streamsize>(
                    static_cast<std::size_t>(header.npoints) * sizeof(CachePoint)
                )
            );

            if (!in) {
                throw std::runtime_error("Unexpected EOF while reading CUDA geodesic cache");
            }
        }

        rays.push_back(ray);
        total_loaded += 1;
        total_points += static_cast<std::size_t>(header.npoints);

        if (static_cast<int>(rays.size()) >= chunk_rays) {
            process_chunk();
        }
    }

    process_chunk();

    const std::string output_stem =
        std::string("output/images/kerr_image_cuda_cache_rho0_torus_")
        + make_compact_scientific_tag(params.torus_rho0)
        + (density_profile == "gaussian"
            ? std::string("")
            : std::string("_Profile_") + density_profile)
        + (params.source_model_kind == SOURCE_INNER_RING
            ? std::string("")
            : std::string("_Source_") + source_model)
        + "_Enu_"
        + make_tag("%.0e", params.Enu_obs_GeV)
        + "_MeVEnu_"
        + make_tag("%.0e", params.mev_Enu_obs_MeV)
        + "_MeVNorm_"
        + make_tag("%.0e", params.mev_norm)
        + "_CamTheta_"
        + make_tag("%.1f", cam_theta_deg);

    const std::string binary_filename = output_stem + ".bin";
    const std::string text_filename = output_stem + ".dat";

    write_image_binary(binary_filename, image, params, cam_theta_deg);

    std::ofstream out(text_filename);
    if (!out) {
        throw std::runtime_error("Could not open output file: " + text_filename);
    }

    out << "# profile_type " << density_profile << "\n"
        << "# rho0 " << params.torus_rho0 << "\n"
        << "# rho0_gcm3 " << params.torus_rho0 << "\n"
        << "# r0 " << params.torus_r0_rg << "\n"
        << "# r0_rg " << params.torus_r0_rg << "\n"
        << "# sigma_r " << params.torus_sigma_rg << "\n"
        << "# sigma_r_rg " << params.torus_sigma_rg << "\n"
        << "# H_over_R " << params.torus_h_over_r << "\n"
        << "# radial_power " << params.torus_radial_power << "\n"
        << "# funnel_depletion " << params.funnel_depletion << "\n"
        << "# funnel_theta_deg " << funnel_theta_deg << "\n"
        << "# rho_floor " << params.rho_floor << "\n"
        << "# rho_floor_gcm3 " << params.rho_floor << "\n"
        << "# envelope_rho0 " << params.envelope_rho0 << "\n"
        << "# envelope_rho0_gcm3 " << params.envelope_rho0 << "\n"
        << "# envelope_alpha " << params.envelope_alpha << "\n"
        << "# R_min " << params.torus_r_min_rg << "\n"
        << "# R_min_rg " << params.torus_r_min_rg << "\n"
        << "# R_max " << params.torus_r_max_rg << "\n"
        << "# R_max_rg " << params.torus_r_max_rg << "\n"
        << "# spin " << params.a_spin << "\n"
        << "# observer_distance " << params.r_obs_rg << "\n"
        << "# observer_distance_rg " << params.r_obs_rg << "\n"
        << "# observer_inclination " << cam_theta_deg << "\n"
        << "# observer_inclination_deg " << cam_theta_deg << "\n"
        << "# DIS_model GBW\n"
        << "# use_F3 1\n"
        << "# source_model " << source_model << "\n"
        << "# source_r_rg " << params.source_r_center_rg << "\n"
        << "# source_sigma_rg " << params.source_sigma_r_rg << "\n"
        << "# source_theta_deg " << params.source_theta_width_rad * 180.0 / PI << "\n"
        << "# source_powerlaw " << params.source_powerlaw << "\n"
        << "# source_emax_GeV " << params.source_emax_GeV << "\n"
        << "# source_norm " << params.source_norm << "\n"
        << "# source_funnel_theta_deg " << params.source_funnel_theta_rad * 180.0 / PI << "\n"
        << "# source_rho_ref " << params.source_rho_ref_gcm3 << "\n"
        << "# rho_ref " << params.source_rho_ref_gcm3 << "\n"
        << "# source_q " << params.source_density_power_q << "\n"
        << "# source_s " << params.source_radial_power_s << "\n"
        << "# source_density_power_q " << params.source_density_power_q << "\n"
        << "# source_radial_power_s " << params.source_radial_power_s << "\n"
        << "# source_cutoff_min " << params.source_cutoff_min << "\n"
        << "# source_cutoff_max " << params.source_cutoff_max << "\n"
        << "# source_gradient_dr_rg " << params.source_gradient_dr_rg << "\n"
        << "# source_gradient_dtheta_deg " << params.source_gradient_dtheta_rad * 180.0 / PI << "\n";

    out << "# i j alpha beta tau_uhe P_surv_uhe I_obs_uhe captured "
        << "tau_mev P_surv_mev I_obs_mev r_neutrinosphere_rg leakage_mev\n";

    for (const OutputPixel& p : image) {
        out << std::scientific << std::setprecision(8)
            << p.pixel_i << " "
            << p.pixel_j << " "
            << p.alpha_rg << " "
            << p.beta_rg << " "
            << p.tau_uhe << " "
            << p.P_surv_uhe << " "
            << p.I_uhe << " "
            << p.captured << " "
            << p.tau_mev << " "
            << p.P_surv_mev << " "
            << p.I_mev << " "
            << p.r_neutrinosphere_rg << " "
            << p.leakage_mev << "\n";
    }

    const auto t_end = std::chrono::steady_clock::now();
    const double elapsed_s =
        std::chrono::duration<double>(t_end - t_start).count();

    CUDA_CHECK(cudaFree(d_rays));
    CUDA_CHECK(cudaFree(d_outputs));
    CUDA_CHECK(cudaFree(d_sigma_E));
    CUDA_CHECK(cudaFree(d_sigma_cm2));

    std::cout << "Loaded rays: " << total_loaded << "\n";
    std::cout << "Loaded points: " << total_points << "\n";
    std::cout << "Saved: " << binary_filename << "\n";
    std::cout << "Saved: " << text_filename << "\n";
    std::cout << "CUDA cache RT timing: " << elapsed_s << " s\n";

    return 0;
}

void print_usage(const char* argv0)
{
    std::cerr
        << "Usage:\n"
        << "  " << argv0
        << " ENU ASPIN MBH_MSUN TORUS_RHO0 TORUS_R0 TORUS_SIGMA TORUS_H_OVER_R"
        << " SOURCE_R SOURCE_SIGMA SOURCE_THETA_DEG SOURCE_POWERLAW SOURCE_EMAX SOURCE_NORM"
        << " MEV_ENU MEV_NORM CAM_THETA_DEG CAM_R_OBS CAM_FOV CAM_NX CAM_NY CAM_R_MAX CAM_STEP"
        << " GPU_TILE GPU_STEPS_PER_LAUNCH GPU_MAX_STEPS GPU_TOL"
        << " [PROFILE RADIAL_POWER FUNNEL_DEPLETION FUNNEL_THETA_DEG"
        << " ENVELOPE_RHO0 ENVELOPE_ALPHA R_MIN R_MAX RHO_FLOOR"
        << " SOURCE_MODEL SOURCE_FUNNEL_THETA_DEG SOURCE_DENSITY_Q"
        << " SOURCE_RADIAL_S SOURCE_GRADIENT_DR SOURCE_GRADIENT_DTHETA_DEG]\n";
}

} // namespace

int main(int argc, char* argv[])
{
    try {
        if (argc > 1 && std::string(argv[1]) == "--geodesic-cache") {
            return run_geodesic_cache_gpu(argc, argv);
        }

        if (argc > 1 && std::string(argv[1]) == "--image-from-cache") {
            return run_image_from_cache_gpu(argc, argv);
        }

        Params params{};
        params.Enu_obs_GeV = 1.0e9;
        params.a_spin = 0.9;
        params.M_bh_msun = 3.0;
        params.torus_rho0 = 1.0e-2;
        params.torus_r0_rg = 10.0;
        params.torus_sigma_rg = 5.0;
        params.torus_h_over_r = 0.25;
        params.density_profile_kind = 0;
        params.torus_radial_power = 2.0;
        params.funnel_depletion = 0.0;
        params.funnel_theta_rad = 15.0 * PI / 180.0;
        params.envelope_rho0 = 0.0;
        params.envelope_alpha = 2.5;
        params.torus_r_min_rg = 4.0;
        params.torus_r_max_rg = 60.0;
        params.rho_floor = 1.0e-99;
        params.source_r_center_rg = 3.5;
        params.source_sigma_r_rg = 1.0;
        params.source_theta_width_rad = 15.0 * PI / 180.0;
        params.source_powerlaw = 2.0;
        params.source_emax_GeV = 1.0e12;
        params.source_norm = 1.0;
        params.source_model_kind = SOURCE_INNER_RING;
        params.source_funnel_theta_rad = 20.0 * PI / 180.0;
        params.source_density_power_q = 1.0;
        params.source_radial_power_s = 2.0;
        params.source_rho_ref_gcm3 = -1.0;
        params.source_cutoff_min = 0.0;
        params.source_cutoff_max = 1.0e2;
        params.source_gradient_dr_rg = 0.1;
        params.source_gradient_dtheta_rad = 1.0 * PI / 180.0;
        params.mev_Enu_obs_MeV = 10.0;
        params.mev_norm = 1.0;
        params.mev_sigma_abs0_cm2 = 9.6e-44;
        params.mev_sigma_scat0_cm2 = 1.7e-44;
        params.mev_neutrinosphere_tau = 2.0 / 3.0;
        params.r_obs_rg = 60.0;
        params.theta_obs_rad = 80.0 * PI / 180.0;
        params.fov_rad = 25.0 * PI / 180.0;
        params.nx = 100;
        params.ny = 100;
        params.r_max_rg = 120.0;
        params.h = 0.001;
        params.tolerance = 1.0e-8;
        params.max_steps = 200000;

        int gpu_tile = 32;
        int gpu_steps_per_launch = 128;
        std::string density_profile = "gaussian";
        std::string source_model = "inner_ring";
        double funnel_theta_deg = 15.0;

        if (argc > 1) params.Enu_obs_GeV = std::atof(argv[1]);
        if (argc > 2) params.a_spin = std::atof(argv[2]);
        if (argc > 3) params.M_bh_msun = std::atof(argv[3]);
        if (argc > 4) params.torus_rho0 = std::atof(argv[4]);
        if (argc > 5) params.torus_r0_rg = std::atof(argv[5]);
        if (argc > 6) params.torus_sigma_rg = std::atof(argv[6]);
        if (argc > 7) params.torus_h_over_r = std::atof(argv[7]);
        if (argc > 8) params.source_r_center_rg = std::atof(argv[8]);
        if (argc > 9) params.source_sigma_r_rg = std::atof(argv[9]);
        if (argc > 10) params.source_theta_width_rad = std::atof(argv[10]) * PI / 180.0;
        if (argc > 11) params.source_powerlaw = std::atof(argv[11]);
        if (argc > 12) params.source_emax_GeV = std::atof(argv[12]);
        if (argc > 13) params.source_norm = std::atof(argv[13]);
        if (argc > 14) params.mev_Enu_obs_MeV = std::atof(argv[14]);
        if (argc > 15) params.mev_norm = std::atof(argv[15]);
        if (argc > 16) params.theta_obs_rad = std::atof(argv[16]) * PI / 180.0;
        if (argc > 17) params.r_obs_rg = std::atof(argv[17]);
        if (argc > 18) params.fov_rad = std::atof(argv[18]) * PI / 180.0;
        if (argc > 19) params.nx = std::atoi(argv[19]);
        if (argc > 20) params.ny = std::atoi(argv[20]);
        if (argc > 21) params.r_max_rg = std::atof(argv[21]);
        if (argc > 22) params.h = std::atof(argv[22]);
        if (argc > 23) gpu_tile = std::atoi(argv[23]);
        if (argc > 24) gpu_steps_per_launch = std::atoi(argv[24]);
        if (argc > 25) params.max_steps = std::atoi(argv[25]);
        if (argc > 26) params.tolerance = std::atof(argv[26]);
        if (argc > 27) {
            density_profile = argv[27];
            params.density_profile_kind =
                parse_density_profile_kind(density_profile);
        }
        if (argc > 28) params.torus_radial_power = std::atof(argv[28]);
        if (argc > 29) params.funnel_depletion = std::atof(argv[29]);
        if (argc > 30) {
            funnel_theta_deg = std::atof(argv[30]);
            params.funnel_theta_rad = funnel_theta_deg * PI / 180.0;
        }
        if (argc > 31) params.envelope_rho0 = std::atof(argv[31]);
        if (argc > 32) params.envelope_alpha = std::atof(argv[32]);
        if (argc > 33) params.torus_r_min_rg = std::atof(argv[33]);
        if (argc > 34) params.torus_r_max_rg = std::atof(argv[34]);
        if (argc > 35) params.rho_floor = std::atof(argv[35]);
        if (argc > 36) {
            source_model = argv[36];
            params.source_model_kind = parse_source_model_kind(source_model);
        } else {
            source_model = source_model_name_from_kind(params.source_model_kind);
        }
        if (argc > 37) params.source_funnel_theta_rad = std::atof(argv[37]) * PI / 180.0;
        if (argc > 38) params.source_density_power_q = std::atof(argv[38]);
        if (argc > 39) params.source_radial_power_s = std::atof(argv[39]);
        if (argc > 40) params.source_gradient_dr_rg = std::atof(argv[40]);
        if (argc > 41) params.source_gradient_dtheta_rad = std::atof(argv[41]) * PI / 180.0;
        if (argc > 42) params.source_rho_ref_gcm3 = std::atof(argv[42]);
        if (argc > 43) params.source_cutoff_min = std::atof(argv[43]);
        if (argc > 44) params.source_cutoff_max = std::atof(argv[44]);

        if (params.nx <= 0 || params.ny <= 0 || gpu_tile <= 0 ||
            gpu_steps_per_launch <= 0 || params.max_steps <= 0) {
            print_usage(argv[0]);
            return 1;
        }

        const HostSigmaTable sigma =
            load_sigma_table("data/sigma/sigma_nuN_CC_GBW.dat");
        params.sigma_n = static_cast<int>(sigma.energy_GeV.size());

        int device_count = 0;
        CUDA_CHECK(cudaGetDeviceCount(&device_count));
        if (device_count <= 0) {
            throw std::runtime_error("No CUDA device found.");
        }

        cudaDeviceProp prop{};
        CUDA_CHECK(cudaGetDeviceProperties(&prop, 0));
        std::cout << "CUDA device 0: " << prop.name << "\n";
        std::cout << "Grid: " << params.nx << " x " << params.ny << "\n";
        std::cout << "Tile: " << gpu_tile << " x " << gpu_tile
                  << ", steps/launch: " << gpu_steps_per_launch
                  << ", max steps/ray: " << params.max_steps
                  << ", RKF45 tol: " << params.tolerance << "\n";

        double* d_sigma_E = nullptr;
        double* d_sigma_cm2 = nullptr;
        CUDA_CHECK(cudaMalloc(&d_sigma_E, sigma.energy_GeV.size() * sizeof(double)));
        CUDA_CHECK(cudaMalloc(&d_sigma_cm2, sigma.sigma_cm2.size() * sizeof(double)));
        CUDA_CHECK(cudaMemcpy(
            d_sigma_E,
            sigma.energy_GeV.data(),
            sigma.energy_GeV.size() * sizeof(double),
            cudaMemcpyHostToDevice
        ));
        CUDA_CHECK(cudaMemcpy(
            d_sigma_cm2,
            sigma.sigma_cm2.data(),
            sigma.sigma_cm2.size() * sizeof(double),
            cudaMemcpyHostToDevice
        ));

        const int max_tile_rays = gpu_tile * gpu_tile;
        RayAccum* d_rays = nullptr;
        int* d_active_count = nullptr;
        CUDA_CHECK(cudaMalloc(&d_rays, max_tile_rays * sizeof(RayAccum)));
        CUDA_CHECK(cudaMalloc(&d_active_count, sizeof(int)));

        std::vector<OutputPixel> image(
            static_cast<std::size_t>(params.nx * params.ny)
        );

        const int threads = 128;
        std::uint64_t total_tiles = 0;
        std::uint64_t total_advance_launches = 0;
        std::uint64_t total_accepted_steps = 0;
        const auto t_start = std::chrono::steady_clock::now();

        for (int j0 = 0; j0 < params.ny; j0 += gpu_tile) {
            for (int i0 = 0; i0 < params.nx; i0 += gpu_tile) {
                total_tiles += 1;
                const int tile_nx = std::min(gpu_tile, params.nx - i0);
                const int tile_ny = std::min(gpu_tile, params.ny - j0);
                const int n_rays = tile_nx * tile_ny;
                const int blocks = (n_rays + threads - 1) / threads;

                CUDA_CHECK(cudaMemset(d_active_count, 0, sizeof(int)));
                init_rays_kernel<<<blocks, threads>>>(
                    params,
                    i0,
                    j0,
                    tile_nx,
                    tile_ny,
                    d_rays,
                    d_active_count
                );
                CUDA_CHECK(cudaGetLastError());
                CUDA_CHECK(cudaDeviceSynchronize());

                int active = n_rays;
                while (active > 0) {
                    total_advance_launches += 1;
                    CUDA_CHECK(cudaMemset(d_active_count, 0, sizeof(int)));
                    advance_rays_kernel<<<blocks, threads>>>(
                        params,
                        d_sigma_E,
                        d_sigma_cm2,
                        d_rays,
                        n_rays,
                        gpu_steps_per_launch,
                        d_active_count
                    );
                    CUDA_CHECK(cudaGetLastError());
                    CUDA_CHECK(cudaDeviceSynchronize());
                    CUDA_CHECK(cudaMemcpy(
                        &active,
                        d_active_count,
                        sizeof(int),
                        cudaMemcpyDeviceToHost
                    ));
                }

                std::vector<RayAccum> tile_rays(
                    static_cast<std::size_t>(n_rays)
                );
                CUDA_CHECK(cudaMemcpy(
                    tile_rays.data(),
                    d_rays,
                    tile_rays.size() * sizeof(RayAccum),
                    cudaMemcpyDeviceToHost
                ));

                write_tile_outputs(tile_rays, image, params.nx);

                for (const RayAccum& ray : tile_rays) {
                    total_accepted_steps +=
                        static_cast<std::uint64_t>(ray.steps_done);
                }

                std::cout << "Finished tile i=[" << i0 << ","
                          << (i0 + tile_nx - 1) << "] j=[" << j0 << ","
                          << (j0 + tile_ny - 1) << "]\n";
            }
        }

        const double cam_theta_deg = params.theta_obs_rad * 180.0 / PI;
        const std::string output_stem =
            std::string("output/images/kerr_image_cuda_rho0_torus_")
            + make_compact_scientific_tag(params.torus_rho0)
            + (density_profile == "gaussian"
                ? std::string("")
                : std::string("_Profile_") + density_profile)
            + (params.source_model_kind == SOURCE_INNER_RING
                ? std::string("")
                : std::string("_Source_") + source_model)
            + "_Enu_"
            + make_tag("%.0e", params.Enu_obs_GeV)
            + "_MeVEnu_"
            + make_tag("%.0e", params.mev_Enu_obs_MeV)
            + "_MeVNorm_"
            + make_tag("%.0e", params.mev_norm)
            + "_CamTheta_"
            + make_tag("%.1f", cam_theta_deg);

        const std::string binary_filename =
            output_stem + ".bin";

        const std::string text_filename =
            output_stem + ".dat";

        write_image_binary(
            binary_filename,
            image,
            params,
            cam_theta_deg
        );

        std::ofstream out(text_filename);
        if (!out) {
            throw std::runtime_error("Could not open output file: " + text_filename);
        }

        out << "# profile_type " << density_profile << "\n"
            << "# rho0 " << params.torus_rho0 << "\n"
            << "# rho0_gcm3 " << params.torus_rho0 << "\n"
            << "# r0 " << params.torus_r0_rg << "\n"
            << "# r0_rg " << params.torus_r0_rg << "\n"
            << "# sigma_r " << params.torus_sigma_rg << "\n"
            << "# sigma_r_rg " << params.torus_sigma_rg << "\n"
            << "# H_over_R " << params.torus_h_over_r << "\n"
            << "# radial_power " << params.torus_radial_power << "\n"
            << "# funnel_depletion " << params.funnel_depletion << "\n"
            << "# funnel_theta_deg " << funnel_theta_deg << "\n"
            << "# rho_floor " << params.rho_floor << "\n"
            << "# rho_floor_gcm3 " << params.rho_floor << "\n"
            << "# envelope_rho0 " << params.envelope_rho0 << "\n"
            << "# envelope_rho0_gcm3 " << params.envelope_rho0 << "\n"
            << "# envelope_alpha " << params.envelope_alpha << "\n"
            << "# R_min " << params.torus_r_min_rg << "\n"
            << "# R_min_rg " << params.torus_r_min_rg << "\n"
            << "# R_max " << params.torus_r_max_rg << "\n"
            << "# R_max_rg " << params.torus_r_max_rg << "\n"
            << "# spin " << params.a_spin << "\n"
            << "# observer_distance " << params.r_obs_rg << "\n"
            << "# observer_distance_rg " << params.r_obs_rg << "\n"
            << "# observer_inclination " << cam_theta_deg << "\n"
            << "# observer_inclination_deg " << cam_theta_deg << "\n"
            << "# DIS_model GBW\n"
            << "# use_F3 1\n"
            << "# source_model " << source_model << "\n"
            << "# source_r_rg " << params.source_r_center_rg << "\n"
            << "# source_sigma_rg " << params.source_sigma_r_rg << "\n"
            << "# source_theta_deg " << params.source_theta_width_rad * 180.0 / PI << "\n"
            << "# source_powerlaw " << params.source_powerlaw << "\n"
            << "# source_emax_GeV " << params.source_emax_GeV << "\n"
            << "# source_norm " << params.source_norm << "\n"
            << "# source_funnel_theta_deg " << params.source_funnel_theta_rad * 180.0 / PI << "\n"
            << "# source_rho_ref " << params.source_rho_ref_gcm3 << "\n"
            << "# rho_ref " << params.source_rho_ref_gcm3 << "\n"
            << "# source_q " << params.source_density_power_q << "\n"
            << "# source_s " << params.source_radial_power_s << "\n"
            << "# source_density_power_q " << params.source_density_power_q << "\n"
            << "# source_radial_power_s " << params.source_radial_power_s << "\n"
            << "# source_cutoff_min " << params.source_cutoff_min << "\n"
            << "# source_cutoff_max " << params.source_cutoff_max << "\n"
            << "# source_gradient_dr_rg " << params.source_gradient_dr_rg << "\n"
            << "# source_gradient_dtheta_deg " << params.source_gradient_dtheta_rad * 180.0 / PI << "\n";

        out << "# i j alpha beta tau_uhe P_surv_uhe I_obs_uhe captured "
            << "tau_mev P_surv_mev I_obs_mev r_neutrinosphere_rg leakage_mev\n";

        for (const OutputPixel& p : image) {
            out << std::scientific << std::setprecision(8)
                << p.pixel_i << " "
                << p.pixel_j << " "
                << p.alpha_rg << " "
                << p.beta_rg << " "
                << p.tau_uhe << " "
                << p.P_surv_uhe << " "
                << p.I_uhe << " "
                << p.captured << " "
                << p.tau_mev << " "
                << p.P_surv_mev << " "
                << p.I_mev << " "
                << p.r_neutrinosphere_rg << " "
                << p.leakage_mev << "\n";
        }

        const auto t_end = std::chrono::steady_clock::now();
        const double elapsed_s =
            std::chrono::duration<double>(t_end - t_start).count();
        const double total_rays =
            static_cast<double>(params.nx) * static_cast<double>(params.ny);

        CUDA_CHECK(cudaFree(d_rays));
        CUDA_CHECK(cudaFree(d_active_count));
        CUDA_CHECK(cudaFree(d_sigma_E));
        CUDA_CHECK(cudaFree(d_sigma_cm2));

        std::cout << "Saved: " << binary_filename << "\n";
        std::cout << "Saved: " << text_filename << "\n";
        std::cout << "CUDA timing: " << elapsed_s << " s, "
                  << (total_rays / std::max(elapsed_s, 1.0e-300))
                  << " rays/s, tiles=" << total_tiles
                  << ", advance launches=" << total_advance_launches
                  << ", accepted steps=" << total_accepted_steps << "\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << e.what() << "\n";
        return 1;
    }
}
