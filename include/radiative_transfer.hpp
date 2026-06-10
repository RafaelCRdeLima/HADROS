#ifndef RADIATIVE_TRANSFER_HPP
#define RADIATIVE_TRANSFER_HPP

#include "ray.hpp"
#include "torus_profile.hpp"
#include "sigma_table.hpp"
#include "mev_neutrino_physics.hpp"

#include <string>

struct RTResult {
    double tau = 0.0;
    double P_surv = 1.0;
    double I_obs = 0.0;

    double tau_mev = 0.0;
    double P_surv_mev = 1.0;
    double I_obs_mev = 0.0;
    double r_neutrinosphere_rg = -1.0;
    double leakage_factor = 1.0;
};

enum class UHESpectralModel {
    Monochromatic,
    PowerLaw,
    PowerLawCutoff
};

struct UHESpectralParams {
    UHESpectralModel model = UHESpectralModel::Monochromatic;
    std::string model_name = "monochromatic";
    double gamma = 2.0;
    double ecut_GeV = 1.0e12;
    double e_min_GeV = 1.0e5;
    double e_max_GeV = 1.0e12;
    int n_bins = 8;
};

enum class UHESourceModel {
    InnerRing,
    FunnelWall,
    JetBase,
    ShockLayer,
    DensityWeighted
};

struct UHESourceParams {
    UHESourceModel model = UHESourceModel::InnerRing;
    std::string model_name = "inner_ring";
    double r_center_rg = 3.5;
    double sigma_r_rg = 1.0;
    double theta_width_rad = 8.0 * 3.141592653589793238462643383279502884 / 180.0;
    double powerlaw = 2.0;
    double emax_GeV = 1.0e12;
    double norm = 1.0;
    double funnel_theta_rad = 20.0 * 3.141592653589793238462643383279502884 / 180.0;
    double density_power_q = 1.0;
    double radial_power_s = 2.0;
    double rho_ref_gcm3 = -1.0;
    double cutoff_min = 0.0;
    double cutoff_max = 1.0e2;
    double gradient_dr_rg = 0.1;
    double gradient_dtheta_rad = 1.0 * 3.141592653589793238462643383279502884 / 180.0;
};

using UHECircularSourceParams = UHESourceParams;

struct MeVThermalParams {
    double Enu_obs_MeV = 10.0;
    double norm = 1.0;
    double sigma_abs0_cm2 = 9.6e-44;
    double sigma_scat0_cm2 = 1.7e-44;
    double neutrinosphere_tau = 2.0 / 3.0;
    mev_neutrino::MeVModel model = mev_neutrino::MeVModel::Physical;
    std::string model_name = "physical";
    mev_neutrino::MeVFlavor flavor = mev_neutrino::MeVFlavor::AntiNuE;
    std::string flavor_name = "anti_nu_e";
    bool include_urca = true;
    bool include_pair = true;
    bool include_brems = true;
    bool include_absorption = true;
    bool include_scattering = true;
    bool use_degeneracy_correction = false;
    bool include_abs_n = true;
    bool include_abs_p = true;
    bool include_scat_n = true;
    bool include_scat_p = true;
    bool include_scat_e = true;
    mev_neutrino::MeVThermalProfile thermal_profile = mev_neutrino::MeVThermalProfile::InnerHotTorus;
    std::string thermal_profile_name = "inner_hot_torus";
    mev_neutrino::MeVYeProfile ye_profile = mev_neutrino::MeVYeProfile::NeutronRichTorus;
    std::string ye_profile_name = "neutron_rich_torus";
    mev_neutrino::MeVSpectralMode spectral_mode = mev_neutrino::MeVSpectralMode::Monochromatic;
    std::string spectral_mode_name = "monochromatic";
    double T0_MeV = 6.0;
    double T_floor_MeV = 0.1;
    double T_power = 0.2;
    double Ye_torus = 0.25;
    double Ye_funnel = 0.55;
    double Ye_envelope = 0.45;
    double Ye_floor = 0.01;
    double Ye_ceil = 0.60;
    double E_min_MeV = 3.0;
    double E_max_MeV = 50.0;
    int n_bins = 8;
};

namespace radiative_transfer {

    RTResult integrate_kerr_ray(
        const RayPath& ray,
        double Enu_obs_GeV,
        double M_bh_msun,
        const TorusProfile& torus,
        const SigmaTable& sigma,
        UHESourceParams source = UHESourceParams{},
        MeVThermalParams mev = MeVThermalParams{},
        UHESpectralParams spectral = UHESpectralParams{}
    );

    RTResult integrate_kerr_ray_spectral(
        const RayPath& ray,
        double Enu_obs_GeV,
        double M_bh_msun,
        const TorusProfile& torus,
        const SigmaTable& sigma,
        UHESourceParams source = UHESourceParams{},
        MeVThermalParams mev = MeVThermalParams{},
        UHESpectralParams spectral = UHESpectralParams{}
    );

    double emissivity_collapsar_ring(
        double r_rg,
        double theta,
        double Enu_local_GeV,
        const UHESourceParams& source,
        const UHESpectralParams& spectral = UHESpectralParams{}
    );

    double emissivity_uhe(
        double r_rg,
        double theta,
        double Enu_local_GeV,
        const TorusProfile& torus,
        const UHESourceParams& source,
        const UHESpectralParams& spectral = UHESpectralParams{}
    );

    UHESpectralModel parse_uhe_spectral_model(const std::string& model_name);

    const char* uhe_spectral_model_name(UHESpectralModel model);

    double uhe_spectral_weight(
        double E_GeV,
        const UHESourceParams& source,
        const UHESpectralParams& spectral
    );

    UHESourceModel parse_uhe_source_model(const std::string& model_name);

    const char* uhe_source_model_name(UHESourceModel model);

    double uhe_source_spatial_weight(
        double r_rg,
        double theta,
        const TorusProfile& torus,
        const UHESourceParams& source
    );

    double emissivity_mev_thermal(
        double rho_gcm3,
        double T_MeV,
        double Ye,
        double Enu_local_MeV,
        const MeVThermalParams& mev
    );

    double opacity_mev_absorption_cm_inv(
        double rho_gcm3,
        double Ye,
        double Enu_local_MeV,
        const MeVThermalParams& mev
    );

    double opacity_mev_scattering_cm_inv(
        double rho_gcm3,
        double Ye,
        double Enu_local_MeV,
        const MeVThermalParams& mev
    );

}

#endif
