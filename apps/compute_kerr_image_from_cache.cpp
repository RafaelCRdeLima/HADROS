#include "ray.hpp"
#include "sigma_table.hpp"
#include "torus_profile.hpp"
#include "radiative_transfer.hpp"

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <vector>
#include <cstdint>
#include <cmath>
#include <cstdio>

#include <algorithm>

#ifdef _OPENMP
#include <omp.h>
#endif

struct BinaryRayHeader {
    std::int32_t ray_id;
    std::int32_t pixel_i;
    std::int32_t pixel_j;
    std::int32_t captured;
    std::int32_t npoints;
    double alpha_rg;
    double beta_rg;
};

int main(int argc, char* argv[])
{
    double Enu_obs_GeV = 1.0e9;
    double a_spin = 0.9;
    double M_bh_msun = 3.0;
    double torus_rho0 = 1.0e10;
    double torus_r0_rg = 10.0;
    double torus_sigma_rg = 5.0;
    double torus_h_over_r = 0.01;
    std::string density_profile = "gaussian";
    double torus_radial_power = 2.0;
    double funnel_depletion = 0.0;
    double funnel_theta_deg = 15.0;
    double envelope_rho0 = 0.0;
    double envelope_alpha = 2.5;
    double torus_r_min_rg = 4.0;
    double torus_r_max_rg = 60.0;
    double rho_floor = 1.0e-99;
    double observer_distance_rg = 60.0;
    int use_F3 = 1;
    std::string requested_cache_path;
    UHECircularSourceParams source;
    UHESpectralParams spectral;
    MeVThermalParams mev;
    double cam_theta_deg = 0.0;

    if (argc > 1) {
        Enu_obs_GeV = std::atof(argv[1]);
    }
    if (argc > 2) {
        a_spin = std::atof(argv[2]);
    }
    if (argc > 3) {
        M_bh_msun = std::atof(argv[3]);
    }
    if (argc > 4) {
        torus_rho0 = std::atof(argv[4]);
    }
    if (argc > 5) {
        torus_r0_rg = std::atof(argv[5]);
    }
    if (argc > 6) {
        torus_sigma_rg = std::atof(argv[6]);
    }
    if (argc > 7) {
        torus_h_over_r = std::atof(argv[7]);
    }
    if (argc > 8) {
        source.r_center_rg = std::atof(argv[8]);
    }
    if (argc > 9) {
        source.sigma_r_rg = std::atof(argv[9]);
    }
    if (argc > 10) {
        source.theta_width_rad = std::atof(argv[10]) * M_PI / 180.0;
    }
    if (argc > 11) {
        source.powerlaw = std::atof(argv[11]);
    }
    if (argc > 12) {
        source.emax_GeV = std::atof(argv[12]);
    }
    if (argc > 13) {
        source.norm = std::atof(argv[13]);
    }
    if (argc > 14) {
        mev.Enu_obs_MeV = std::atof(argv[14]);
    }
    if (argc > 15) {
        mev.norm = std::atof(argv[15]);
    }
    if (argc > 16) {
        cam_theta_deg = std::atof(argv[16]);
    }

    std::string sigma_path = "data/sigma/sigma_nuN_CC_GBW.dat";
    if (argc > 17) {
        sigma_path = argv[17];
    }
    if (argc > 18) {
        density_profile = argv[18];
    }
    if (argc > 19) {
        torus_radial_power = std::atof(argv[19]);
    }
    if (argc > 20) {
        funnel_depletion = std::atof(argv[20]);
    }
    if (argc > 21) {
        funnel_theta_deg = std::atof(argv[21]);
    }
    if (argc > 22) {
        envelope_rho0 = std::atof(argv[22]);
    }
    if (argc > 23) {
        envelope_alpha = std::atof(argv[23]);
    }
    if (argc > 24) {
        torus_r_min_rg = std::atof(argv[24]);
    }
    if (argc > 25) {
        torus_r_max_rg = std::atof(argv[25]);
    }
    if (argc > 26) {
        rho_floor = std::atof(argv[26]);
    }
    if (argc > 27) {
        observer_distance_rg = std::atof(argv[27]);
    }
    if (argc > 28) {
        use_F3 = std::atoi(argv[28]);
    }
    if (argc > 29) {
        requested_cache_path = argv[29];
    }
    if (argc > 30) {
        source.model_name = argv[30];
        source.model = radiative_transfer::parse_uhe_source_model(source.model_name);
    } else {
        source.model_name =
            radiative_transfer::uhe_source_model_name(source.model);
    }
    if (argc > 31) {
        source.funnel_theta_rad = std::atof(argv[31]) * M_PI / 180.0;
    }
    if (argc > 32) {
        source.density_power_q = std::atof(argv[32]);
    }
    if (argc > 33) {
        source.radial_power_s = std::atof(argv[33]);
    }
    if (argc > 34) {
        source.gradient_dr_rg = std::atof(argv[34]);
    }
    if (argc > 35) {
        source.gradient_dtheta_rad = std::atof(argv[35]) * M_PI / 180.0;
    }
    if (argc > 36) {
        source.rho_ref_gcm3 = std::atof(argv[36]);
    }
    if (argc > 37) {
        source.cutoff_min = std::atof(argv[37]);
    }
    if (argc > 38) {
        source.cutoff_max = std::atof(argv[38]);
    }
    if (argc > 39) {
        spectral.model_name = argv[39];
        spectral.model = radiative_transfer::parse_uhe_spectral_model(
            spectral.model_name
        );
    } else {
        spectral.model_name =
            radiative_transfer::uhe_spectral_model_name(spectral.model);
    }
    if (argc > 40) {
        spectral.gamma = std::atof(argv[40]);
    }
    if (argc > 41) {
        spectral.ecut_GeV = std::atof(argv[41]);
    }
    if (argc > 42) {
        spectral.e_min_GeV = std::atof(argv[42]);
    }
    if (argc > 43) {
        spectral.e_max_GeV = std::atof(argv[43]);
    }
    if (argc > 44) {
        spectral.n_bins = std::atoi(argv[44]);
    }
    if (argc > 45) {
        mev.model_name = argv[45];
        mev.model = mev_neutrino::parse_mev_model(mev.model_name);
    } else {
        mev.model_name = mev_neutrino::mev_model_name(mev.model);
    }
    if (argc > 46) {
        mev.flavor_name = argv[46];
        mev.flavor = mev_neutrino::parse_mev_flavor(mev.flavor_name);
    } else {
        mev.flavor_name = mev_neutrino::mev_flavor_name(mev.flavor);
    }
    if (argc > 47) {
        mev.include_urca = std::atoi(argv[47]) != 0;
    }
    if (argc > 48) {
        mev.include_pair = std::atoi(argv[48]) != 0;
    }
    if (argc > 49) {
        mev.include_brems = std::atoi(argv[49]) != 0;
    }
    if (argc > 50) {
        mev.include_absorption = std::atoi(argv[50]) != 0;
    }
    if (argc > 51) {
        mev.include_scattering = std::atoi(argv[51]) != 0;
    }
    if (argc > 52) {
        mev.thermal_profile_name = argv[52];
        mev.thermal_profile = mev_neutrino::parse_mev_thermal_profile(
            mev.thermal_profile_name
        );
    } else {
        mev.thermal_profile_name =
            mev_neutrino::mev_thermal_profile_name(mev.thermal_profile);
    }
    if (argc > 53) {
        mev.ye_profile_name = argv[53];
        mev.ye_profile = mev_neutrino::parse_mev_ye_profile(
            mev.ye_profile_name
        );
    } else {
        mev.ye_profile_name =
            mev_neutrino::mev_ye_profile_name(mev.ye_profile);
    }
    if (argc > 54) {
        mev.T0_MeV = std::atof(argv[54]);
    }
    if (argc > 55) {
        mev.T_floor_MeV = std::atof(argv[55]);
    }
    if (argc > 56) {
        mev.T_power = std::atof(argv[56]);
    }
    if (argc > 57) {
        mev.Ye_torus = std::atof(argv[57]);
    }
    if (argc > 58) {
        mev.Ye_funnel = std::atof(argv[58]);
    }
    if (argc > 59) {
        mev.Ye_envelope = std::atof(argv[59]);
    }
    if (argc > 60) {
        mev.Ye_floor = std::atof(argv[60]);
    }
    if (argc > 61) {
        mev.Ye_ceil = std::atof(argv[61]);
    }
    if (argc > 62) {
        mev.spectral_mode_name = argv[62];
        mev.spectral_mode = mev_neutrino::parse_mev_spectral_mode(
            mev.spectral_mode_name
        );
    } else {
        mev.spectral_mode_name =
            mev_neutrino::mev_spectral_mode_name(mev.spectral_mode);
    }
    if (argc > 63) {
        mev.E_min_MeV = std::atof(argv[63]);
    }
    if (argc > 64) {
        mev.E_max_MeV = std::atof(argv[64]);
    }
    if (argc > 65) {
        mev.n_bins = std::atoi(argv[65]);
    }
    if (argc > 66) {
        mev.use_degeneracy_correction = std::atoi(argv[66]) != 0;
    }
    if (argc > 67) {
        mev.include_abs_n = std::atoi(argv[67]) != 0;
    }
    if (argc > 68) {
        mev.include_abs_p = std::atoi(argv[68]) != 0;
    }
    if (argc > 69) {
        mev.include_scat_n = std::atoi(argv[69]) != 0;
    }
    if (argc > 70) {
        mev.include_scat_p = std::atoi(argv[70]) != 0;
    }
    if (argc > 71) {
        mev.include_scat_e = std::atoi(argv[71]) != 0;
    }

    // Extract short model tag from sigma filename for metadata/output labels.
    std::string sigma_model = "GBW";
    if (
        sigma_path.find("PDF_reference") != std::string::npos
        || sigma_path.find("CTW_reference") != std::string::npos
    ) {
        sigma_model = "PDF_reference";
    } else if (sigma_path.find("IIM") != std::string::npos) {
        sigma_model = "IIM";
    } else if (sigma_path.find("GBW") != std::string::npos) {
        sigma_model = "GBW";
    }

#ifdef _OPENMP
    std::cout
        << "OpenMP max threads = "
        << omp_get_max_threads()
        << "\n";
#else
    std::cout << "OpenMP disabled\n";
#endif

    const std::size_t CHUNK_SIZE = 600;

    SigmaTable sigma(sigma_path);

    TorusProfile torus(
        torus_rho0,
        torus_r0_rg,
        torus_sigma_rg,
        torus_h_over_r,
        density_profile,
        torus_radial_power,
        funnel_depletion,
        funnel_theta_deg * M_PI / 180.0,
        envelope_rho0,
        envelope_alpha,
        torus_r_min_rg,
        torus_r_max_rg,
        rho_floor
    );

    // Use the CUDA geodesic cache (kerr_geodesics_cuda.bin) if it exists,
    // falling back to the CPU cache (kerr_geodesics.bin) otherwise.
    const std::string cuda_cache_path = "output/rays/kerr_geodesics_cuda.bin";
    const std::string cpu_cache_path  = "output/rays/kerr_geodesics.bin";
    const std::string cache_path =
        !requested_cache_path.empty()
        ? requested_cache_path
        : (std::ifstream(cuda_cache_path, std::ios::binary) ? cuda_cache_path : cpu_cache_path);

    std::ifstream in(cache_path, std::ios::binary);

    if (!in) {
        std::cerr << "Could not open geodesic cache (tried "
                  << cuda_cache_path << " and " << cpu_cache_path << ")\n";
        return 1;
    }

    std::cout << "Using geodesic cache: " << cache_path << "\n";
    std::cout << "Using sigma model: " << sigma_model
              << " (" << sigma_path << ")\n";

    std::int32_t magic = 0;
    std::int32_t version = 0;
    std::int32_t nx = 0;
    std::int32_t ny = 0;
    double a_spin_cache = 0.0;

    in.read(reinterpret_cast<char*>(&magic), sizeof(magic));
    in.read(reinterpret_cast<char*>(&version), sizeof(version));
    in.read(reinterpret_cast<char*>(&nx), sizeof(nx));
    in.read(reinterpret_cast<char*>(&ny), sizeof(ny));
    in.read(reinterpret_cast<char*>(&a_spin_cache), sizeof(a_spin_cache));

    const std::int32_t expected_magic = 0x4B47454F; // "KGEO"

    if (magic != expected_magic) {
        std::cerr << "Invalid binary geodesic cache: bad magic number\n";
        return 1;
    }

    if (version != 1) {
        std::cerr << "Unsupported binary geodesic cache version: "
                  << version << "\n";
        return 1;
    }

    if (std::abs(a_spin_cache - a_spin) > 1.0e-12) {
        std::cerr
            << "Warning: cache spin a="
            << a_spin_cache
            << " differs from requested ASPIN="
            << a_spin
            << "\n";
    }

    auto make_tag = [](const char* fmt, double value) {
        char buffer[64];

        std::snprintf(
            buffer,
            sizeof(buffer),
            fmt,
            value
        );

        std::string tag(buffer);

        std::replace(tag.begin(), tag.end(), '+', '_');
        std::replace(tag.begin(), tag.end(), '-', 'm');
        std::replace(tag.begin(), tag.end(), '.', 'p');

        return tag;
    };

    auto make_compact_scientific_tag = [](double value) {
        char buffer[64];

        std::snprintf(
            buffer,
            sizeof(buffer),
            "%.0e",
            value
        );

        std::string tag(buffer);
        const std::size_t epos = tag.find('e');

        if (epos == std::string::npos) {
            return tag;
        }

        const std::string mantissa =
            tag.substr(0, epos);

        const int exponent =
            std::stoi(tag.substr(epos + 1));

        return mantissa
            + "e"
            + std::to_string(exponent);
    };

    const std::string energy_tag =
        make_tag("%.0e", Enu_obs_GeV);

    const std::string mev_energy_tag =
        make_tag("%.0e", mev.Enu_obs_MeV);

    const std::string mev_norm_tag =
        make_tag("%.0e", mev.norm);

    const std::string mev_band_tag =
        mev.spectral_mode == mev_neutrino::MeVSpectralMode::Monochromatic
        ? std::string("")
        : std::string("_MeVSpectrum_")
            + mev.spectral_mode_name
            + "_E"
            + make_tag("%.0e", mev.E_min_MeV)
            + "_"
            + make_tag("%.0e", mev.E_max_MeV);

    const std::string cam_theta_tag =
        make_tag("%.1f", cam_theta_deg);

    std::string image_filename =
        std::string("output/images/kerr_image_cuda_cache_")
        + sigma_model
        + "_rho0_torus_"
        + make_compact_scientific_tag(torus_rho0)
        + (density_profile == "gaussian"
            ? std::string("")
            : std::string("_Profile_") + density_profile)
        + (source.model == UHESourceModel::InnerRing
            ? std::string("")
            : std::string("_Source_") + source.model_name)
        + (spectral.model == UHESpectralModel::Monochromatic
            ? std::string("")
            : std::string("_Spectrum_") + spectral.model_name)
        + "_Enu_"
        + energy_tag
        + "_MeVEnu_"
        + mev_energy_tag
        + "_MeVNorm_"
        + mev_norm_tag
        + mev_band_tag
        + "_CamTheta_"
        + cam_theta_tag
        + ".dat";

    std::ofstream out(image_filename);

    if (!out) {
        std::cerr
            << "Could not open "
            << image_filename
            << "\n";

        return 1;
    }
    out << "# profile_type " << density_profile << "\n"
        << "# rho0 " << torus_rho0 << "\n"
        << "# rho0_gcm3 " << torus_rho0 << "\n"
        << "# r0 " << torus_r0_rg << "\n"
        << "# r0_rg " << torus_r0_rg << "\n"
        << "# sigma_r " << torus_sigma_rg << "\n"
        << "# sigma_r_rg " << torus_sigma_rg << "\n"
        << "# H_over_R " << torus_h_over_r << "\n"
        << "# radial_power " << torus_radial_power << "\n"
        << "# funnel_depletion " << funnel_depletion << "\n"
        << "# funnel_theta_deg " << funnel_theta_deg << "\n"
        << "# rho_floor " << rho_floor << "\n"
        << "# rho_floor_gcm3 " << rho_floor << "\n"
        << "# envelope_rho0 " << envelope_rho0 << "\n"
        << "# envelope_rho0_gcm3 " << envelope_rho0 << "\n"
        << "# envelope_alpha " << envelope_alpha << "\n"
        << "# R_min " << torus_r_min_rg << "\n"
        << "# R_min_rg " << torus_r_min_rg << "\n"
        << "# R_max " << torus_r_max_rg << "\n"
        << "# R_max_rg " << torus_r_max_rg << "\n"
        << "# spin " << a_spin << "\n"
        << "# observer_distance " << observer_distance_rg << "\n"
        << "# observer_distance_rg " << observer_distance_rg << "\n"
        << "# observer_inclination " << cam_theta_deg << "\n"
        << "# observer_inclination_deg " << cam_theta_deg << "\n"
        << "# DIS_model " << sigma_model << "\n"
        << "# use_F3 " << use_F3 << "\n"
        << "# source_model " << source.model_name << "\n"
        << "# source_r_rg " << source.r_center_rg << "\n"
        << "# source_sigma_rg " << source.sigma_r_rg << "\n"
        << "# source_theta_deg " << source.theta_width_rad * 180.0 / M_PI << "\n"
        << "# source_powerlaw " << source.powerlaw << "\n"
        << "# source_emax_GeV " << source.emax_GeV << "\n"
        << "# source_norm " << source.norm << "\n"
        << "# source_funnel_theta_deg " << source.funnel_theta_rad * 180.0 / M_PI << "\n"
        << "# source_rho_ref " << source.rho_ref_gcm3 << "\n"
        << "# rho_ref " << source.rho_ref_gcm3 << "\n"
        << "# source_q " << source.density_power_q << "\n"
        << "# source_s " << source.radial_power_s << "\n"
        << "# source_density_power_q " << source.density_power_q << "\n"
        << "# source_radial_power_s " << source.radial_power_s << "\n"
        << "# source_cutoff_min " << source.cutoff_min << "\n"
        << "# source_cutoff_max " << source.cutoff_max << "\n"
        << "# source_gradient_dr_rg " << source.gradient_dr_rg << "\n"
        << "# source_gradient_dtheta_deg " << source.gradient_dtheta_rad * 180.0 / M_PI << "\n"
        << "# spectral_model " << spectral.model_name << "\n"
        << "# spectral_gamma " << spectral.gamma << "\n"
        << "# spectral_ecut_GeV " << spectral.ecut_GeV << "\n"
        << "# spectral_E_min_GeV " << spectral.e_min_GeV << "\n"
        << "# spectral_E_max_GeV " << spectral.e_max_GeV << "\n"
        << "# spectral_N_bins " << spectral.n_bins << "\n"
        << "# MEV_MODEL " << mev.model_name << "\n"
        << "# MEV_CPU_MODEL " << mev.model_name << "\n"
        << "# MEV_CUDA_STATUS legacy_toy_not_equivalent\n"
        << "# MEV_FLAVOR " << mev.flavor_name << "\n"
        << "# MEV_ENERGY_MEV " << mev.Enu_obs_MeV << "\n"
        << "# MEV_INCLUDE_URCA " << (mev.include_urca ? 1 : 0) << "\n"
        << "# MEV_INCLUDE_PAIR " << (mev.include_pair ? 1 : 0) << "\n"
        << "# MEV_INCLUDE_BREMS " << (mev.include_brems ? 1 : 0) << "\n"
        << "# MEV_INCLUDE_ABSORPTION " << (mev.include_absorption ? 1 : 0) << "\n"
        << "# MEV_INCLUDE_SCATTERING " << (mev.include_scattering ? 1 : 0) << "\n"
        << "# MEV_THERMAL_PROFILE " << mev.thermal_profile_name << "\n"
        << "# MEV_YE_PROFILE " << mev.ye_profile_name << "\n"
        << "# MEV_T0_MEV " << mev.T0_MeV << "\n"
        << "# MEV_T_FLOOR_MEV " << mev.T_floor_MeV << "\n"
        << "# MEV_T_POWER " << mev.T_power << "\n"
        << "# MEV_YE_TORUS " << mev.Ye_torus << "\n"
        << "# MEV_YE_FUNNEL " << mev.Ye_funnel << "\n"
        << "# MEV_YE_ENVELOPE " << mev.Ye_envelope << "\n"
        << "# MEV_YE_FLOOR " << mev.Ye_floor << "\n"
        << "# MEV_YE_CEIL " << mev.Ye_ceil << "\n"
        << "# MEV_SPECTRAL_MODE " << mev.spectral_mode_name << "\n"
        << "# MEV_E_MIN_MEV " << mev.E_min_MeV << "\n"
        << "# MEV_E_MAX_MEV " << mev.E_max_MeV << "\n"
        << "# MEV_N_BINS " << mev.n_bins << "\n"
        << "# MEV_USE_DEGENERACY_CORRECTION " << (mev.use_degeneracy_correction ? 1 : 0) << "\n"
        << "# MEV_INCLUDE_ABS_N " << (mev.include_abs_n ? 1 : 0) << "\n"
        << "# MEV_INCLUDE_ABS_P " << (mev.include_abs_p ? 1 : 0) << "\n"
        << "# MEV_INCLUDE_SCAT_N " << (mev.include_scat_n ? 1 : 0) << "\n"
        << "# MEV_INCLUDE_SCAT_P " << (mev.include_scat_p ? 1 : 0) << "\n"
        << "# MEV_INCLUDE_SCAT_E " << (mev.include_scat_e ? 1 : 0) << "\n"
        << "# mev_emissivity_units local_spectral_emissivity_proxy_per_MeV\n"
        << "# mev_opacity_units cm^-1\n";

    out << "# i j alpha beta tau_uhe P_surv_uhe I_obs_uhe captured "
        << "tau_mev P_surv_mev I_obs_mev r_neutrinosphere_rg leakage_mev\n";

    std::vector<RayPath> chunk;
    chunk.reserve(CHUNK_SIZE);

    std::size_t total_loaded = 0;
    std::size_t total_processed = 0;

    auto process_chunk = [&]() {
        if (chunk.empty()) {
            return;
        }

        std::vector<RTResult> results(chunk.size());

        for (auto& rt : results) {
            rt.tau = 0.0;
            rt.P_surv = 1.0;
            rt.I_obs = 0.0;
        }

#pragma omp parallel for schedule(static)
        for (int k = 0; k < static_cast<int>(chunk.size()); ++k) {
            if (!chunk[k].points.empty()) {
                results[k] =
                    radiative_transfer::integrate_kerr_ray_spectral(
                        chunk[k],
                        Enu_obs_GeV,
                        M_bh_msun,
                        torus,
                        sigma,
                        source,
                        mev,
                        spectral
                    );
            }
        }

        for (std::size_t k = 0; k < chunk.size(); ++k) {
            const RayPath& ray_out = chunk[k];
            const RTResult& rt = results[k];

            out << std::scientific
                << std::setprecision(8)
                << ray_out.pixel_i << " "
                << ray_out.pixel_j << " "
                << ray_out.alpha_rg << " "
                << ray_out.beta_rg << " "
                << rt.tau << " "
                << rt.P_surv << " "
                << rt.I_obs << " "
                << ray_out.captured << " "
                << rt.tau_mev << " "
                << rt.P_surv_mev << " "
                << rt.I_obs_mev << " "
                << rt.r_neutrinosphere_rg << " "
                << rt.leakage_factor << "\n";
        }

        total_processed += chunk.size();
        chunk.clear();
    };

    while (true) {
        BinaryRayHeader header;

        in.read(
            reinterpret_cast<char*>(&header),
            sizeof(header)
        );

        if (!in) {
            break;
        }

        if (header.npoints < 0) {
            std::cerr << "Invalid ray with negative number of points\n";
            return 1;
        }

        RayPath ray;
        ray.a_bh = a_spin_cache;
        ray.pixel_i = header.pixel_i;
        ray.pixel_j = header.pixel_j;
        ray.alpha_rg = header.alpha_rg;
        ray.beta_rg = header.beta_rg;
        ray.captured = (header.captured != 0);

        ray.points.resize(
            static_cast<std::size_t>(header.npoints)
        );

        if (header.npoints > 0) {
            in.read(
                reinterpret_cast<char*>(ray.points.data()),
                static_cast<std::streamsize>(
                    ray.points.size() * sizeof(PathPoint)
                )
            );

            if (!in) {
                std::cerr << "Unexpected end of file while reading points\n";
                return 1;
            }
        }

        chunk.push_back(std::move(ray));
        total_loaded += 1;

        if (chunk.size() >= CHUNK_SIZE) {
            process_chunk();
        }
    }

    process_chunk();

    std::cout << "Cache grid: " << nx << " x " << ny << "\n";
    std::cout << "Loaded rays: " << total_loaded << "\n";
    std::cout << "Processed rays: " << total_processed << "\n";
    std::cout
        << "Saved: "
        << image_filename
        << "\n";

    return 0;
}
