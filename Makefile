# =========================================================
# HADROS Makefile
# Kerr ray tracing + DIS optical depth + radiative transfer
# =========================================================


# =========================================================
# Compiler and flags
# =========================================================

CXX := g++

CXXFLAGS := -O3 -std=c++17 -Wall -Wextra -Iinclude -MMD -MP -fopenmp

LDFLAGS := -fopenmp

NVCC ?= nvcc

NVCCFLAGS := -O3 -std=c++17

PYTHON ?= micromamba run -n dis python


# =========================================================
# Directories
# =========================================================

SRC_DIR    := src
APP_DIR    := apps
SCRIPT_DIR := scripts
BUILD_DIR  := build
OUTPUT_DIR := output
PLOT_DIR   := plots


# =========================================================
# Runtime parameters
# =========================================================
# Energia do Neutrino
ENU ?= 1e5

# Momento angular
ASPIN ?= 0.0001

# Massa do buraco negro em massas solares
MBH_MSUN ?= 3.0

# Parametros da camera Kerr
CAM_R_OBS_RG ?= 60.0
CAM_THETA_DEG ?= 80.0
CAM_FOV_DEG ?= 25.0
CAM_NX ?= 100
CAM_NY ?= 100
CAM_R_MAX_RG ?= 120.0
CAM_STEP ?= 0.001

# Parametros do toro
# Densidade maxima/central do toro em g/cm^3
TORUS_RHO0 ?= 1.0e-2

# Raio central do toro em unidades de r_g = GM/c^2
TORUS_R0_RG ?= 10.0

# Largura radial gaussiana do toro em unidades de r_g
TORUS_SIGMA_RG ?= 5.0

# Espessura angular vertical do toro, aproximadamente H/R
TORUS_H_OVER_R ?= 0.25

# Familia analitica de background:
# gaussian, powerlaw, gaussian_funnel, powerlaw_funnel,
# gaussian_envelope, powerlaw_envelope, powerlaw_funnel_envelope
DENSITY_PROFILE ?= gaussian

# Expoente radial do disco/toro power-law: rho ~ r^(-TORUS_RADIAL_POWER)
TORUS_RADIAL_POWER ?= 2.0

# Evacuacao polar do funil: 0 sem funil, 1 funil quase vazio
FUNNEL_DEPLETION ?= 0.0

# Abertura angular caracteristica do funil, em graus
FUNNEL_THETA_DEG ?= 15.0

# Envelope externo tipo collapsar, em g/cm^3 no raio TORUS_R0_RG
ENVELOPE_RHO0 ?= 0.0

# Expoente radial do envelope: rho_env ~ r^(-ENVELOPE_ALPHA)
ENVELOPE_ALPHA ?= 2.5

# Dominio radial dos perfis power-law/envelope
TORUS_R_MIN_RG ?= 4.0
TORUS_R_MAX_RG ?= 60.0

# Piso de densidade para evitar regioes artificialmente zeradas no funil
RHO_FLOOR ?= 1.0e-99

# Metadata da tabela DIS usada para gerar as imagens
USE_F3 ?= 1

# Parametros da fonte circular UHE interna do collapsar
# Prescricao fenomenologica da fonte UHE:
# inner_ring, funnel_wall, jet_base, shock_layer, density_weighted
SOURCE_MODEL ?= inner_ring

# Raio central da fonte circular/anel de neutrinos UHE, em r_g.
# Use valor menor que TORUS_R0_RG para deixar a fonte dentro do toro.
SOURCE_R_RG ?= 3.5

# Largura radial gaussiana da fonte circular, em r_g
SOURCE_SIGMA_RG ?= 1.0

# Espessura angular da fonte circular ao redor do plano equatorial, em graus
SOURCE_THETA_DEG ?= 15.0

# Indice espectral da emissividade UHE: j_E proporcional a E^(-SOURCE_POWERLAW)
SOURCE_POWERLAW ?= 2.0

# Energia de corte exponencial do espectro UHE da fonte, em GeV
SOURCE_EMAX_GEV ?= 1.0e12

# Normalizacao arbitraria da emissividade UHE da fonte circular
SOURCE_NORM ?= 1.0

# Abertura da parede do funil para SOURCE_MODEL=funnel_wall
SOURCE_FUNNEL_THETA_DEG ?= 20.0

# Expoente q para SOURCE_MODEL=density_weighted: j_UHE ~ rho^q
SOURCE_DENSITY_Q ?= 1.0

# Expoente s para SOURCE_MODEL=density_weighted: j_UHE ~ r^(-s)
SOURCE_RADIAL_S ?= 2.0

# Referencia e limites numericos para SOURCE_MODEL=density_weighted
# SOURCE_RHO_REF < 0 usa uma referencia automatica; por seguranca, o default
# usa a normalizacao do background para evitar referencia no piso de densidade.
SOURCE_RHO_REF ?= $(TORUS_RHO0)
SOURCE_CUTOFF_MIN ?= 0.0
SOURCE_CUTOFF_MAX ?= 1.0e2

# Passos de diferenca finita para SOURCE_MODEL=shock_layer
SOURCE_GRADIENT_DR_RG ?= 0.1
SOURCE_GRADIENT_DTHETA_DEG ?= 1.0

SOURCE_ARGS := $(SOURCE_MODEL) $(SOURCE_FUNNEL_THETA_DEG) $(SOURCE_DENSITY_Q) $(SOURCE_RADIAL_S) $(SOURCE_GRADIENT_DR_RG) $(SOURCE_GRADIENT_DTHETA_DEG) $(SOURCE_RHO_REF) $(SOURCE_CUTOFF_MIN) $(SOURCE_CUTOFF_MAX)

# Modelo espectral UHE. O default preserva o comportamento antigo de energia
# unica; powerlaw e powerlaw_cutoff ativam integracao logaritmica em energia.
SPECTRAL_MODEL ?= monochromatic
SPECTRAL_GAMMA ?= $(SOURCE_POWERLAW)
SPECTRAL_ECUT_GEV ?= $(SOURCE_EMAX_GEV)
SPECTRAL_E_MIN_GEV ?= 1.0e5
SPECTRAL_E_MAX_GEV ?= 1.0e12
SPECTRAL_N_BINS ?= 8

SPECTRAL_ARGS := $(SPECTRAL_MODEL) $(SPECTRAL_GAMMA) $(SPECTRAL_ECUT_GEV) $(SPECTRAL_E_MIN_GEV) $(SPECTRAL_E_MAX_GEV) $(SPECTRAL_N_BINS)

# Energia observada do canal termico MeV
MEV_ENERGY_MEV ?= 10.0
MEV_ENU ?= $(MEV_ENERGY_MEV)

# Normalizacao arbitraria da emissividade termica MeV
MEV_NORM ?= 1.0

# Modelo fisico MeV de baixa energia. O default usa canais fisicamente
# motivados; MEV_MODEL=toy preserva o emissor fenomenologico antigo.
MEV_MODEL ?= physical
MEV_INCLUDE_URCA ?= 1
MEV_INCLUDE_PAIR ?= 1
MEV_INCLUDE_BREMS ?= 1
MEV_INCLUDE_ABSORPTION ?= 1
MEV_INCLUDE_SCATTERING ?= 1
MEV_FLAVOR ?= anti_nu_e
MEV_VALIDATION_TORUS_RHO0 ?= 1.0e10
MEV_VALIDATION_RHO_FLOOR ?= 1.0e-20
MEV_VALIDATION_FOV_DEG ?= 120
MEV_VALIDATION_R_MAX_RG ?= 120
MEV_VALIDATION_CAM_STEP ?= 0.02
MEV_THERMAL_PROFILE ?= inner_hot_torus
MEV_YE_PROFILE ?= neutron_rich_torus
MEV_T0_MEV ?= 6.0
MEV_T_FLOOR_MEV ?= 0.1
MEV_T_POWER ?= 0.2
MEV_YE_TORUS ?= 0.25
MEV_YE_FUNNEL ?= 0.55
MEV_YE_ENVELOPE ?= 0.45
MEV_YE_FLOOR ?= 0.01
MEV_YE_CEIL ?= 0.60
MEV_SPECTRAL_MODE ?= monochromatic
MEV_E_MIN_MEV ?= 3.0
MEV_E_MAX_MEV ?= 50.0
MEV_N_BINS ?= 8
MEV_USE_DEGENERACY_CORRECTION ?= 0
MEV_INCLUDE_ABS_N ?= 1
MEV_INCLUDE_ABS_P ?= 1
MEV_INCLUDE_SCAT_N ?= 1
MEV_INCLUDE_SCAT_P ?= 1
MEV_INCLUDE_SCAT_E ?= 1
MEV_LUMINOSITY_NR ?= 96
MEV_LUMINOSITY_NTH ?= 72
MEV_LUMINOSITY_E_MIN_MEV ?= 1.0
MEV_LUMINOSITY_E_MAX_MEV ?= 80.0
MEV_LUMINOSITY_E_BINS ?= 24

MEV_ARGS := $(MEV_MODEL) $(MEV_FLAVOR) $(MEV_INCLUDE_URCA) $(MEV_INCLUDE_PAIR) $(MEV_INCLUDE_BREMS) $(MEV_INCLUDE_ABSORPTION) $(MEV_INCLUDE_SCATTERING) $(MEV_THERMAL_PROFILE) $(MEV_YE_PROFILE) $(MEV_T0_MEV) $(MEV_T_FLOOR_MEV) $(MEV_T_POWER) $(MEV_YE_TORUS) $(MEV_YE_FUNNEL) $(MEV_YE_ENVELOPE) $(MEV_YE_FLOOR) $(MEV_YE_CEIL) $(MEV_SPECTRAL_MODE) $(MEV_E_MIN_MEV) $(MEV_E_MAX_MEV) $(MEV_N_BINS) $(MEV_USE_DEGENERACY_CORRECTION) $(MEV_INCLUDE_ABS_N) $(MEV_INCLUDE_ABS_P) $(MEV_INCLUDE_SCAT_N) $(MEV_INCLUDE_SCAT_P) $(MEV_INCLUDE_SCAT_E)

# Parametros do backend CUDA.
# Use tiles e poucos passos por lancamento para nao monopolizar a GPU.
GPU_TILE ?= 32
GPU_STEPS_PER_LAUNCH ?= 128
GPU_MAX_STEPS ?= 200000
GPU_TOL ?= 1.0e-8
GPU_CACHE_RAYS ?= 32
GPU_CACHE_MAX_POINTS ?= 50000

# =========================================================
# Numero seguro de processadores para calculos OpenMP
# =========================================================

NTHREADS ?= 4
OMP_NUM_THREADS ?= $(NTHREADS)
export OMP_NUM_THREADS

# Parametros baratos para validacao pequena. Producao deve ser explicita.
SMALL_NX ?= 32
SMALL_NY ?= 32
SMALL_R_MAX_RG ?= 40.0
SMALL_CAM_STEP ?= 0.01
SMALL_TORUS_NR ?= 80
SMALL_TORUS_NTH ?= 48
SMALL_CACHE_PATH ?= $(OUTPUT_DIR)/rays/kerr_geodesics_small.bin
SMALL_CACHE_ENU ?= 1e11
SMALL_CACHE_TORUS_RHO0 ?= 1.0e0

# Parametros medios para validacao CPU/CUDA. Ainda nao sao producao.
MEDIUM_NX ?= 32
MEDIUM_NY ?= 32
MEDIUM_R_MAX_RG ?= 80.0
MEDIUM_CAM_STEP ?= 0.01
MEDIUM_CACHE_PATH ?= $(OUTPUT_DIR)/rays/kerr_geodesics_medium_cuda_cpu.bin

# Parametros de superficies de opacidade UHE
TAU_SURFACE_VALUE ?= 1.0
OPACITY_SURFACE_NTHETA ?= 181
OPACITY_SURFACE_NR ?= 1200

# Parametros da exportacao ParaView. Isto apenas amostra o modelo
# axisimetrico atual em uma grade cartesiana 3D.
PARAVIEW_NX ?= 64
PARAVIEW_NY ?= 64
PARAVIEW_NZ ?= 64
PARAVIEW_BOX_RG ?= 80


# =========================================================
# Source groups
# =========================================================

KERR_SRC := \
	$(SRC_DIR)/kerr_metric.cpp \
	$(SRC_DIR)/kerr_geodesic.cpp \
	$(SRC_DIR)/kerr_camera.cpp

SCHW_SRC := \
	$(SRC_DIR)/schwarzschild_raytracer.cpp

SCHW_CAM_SRC := \
	$(SRC_DIR)/schwarzschild_camera.cpp

COMMON_SRC := \
	$(SRC_DIR)/sigma_table.cpp \
	$(SRC_DIR)/torus_profile.cpp

OPTICAL_SRC := \
	$(SRC_DIR)/optical_depth.cpp

RT_SRC := \
	$(SRC_DIR)/radiative_transfer.cpp \
	$(SRC_DIR)/mev_neutrino_physics.cpp


# =========================================================
# Object groups
# =========================================================

KERR_OBJS := \
	$(BUILD_DIR)/kerr_metric.o \
	$(BUILD_DIR)/kerr_geodesic.o \
	$(BUILD_DIR)/kerr_camera.o

SCHW_OBJS := \
	$(BUILD_DIR)/schwarzschild_raytracer.o

SCHW_CAM_OBJS := \
	$(BUILD_DIR)/schwarzschild_camera.o

COMMON_OBJS := \
	$(BUILD_DIR)/sigma_table.o \
	$(BUILD_DIR)/torus_profile.o

OPTICAL_OBJS := \
	$(BUILD_DIR)/optical_depth.o

RT_OBJS := \
	$(BUILD_DIR)/radiative_transfer.o \
	$(BUILD_DIR)/mev_neutrino_physics.o


# =========================================================
# Executables
# =========================================================

EXECS := \
	dump_torus_grid \
	scan_tau_straight \
	dump_schwarzschild_rays \
	dump_camera_rays \
	compute_tau_camera \
	dump_kerr_camera_rays \
	compute_kerr_image \
	compute_kerr_geodesics \
	compute_kerr_image_from_cache \
	compute_kerr_image_cuda \
	extract_opacity_surface \
	export_paraview_fields \
	validate_gaussian_compatibility \
	test_kerr_metric


# =========================================================
# Phony targets
# =========================================================

.PHONY: all help build validate_small validate_small_cuda validate_medium_cuda validate_density_plots validate_source_plots validate_spectra validate_dis_models validate_mev_physics validate_mev_degeneracy validate_mev_opacity_components validate_mev_upgrades audit_mev_luminosity audit_torus_regime audit_collapsar_ndaf_like mev_multiband_image mev_neutrinosphere mev_luminosity dashboard plot_spectrum_observed plot_multiband_image paper_schematics robustness_scans opacity_surfaces paraview_fields image-from-small-cache validate_gaussian_compatibility_run run_production clean clean-build clean-output dirs \
	torus plot-torus \
	tau-straight \
	dump_schwarzschild_rays rays-schwarzschild \
	dump_camera_rays rays-camera rays-on-torus compute_tau_camera tau-camera \
	kerr-rays kerr-rays-3d kerr-rays-on-torus-3d kerr-image \
	compute_kerr_image kerr-geodesics kerr-criar-cache image-from-cache plot-kerr-image \
	compute_kerr_image_cuda kerr-image-gpu kerr-image-gpu-plot \
	kerr-geodesics-to-cache-gpu image-from-gpu-cache \
	plot-energy-diagnostics \
	plot-geometry-schematic \
	tests


# =========================================================
# Default target
# =========================================================

.DEFAULT_GOAL := help

all: help

help:
	@echo "Safe targets:"
	@echo "  make build                         # compile CPU binaries only"
	@echo "  make validate_small NTHREADS=4      # cheap 32x32 validation run"
	@echo "  make validate_small_cuda NTHREADS=2 # tiny CPU/CUDA comparison"
	@echo "  make validate_medium_cuda NTHREADS=2 # medium CPU/CUDA image comparison"
	@echo "  make validate_source_plots          # UHE source morphology diagnostics"
	@echo "  make validate_spectra               # UHE spectral model validation"
	@echo "  make validate_dis_models NTHREADS=2 # GBW/IIM/PDF_reference DIS validation"
	@echo "  make validate_mev_physics           # MeV emissivity/opacity validation"
	@echo "  make validate_mev_degeneracy        # MeV electron degeneracy diagnostics"
	@echo "  make validate_mev_opacity_components # MeV opacity decomposition diagnostics"
	@echo "  make validate_mev_upgrades          # run all next-generation MeV diagnostics"
	@echo "  make mev_multiband_image            # MeV false-color diagnostic bands"
	@echo "  make mev_neutrinosphere             # MeV tau surfaces and diagnostics"
	@echo "  make mev_luminosity                 # diagnostic MeV luminosity proxy"
	@echo "  make audit_mev_luminosity           # audit MeV luminosity budget/units/gaps"
	@echo "  make audit_torus_regime             # classify current torus physical regime"
	@echo "  make audit_collapsar_ndaf_like      # compare fiducial and collapsar/NDAF-like preset"
	@echo "  make dashboard                      # index existing plots/output products"
	@echo "  make plot_spectrum_observed         # emitted vs observed UHE spectrum"
	@echo "  make plot_multiband_image           # observed image split into energy bands"
	@echo "  make paper_schematics               # updated paper Figures 1 and 2"
	@echo "  make robustness_scans NTHREADS=2    # Point-3 robustness scans"
	@echo "  make opacity_surfaces               # Point-4 tau surfaces and opacity maps"
	@echo "  make paraview_fields NTHREADS=4     # export 3D Cartesian VTK fields"
	@echo "  make image-from-small-cache         # reuse small geodesics for a cheap image"
	@echo "  make run_production NTHREADS=4 ...  # explicit full CPU production run"
	@echo "  make torus NTHREADS=4               # generate semi-analytic density grid"
	@echo ""
	@echo "Use make -j2 build for limited make parallelism."
	@echo "CUDA validation requires nvcc plus an accessible NVIDIA driver/device and must be run separately."

build: dirs dump_torus_grid compute_kerr_geodesics compute_kerr_image_from_cache extract_opacity_surface validate_gaussian_compatibility tests


# =========================================================
# Directory setup
# =========================================================

dirs:
	mkdir -p $(BUILD_DIR)
	mkdir -p $(OUTPUT_DIR)/profiles
	mkdir -p $(OUTPUT_DIR)/tau
	mkdir -p $(OUTPUT_DIR)/rays
	mkdir -p $(OUTPUT_DIR)/images
	mkdir -p $(OUTPUT_DIR)/paraview
	mkdir -p $(PLOT_DIR)


# =========================================================
# Generic compilation rules
# =========================================================

$(BUILD_DIR)/%.o: $(SRC_DIR)/%.cpp | dirs
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(BUILD_DIR)/%.o: $(APP_DIR)/%.cpp | dirs
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(BUILD_DIR)/%.cuda.o: $(APP_DIR)/%.cu | dirs
	$(NVCC) $(NVCCFLAGS) -c $< -o $@


# =========================================================
# Torus profile generation
# =========================================================

dump_torus_grid: \
	$(BUILD_DIR)/torus_profile.o \
	$(BUILD_DIR)/dump_torus_grid.o
	$(CXX) $(LDFLAGS) $^ -o $@

torus: dirs dump_torus_grid
	OMP_NUM_THREADS=$(NTHREADS) ./dump_torus_grid $(TORUS_RHO0) $(TORUS_R0_RG) $(TORUS_SIGMA_RG) $(TORUS_H_OVER_R) $(DENSITY_PROFILE) $(TORUS_RADIAL_POWER) $(FUNNEL_DEPLETION) $(FUNNEL_THETA_DEG) $(ENVELOPE_RHO0) $(ENVELOPE_ALPHA) $(TORUS_R_MIN_RG) $(TORUS_R_MAX_RG) $(RHO_FLOOR) $(SMALL_TORUS_NR) $(SMALL_TORUS_NTH)

plot-torus: torus
	$(PYTHON) $(SCRIPT_DIR)/plot_torus.py

plot-geometry-schematic: dirs
	$(PYTHON) $(SCRIPT_DIR)/plot_geometry_schematic.py $(ASPIN) $(MBH_MSUN) $(TORUS_R0_RG) $(TORUS_SIGMA_RG) $(TORUS_H_OVER_R) $(SOURCE_R_RG) $(SOURCE_SIGMA_RG) $(SOURCE_THETA_DEG) $(CAM_THETA_DEG)

paper_schematics: dirs
	$(PYTHON) $(SCRIPT_DIR)/plot_paper_schematics_point123.py


# =========================================================
# Straight-ray DIS optical depth
# =========================================================

scan_tau_straight: \
	$(COMMON_OBJS) \
	$(OPTICAL_OBJS) \
	$(BUILD_DIR)/scan_tau_straight.o
	$(CXX) $(LDFLAGS) $^ -o $@

tau-straight: dirs scan_tau_straight
	./scan_tau_straight $(ENU)
	$(PYTHON) $(SCRIPT_DIR)/plot_tau_straight.py


# =========================================================
# Schwarzschild flows disabled
# =========================================================

dump_schwarzschild_rays rays-schwarzschild dump_camera_rays rays-camera rays-on-torus compute_tau_camera tau-camera:
	@echo "Fluxo Schwarzschild desativado. Este projeto agora usa apenas Kerr."
	@false


# =========================================================
# Kerr 3D ray tracing
# =========================================================

dump_kerr_camera_rays: \
	$(KERR_OBJS) \
	$(BUILD_DIR)/dump_kerr_camera_rays.o
	$(CXX) $(LDFLAGS) $^ -o $@

kerr-rays: dirs dump_kerr_camera_rays
	./dump_kerr_camera_rays $(ASPIN) $(CAM_R_OBS_RG) $(CAM_THETA_DEG) $(CAM_FOV_DEG) $(CAM_NX) $(CAM_NY) $(CAM_R_MAX_RG) $(CAM_STEP)
	$(PYTHON) $(SCRIPT_DIR)/plot_kerr_camera_rays.py

kerr-rays-3d: dirs dump_kerr_camera_rays
	./dump_kerr_camera_rays $(ASPIN) $(CAM_R_OBS_RG) $(CAM_THETA_DEG) $(CAM_FOV_DEG) $(CAM_NX) $(CAM_NY) $(CAM_R_MAX_RG) $(CAM_STEP)
	$(PYTHON) $(SCRIPT_DIR)/plot_kerr_camera_rays_3d.py

kerr-rays-on-torus-3d: dirs torus
	$(PYTHON) $(SCRIPT_DIR)/plot_kerr_rays_on_torus_3d.py


# =========================================================
# Kerr 3D radiative transfer + DIS
# =========================================================

compute_kerr_image:
	@echo "Fluxo sem cache desativado. Use: make kerr-criar-cache ENU=$(ENU) ASPIN=$(ASPIN)"
	@false

kerr-image: kerr-criar-cache



# =========================================================
# Tests
# =========================================================

test_kerr_metric: \
	$(BUILD_DIR)/kerr_metric.o \
	$(BUILD_DIR)/test_kerr_metric.o
	$(CXX) $(LDFLAGS) $^ -o $@

tests: test_kerr_metric
	./test_kerr_metric


# =========================================================
# Cleaning
# =========================================================

clean-build:
	rm -rf $(BUILD_DIR)
	rm -f $(EXECS)

clean-output:
	rm -rf $(OUTPUT_DIR)/profiles/*
	rm -rf $(OUTPUT_DIR)/tau/*
	rm -rf $(OUTPUT_DIR)/rays/*
	rm -rf $(OUTPUT_DIR)/images/*
	rm -rf $(PLOT_DIR)/*

clean: clean-build clean-output

# =========================================================
# Kerr geodesic cache
# =========================================================

compute_kerr_geodesics: \
	$(KERR_OBJS) \
	$(BUILD_DIR)/compute_kerr_geodesics.o
	$(CXX) $(LDFLAGS) $^ -o $@

kerr-geodesics: dirs compute_kerr_geodesics
	OMP_NUM_THREADS=$(NTHREADS) ./compute_kerr_geodesics $(ASPIN) $(CAM_R_OBS_RG) $(CAM_THETA_DEG) $(CAM_FOV_DEG) $(CAM_NX) $(CAM_NY) $(CAM_R_MAX_RG) $(CAM_STEP) $(OUTPUT_DIR)/rays/kerr_geodesics.bin

# =========================================================
# Kerr image from cached geodesics
# =========================================================

compute_kerr_image_from_cache: \
	$(COMMON_OBJS) \
	$(RT_OBJS) \
	$(BUILD_DIR)/compute_kerr_image_from_cache.o
	$(CXX) $(LDFLAGS) $^ -o $@

extract_opacity_surface: \
	$(COMMON_OBJS) \
	$(BUILD_DIR)/extract_opacity_surface.o
	$(CXX) $(LDFLAGS) $^ -o $@

export_paraview_fields: \
	$(COMMON_OBJS) \
	$(RT_OBJS) \
	$(BUILD_DIR)/export_paraview_fields.o
	$(CXX) $(LDFLAGS) $^ -o $@

validate_gaussian_compatibility: \
	$(COMMON_OBJS) \
	$(RT_OBJS) \
	$(BUILD_DIR)/validate_gaussian_compatibility.o
	$(CXX) $(LDFLAGS) $^ -o $@

validate_gaussian_compatibility_run: dirs validate_gaussian_compatibility
	mkdir -p $(PLOT_DIR)/validation_density_backgrounds
	OMP_NUM_THREADS=$(NTHREADS) ./validate_gaussian_compatibility $(SMALL_CACHE_PATH) $(PLOT_DIR)/validation_density_backgrounds/gaussian_backward_compatibility.txt

validate_density_plots: dirs
	$(PYTHON) $(SCRIPT_DIR)/validate_density_backgrounds.py

validate_source_plots: dirs
	$(PYTHON) $(SCRIPT_DIR)/validate_uhe_source_morphologies.py

validate_spectra: dirs compute_kerr_image_from_cache
	OMP_NUM_THREADS=$(NTHREADS) $(PYTHON) $(SCRIPT_DIR)/validate_uhe_spectral_models.py

validate_mev_physics: dirs compute_kerr_geodesics compute_kerr_image_from_cache
	mkdir -p $(OUTPUT_DIR)/validation
	mkdir -p $(PLOT_DIR)/mev_physics
	OMP_NUM_THREADS=$(NTHREADS) ./compute_kerr_geodesics $(ASPIN) $(CAM_R_OBS_RG) $(CAM_THETA_DEG) $(MEV_VALIDATION_FOV_DEG) $(SMALL_NX) $(SMALL_NY) $(MEV_VALIDATION_R_MAX_RG) $(MEV_VALIDATION_CAM_STEP) $(SMALL_CACHE_PATH)
	OMP_NUM_THREADS=$(NTHREADS) SMALL_CACHE_PATH=$(SMALL_CACHE_PATH) SMALL_CACHE_ENU=$(SMALL_CACHE_ENU) SMALL_CACHE_TORUS_RHO0=$(MEV_VALIDATION_TORUS_RHO0) RHO_FLOOR=$(MEV_VALIDATION_RHO_FLOOR) MEV_ENU=$(MEV_ENU) MEV_ENERGY_MEV=$(MEV_ENERGY_MEV) MEV_NORM=$(MEV_NORM) MEV_FLAVOR=$(MEV_FLAVOR) MEV_INCLUDE_URCA=$(MEV_INCLUDE_URCA) MEV_INCLUDE_PAIR=$(MEV_INCLUDE_PAIR) MEV_INCLUDE_BREMS=$(MEV_INCLUDE_BREMS) MEV_INCLUDE_ABSORPTION=$(MEV_INCLUDE_ABSORPTION) MEV_INCLUDE_SCATTERING=$(MEV_INCLUDE_SCATTERING) MEV_THERMAL_PROFILE=$(MEV_THERMAL_PROFILE) MEV_YE_PROFILE=$(MEV_YE_PROFILE) MEV_T0_MEV=$(MEV_T0_MEV) MEV_T_FLOOR_MEV=$(MEV_T_FLOOR_MEV) MEV_T_POWER=$(MEV_T_POWER) MEV_YE_TORUS=$(MEV_YE_TORUS) MEV_YE_FUNNEL=$(MEV_YE_FUNNEL) MEV_YE_ENVELOPE=$(MEV_YE_ENVELOPE) MEV_YE_FLOOR=$(MEV_YE_FLOOR) MEV_YE_CEIL=$(MEV_YE_CEIL) MEV_SPECTRAL_MODE=$(MEV_SPECTRAL_MODE) MEV_E_MIN_MEV=$(MEV_E_MIN_MEV) MEV_E_MAX_MEV=$(MEV_E_MAX_MEV) MEV_N_BINS=$(MEV_N_BINS) MEV_USE_DEGENERACY_CORRECTION=$(MEV_USE_DEGENERACY_CORRECTION) MEV_INCLUDE_ABS_N=$(MEV_INCLUDE_ABS_N) MEV_INCLUDE_ABS_P=$(MEV_INCLUDE_ABS_P) MEV_INCLUDE_SCAT_N=$(MEV_INCLUDE_SCAT_N) MEV_INCLUDE_SCAT_P=$(MEV_INCLUDE_SCAT_P) MEV_INCLUDE_SCAT_E=$(MEV_INCLUDE_SCAT_E) DENSITY_PROFILE=$(DENSITY_PROFILE) TORUS_RHO0=$(TORUS_RHO0) TORUS_R0_RG=$(TORUS_R0_RG) TORUS_SIGMA_RG=$(TORUS_SIGMA_RG) TORUS_H_OVER_R=$(TORUS_H_OVER_R) SOURCE_MODEL=$(SOURCE_MODEL) SOURCE_RHO_REF=$(SOURCE_RHO_REF) $(PYTHON) $(SCRIPT_DIR)/validate_mev_physics.py

validate_mev_degeneracy: dirs
	mkdir -p $(OUTPUT_DIR)/validation
	mkdir -p $(PLOT_DIR)/mev_physics
	$(PYTHON) $(SCRIPT_DIR)/validate_mev_degeneracy.py

validate_mev_opacity_components: dirs
	mkdir -p $(OUTPUT_DIR)/validation
	mkdir -p $(PLOT_DIR)/mev_physics
	$(PYTHON) $(SCRIPT_DIR)/validate_mev_opacity_components.py

mev_luminosity: dirs
	mkdir -p $(OUTPUT_DIR)/mev_luminosity
	mkdir -p $(OUTPUT_DIR)/validation
	mkdir -p $(PLOT_DIR)/mev_physics
	OMP_NUM_THREADS=$(NTHREADS) MEV_VALIDATION_TORUS_RHO0=$(MEV_VALIDATION_TORUS_RHO0) MEV_VALIDATION_RHO_FLOOR=$(MEV_VALIDATION_RHO_FLOOR) TORUS_RHO0=$(TORUS_RHO0) TORUS_R0_RG=$(TORUS_R0_RG) TORUS_SIGMA_RG=$(TORUS_SIGMA_RG) TORUS_H_OVER_R=$(TORUS_H_OVER_R) TORUS_R_MIN_RG=$(TORUS_R_MIN_RG) TORUS_R_MAX_RG=$(TORUS_R_MAX_RG) MBH_MSUN=$(MBH_MSUN) MEV_T0_MEV=$(MEV_T0_MEV) MEV_T_FLOOR_MEV=$(MEV_T_FLOOR_MEV) MEV_T_POWER=$(MEV_T_POWER) MEV_YE_TORUS=$(MEV_YE_TORUS) MEV_YE_FUNNEL=$(MEV_YE_FUNNEL) MEV_YE_ENVELOPE=$(MEV_YE_ENVELOPE) MEV_YE_FLOOR=$(MEV_YE_FLOOR) MEV_YE_CEIL=$(MEV_YE_CEIL) MEV_LUMINOSITY_NR=$(MEV_LUMINOSITY_NR) MEV_LUMINOSITY_NTH=$(MEV_LUMINOSITY_NTH) MEV_LUMINOSITY_E_MIN_MEV=$(MEV_LUMINOSITY_E_MIN_MEV) MEV_LUMINOSITY_E_MAX_MEV=$(MEV_LUMINOSITY_E_MAX_MEV) MEV_LUMINOSITY_E_BINS=$(MEV_LUMINOSITY_E_BINS) $(PYTHON) $(SCRIPT_DIR)/compute_mev_luminosity.py

validate_mev_upgrades: validate_mev_degeneracy validate_mev_opacity_components mev_luminosity

audit_mev_luminosity: dirs
	mkdir -p $(OUTPUT_DIR)/mev_luminosity
	mkdir -p $(PLOT_DIR)/mev_physics
	OMP_NUM_THREADS=$(NTHREADS) MEV_VALIDATION_TORUS_RHO0=$(MEV_VALIDATION_TORUS_RHO0) MEV_VALIDATION_RHO_FLOOR=$(MEV_VALIDATION_RHO_FLOOR) TORUS_RHO0=$(TORUS_RHO0) TORUS_R0_RG=$(TORUS_R0_RG) TORUS_SIGMA_RG=$(TORUS_SIGMA_RG) TORUS_H_OVER_R=$(TORUS_H_OVER_R) TORUS_R_MIN_RG=$(TORUS_R_MIN_RG) TORUS_R_MAX_RG=$(TORUS_R_MAX_RG) MBH_MSUN=$(MBH_MSUN) MEV_T0_MEV=$(MEV_T0_MEV) MEV_T_FLOOR_MEV=$(MEV_T_FLOOR_MEV) MEV_T_POWER=$(MEV_T_POWER) MEV_YE_TORUS=$(MEV_YE_TORUS) MEV_YE_FLOOR=$(MEV_YE_FLOOR) MEV_YE_CEIL=$(MEV_YE_CEIL) MEV_LUMINOSITY_NR=$(MEV_LUMINOSITY_NR) MEV_LUMINOSITY_NTH=$(MEV_LUMINOSITY_NTH) MEV_LUMINOSITY_E_MIN_MEV=$(MEV_LUMINOSITY_E_MIN_MEV) MEV_LUMINOSITY_E_MAX_MEV=$(MEV_LUMINOSITY_E_MAX_MEV) MEV_LUMINOSITY_E_BINS=$(MEV_LUMINOSITY_E_BINS) $(PYTHON) $(SCRIPT_DIR)/audit_mev_luminosity.py

audit_torus_regime: dirs
	mkdir -p $(OUTPUT_DIR)/torus_regime
	mkdir -p $(PLOT_DIR)/torus_regime
	OMP_NUM_THREADS=$(NTHREADS) MEV_VALIDATION_TORUS_RHO0=$(MEV_VALIDATION_TORUS_RHO0) MEV_VALIDATION_RHO_FLOOR=$(MEV_VALIDATION_RHO_FLOOR) TORUS_RHO0=$(TORUS_RHO0) TORUS_R0_RG=$(TORUS_R0_RG) TORUS_SIGMA_RG=$(TORUS_SIGMA_RG) TORUS_H_OVER_R=$(TORUS_H_OVER_R) TORUS_R_MIN_RG=$(TORUS_R_MIN_RG) TORUS_R_MAX_RG=$(TORUS_R_MAX_RG) MBH_MSUN=$(MBH_MSUN) MEV_T0_MEV=$(MEV_T0_MEV) MEV_T_FLOOR_MEV=$(MEV_T_FLOOR_MEV) MEV_T_POWER=$(MEV_T_POWER) MEV_YE_TORUS=$(MEV_YE_TORUS) MEV_YE_FLOOR=$(MEV_YE_FLOOR) MEV_YE_CEIL=$(MEV_YE_CEIL) MEV_LUMINOSITY_E_MIN_MEV=$(MEV_LUMINOSITY_E_MIN_MEV) MEV_LUMINOSITY_E_MAX_MEV=$(MEV_LUMINOSITY_E_MAX_MEV) MEV_LUMINOSITY_E_BINS=$(MEV_LUMINOSITY_E_BINS) $(PYTHON) $(SCRIPT_DIR)/audit_torus_regime.py

audit_collapsar_ndaf_like: dirs
	mkdir -p $(OUTPUT_DIR)/collapsar_ndaf_like
	mkdir -p $(PLOT_DIR)/collapsar_ndaf_like
	$(PYTHON) $(SCRIPT_DIR)/run_collapsar_ndaf_like_audit.py

mev_multiband_image: dirs compute_kerr_geodesics compute_kerr_image_from_cache
	mkdir -p $(OUTPUT_DIR)/validation
	mkdir -p $(PLOT_DIR)/mev_physics
	OMP_NUM_THREADS=$(NTHREADS) ./compute_kerr_geodesics $(ASPIN) $(CAM_R_OBS_RG) $(CAM_THETA_DEG) $(MEV_VALIDATION_FOV_DEG) $(SMALL_NX) $(SMALL_NY) $(MEV_VALIDATION_R_MAX_RG) $(MEV_VALIDATION_CAM_STEP) $(SMALL_CACHE_PATH)
	OMP_NUM_THREADS=$(NTHREADS) $(PYTHON) $(SCRIPT_DIR)/make_mev_multiband_image.py

mev_neutrinosphere: dirs
	mkdir -p $(OUTPUT_DIR)/mev_neutrinosphere
	mkdir -p $(PLOT_DIR)/mev_physics
	$(PYTHON) $(SCRIPT_DIR)/run_mev_neutrinosphere.py

dashboard:
	mkdir -p dashboard
	$(PYTHON) $(SCRIPT_DIR)/build_plot_dashboard.py

plot_spectrum_observed: dirs compute_kerr_image_from_cache
	OMP_NUM_THREADS=$(NTHREADS) SPECTRAL_MODEL=$(SPECTRAL_MODEL) SPECTRAL_GAMMA=$(SPECTRAL_GAMMA) SPECTRAL_ECUT_GEV=$(SPECTRAL_ECUT_GEV) SPECTRAL_E_MIN_GEV=$(SPECTRAL_E_MIN_GEV) SPECTRAL_E_MAX_GEV=$(SPECTRAL_E_MAX_GEV) SPECTRAL_N_BINS=$(SPECTRAL_N_BINS) DENSITY_PROFILE=$(DENSITY_PROFILE) SOURCE_MODEL=$(SOURCE_MODEL) SMALL_CACHE_TORUS_RHO0=$(SMALL_CACHE_TORUS_RHO0) FUNNEL_DEPLETION=$(FUNNEL_DEPLETION) FUNNEL_THETA_DEG=$(FUNNEL_THETA_DEG) ENVELOPE_RHO0=$(ENVELOPE_RHO0) $(PYTHON) $(SCRIPT_DIR)/plot_emitted_vs_observed_spectrum.py

plot_multiband_image: dirs compute_kerr_image_from_cache
	OMP_NUM_THREADS=$(NTHREADS) SPECTRAL_MODEL=$(SPECTRAL_MODEL) SPECTRAL_GAMMA=$(SPECTRAL_GAMMA) SPECTRAL_ECUT_GEV=$(SPECTRAL_ECUT_GEV) SPECTRAL_E_MIN_GEV=$(SPECTRAL_E_MIN_GEV) SPECTRAL_E_MAX_GEV=$(SPECTRAL_E_MAX_GEV) SPECTRAL_N_BINS=$(SPECTRAL_N_BINS) DENSITY_PROFILE=$(DENSITY_PROFILE) SOURCE_MODEL=$(SOURCE_MODEL) SMALL_CACHE_TORUS_RHO0=$(SMALL_CACHE_TORUS_RHO0) FUNNEL_DEPLETION=$(FUNNEL_DEPLETION) FUNNEL_THETA_DEG=$(FUNNEL_THETA_DEG) ENVELOPE_RHO0=$(ENVELOPE_RHO0) $(PYTHON) $(SCRIPT_DIR)/plot_multiband_observed_image.py

robustness_scans: dirs build
	OMP_NUM_THREADS=$(NTHREADS) $(PYTHON) $(SCRIPT_DIR)/run_robustness_scans.py

opacity_surfaces: dirs extract_opacity_surface
	$(PYTHON) $(SCRIPT_DIR)/run_opacity_surfaces.py --tau-surface $(TAU_SURFACE_VALUE) --ntheta $(OPACITY_SURFACE_NTHETA) --nr $(OPACITY_SURFACE_NR) --enu $(ENU) --mbh $(MBH_MSUN) --sigma-path data/sigma/sigma_nuN_CC_GBW.dat --profile $(DENSITY_PROFILE) --rho0 $(TORUS_RHO0) --r0 $(TORUS_R0_RG) --sigma-r $(TORUS_SIGMA_RG) --h-over-r $(TORUS_H_OVER_R) --radial-power $(TORUS_RADIAL_POWER) --funnel-depletion $(FUNNEL_DEPLETION) --funnel-theta $(FUNNEL_THETA_DEG) --envelope-rho0 $(ENVELOPE_RHO0) --envelope-alpha $(ENVELOPE_ALPHA) --r-min $(TORUS_R_MIN_RG) --r-max $(TORUS_R_MAX_RG) --rho-floor $(RHO_FLOOR)

paraview_fields: dirs export_paraview_fields
	OMP_NUM_THREADS=$(NTHREADS) ./export_paraview_fields $(PARAVIEW_NX) $(PARAVIEW_NY) $(PARAVIEW_NZ) $(PARAVIEW_BOX_RG) $(ENU) $(TORUS_RHO0) $(TORUS_R0_RG) $(TORUS_SIGMA_RG) $(TORUS_H_OVER_R) $(DENSITY_PROFILE) $(TORUS_RADIAL_POWER) $(FUNNEL_DEPLETION) $(FUNNEL_THETA_DEG) $(ENVELOPE_RHO0) $(ENVELOPE_ALPHA) $(TORUS_R_MIN_RG) $(TORUS_R_MAX_RG) $(RHO_FLOOR) $(SOURCE_R_RG) $(SOURCE_SIGMA_RG) $(SOURCE_THETA_DEG) $(SOURCE_POWERLAW) $(SOURCE_EMAX_GEV) $(SOURCE_NORM) $(SOURCE_ARGS) $(OUTPUT_DIR)/paraview/bh_torus_fields.vtk

kerr-criar-cache: dirs kerr-geodesics compute_kerr_image_from_cache
	OMP_NUM_THREADS=$(NTHREADS) ./compute_kerr_image_from_cache $(ENU) $(ASPIN) $(MBH_MSUN) $(TORUS_RHO0) $(TORUS_R0_RG) $(TORUS_SIGMA_RG) $(TORUS_H_OVER_R) $(SOURCE_R_RG) $(SOURCE_SIGMA_RG) $(SOURCE_THETA_DEG) $(SOURCE_POWERLAW) $(SOURCE_EMAX_GEV) $(SOURCE_NORM) $(MEV_ENU) $(MEV_NORM) $(CAM_THETA_DEG) data/sigma/sigma_nuN_CC_GBW.dat $(DENSITY_PROFILE) $(TORUS_RADIAL_POWER) $(FUNNEL_DEPLETION) $(FUNNEL_THETA_DEG) $(ENVELOPE_RHO0) $(ENVELOPE_ALPHA) $(TORUS_R_MIN_RG) $(TORUS_R_MAX_RG) $(RHO_FLOOR) $(CAM_R_OBS_RG) $(USE_F3) $(OUTPUT_DIR)/rays/kerr_geodesics.bin $(SOURCE_ARGS) $(SPECTRAL_ARGS) $(MEV_ARGS)
	$(PYTHON) $(SCRIPT_DIR)/plot_kerr_image.py $(ENU) $(MEV_ENU) $(MEV_NORM) $(CAM_THETA_DEG) $(TORUS_RHO0)

image-from-cache: dirs compute_kerr_image_from_cache
	OMP_NUM_THREADS=$(NTHREADS) ./compute_kerr_image_from_cache $(ENU) $(ASPIN) $(MBH_MSUN) $(TORUS_RHO0) $(TORUS_R0_RG) $(TORUS_SIGMA_RG) $(TORUS_H_OVER_R) $(SOURCE_R_RG) $(SOURCE_SIGMA_RG) $(SOURCE_THETA_DEG) $(SOURCE_POWERLAW) $(SOURCE_EMAX_GEV) $(SOURCE_NORM) $(MEV_ENU) $(MEV_NORM) $(CAM_THETA_DEG) data/sigma/sigma_nuN_CC_GBW.dat $(DENSITY_PROFILE) $(TORUS_RADIAL_POWER) $(FUNNEL_DEPLETION) $(FUNNEL_THETA_DEG) $(ENVELOPE_RHO0) $(ENVELOPE_ALPHA) $(TORUS_R_MIN_RG) $(TORUS_R_MAX_RG) $(RHO_FLOOR) $(CAM_R_OBS_RG) $(USE_F3) $(OUTPUT_DIR)/rays/kerr_geodesics.bin $(SOURCE_ARGS) $(SPECTRAL_ARGS) $(MEV_ARGS)
	$(PYTHON) $(SCRIPT_DIR)/plot_kerr_image.py $(ENU) $(MEV_ENU) $(MEV_NORM) $(CAM_THETA_DEG) $(TORUS_RHO0)

image-from-small-cache: dirs compute_kerr_image_from_cache
	@echo "Reusing small geodesic cache: $(SMALL_CACHE_PATH)"
	@echo "Small-cache image defaults: SMALL_CACHE_ENU=$(SMALL_CACHE_ENU), SMALL_CACHE_TORUS_RHO0=$(SMALL_CACHE_TORUS_RHO0)"
	OMP_NUM_THREADS=$(NTHREADS) ./compute_kerr_image_from_cache $(SMALL_CACHE_ENU) $(ASPIN) $(MBH_MSUN) $(SMALL_CACHE_TORUS_RHO0) $(TORUS_R0_RG) $(TORUS_SIGMA_RG) $(TORUS_H_OVER_R) $(SOURCE_R_RG) $(SOURCE_SIGMA_RG) $(SOURCE_THETA_DEG) $(SOURCE_POWERLAW) $(SOURCE_EMAX_GEV) $(SOURCE_NORM) $(MEV_ENU) $(MEV_NORM) $(CAM_THETA_DEG) data/sigma/sigma_nuN_CC_GBW.dat $(DENSITY_PROFILE) $(TORUS_RADIAL_POWER) $(FUNNEL_DEPLETION) $(FUNNEL_THETA_DEG) $(ENVELOPE_RHO0) $(ENVELOPE_ALPHA) $(TORUS_R_MIN_RG) $(TORUS_R_MAX_RG) $(RHO_FLOOR) $(CAM_R_OBS_RG) $(USE_F3) $(SMALL_CACHE_PATH) $(SOURCE_ARGS) $(SPECTRAL_ARGS) $(MEV_ARGS)
	DENSITY_PROFILE=$(DENSITY_PROFILE) SOURCE_MODEL=$(SOURCE_MODEL) SPECTRAL_MODEL=$(SPECTRAL_MODEL) $(PYTHON) $(SCRIPT_DIR)/plot_kerr_image.py $(SMALL_CACHE_ENU) $(MEV_ENU) $(MEV_NORM) $(CAM_THETA_DEG) $(SMALL_CACHE_TORUS_RHO0)

plot-kerr-image: dirs
	$(PYTHON) $(SCRIPT_DIR)/plot_kerr_image.py $(ENU) $(MEV_ENU) $(MEV_NORM) $(CAM_THETA_DEG) $(TORUS_RHO0)

plot-energy-diagnostics: dirs
	ASPIN=$(ASPIN) MBH_MSUN=$(MBH_MSUN) TORUS_R0_RG=$(TORUS_R0_RG) TORUS_SIGMA_RG=$(TORUS_SIGMA_RG) TORUS_H_OVER_R=$(TORUS_H_OVER_R) MEV_ENU=$(MEV_ENU) $(PYTHON) $(SCRIPT_DIR)/plot_energy_diagnostics.py $(if $(filter command line,$(origin TORUS_RHO0)),$(TORUS_RHO0),)

validate_small: dirs build
	OMP_NUM_THREADS=$(NTHREADS) ./dump_torus_grid $(TORUS_RHO0) $(TORUS_R0_RG) $(TORUS_SIGMA_RG) $(TORUS_H_OVER_R) $(DENSITY_PROFILE) $(TORUS_RADIAL_POWER) $(FUNNEL_DEPLETION) $(FUNNEL_THETA_DEG) $(ENVELOPE_RHO0) $(ENVELOPE_ALPHA) $(TORUS_R_MIN_RG) $(TORUS_R_MAX_RG) $(RHO_FLOOR) $(SMALL_TORUS_NR) $(SMALL_TORUS_NTH)
	OMP_NUM_THREADS=$(NTHREADS) ./compute_kerr_geodesics $(ASPIN) $(CAM_R_OBS_RG) $(CAM_THETA_DEG) $(CAM_FOV_DEG) $(SMALL_NX) $(SMALL_NY) $(SMALL_R_MAX_RG) $(SMALL_CAM_STEP) $(SMALL_CACHE_PATH)
	OMP_NUM_THREADS=$(NTHREADS) ./compute_kerr_image_from_cache $(ENU) $(ASPIN) $(MBH_MSUN) $(TORUS_RHO0) $(TORUS_R0_RG) $(TORUS_SIGMA_RG) $(TORUS_H_OVER_R) $(SOURCE_R_RG) $(SOURCE_SIGMA_RG) $(SOURCE_THETA_DEG) $(SOURCE_POWERLAW) $(SOURCE_EMAX_GEV) $(SOURCE_NORM) $(MEV_ENU) $(MEV_NORM) $(CAM_THETA_DEG) data/sigma/sigma_nuN_CC_GBW.dat $(DENSITY_PROFILE) $(TORUS_RADIAL_POWER) $(FUNNEL_DEPLETION) $(FUNNEL_THETA_DEG) $(ENVELOPE_RHO0) $(ENVELOPE_ALPHA) $(TORUS_R_MIN_RG) $(TORUS_R_MAX_RG) $(RHO_FLOOR) $(CAM_R_OBS_RG) $(USE_F3) $(SMALL_CACHE_PATH) $(SOURCE_ARGS) $(SPECTRAL_ARGS) $(MEV_ARGS)
	$(MAKE) validate_gaussian_compatibility_run NTHREADS=$(NTHREADS) SMALL_CACHE_PATH=$(SMALL_CACHE_PATH)

validate_dis_models: dirs compute_kerr_geodesics compute_kerr_image_from_cache
	OMP_NUM_THREADS=$(NTHREADS) NTHREADS=$(NTHREADS) $(PYTHON) $(SCRIPT_DIR)/validate_dis_models.py

validate_small_cuda: dirs build compute_kerr_image_cuda
	mkdir -p $(OUTPUT_DIR)/validation
	OMP_NUM_THREADS=$(NTHREADS) ./compute_kerr_geodesics $(ASPIN) $(CAM_R_OBS_RG) $(CAM_THETA_DEG) $(CAM_FOV_DEG) $(SMALL_NX) $(SMALL_NY) $(SMALL_R_MAX_RG) $(SMALL_CAM_STEP) $(OUTPUT_DIR)/rays/kerr_geodesics_small_cuda_cpu.bin
	./compute_kerr_image_cuda --geodesic-cache $(ASPIN) $(CAM_R_OBS_RG) $(CAM_THETA_DEG) $(CAM_FOV_DEG) $(SMALL_NX) $(SMALL_NY) $(SMALL_R_MAX_RG) $(SMALL_CAM_STEP) $(GPU_CACHE_RAYS) $(GPU_CACHE_MAX_POINTS) $(GPU_TOL)
	OMP_NUM_THREADS=$(NTHREADS) ./compute_kerr_image_from_cache $(ENU) $(ASPIN) $(MBH_MSUN) $(TORUS_RHO0) $(TORUS_R0_RG) $(TORUS_SIGMA_RG) $(TORUS_H_OVER_R) $(SOURCE_R_RG) $(SOURCE_SIGMA_RG) $(SOURCE_THETA_DEG) $(SOURCE_POWERLAW) $(SOURCE_EMAX_GEV) $(SOURCE_NORM) $(MEV_ENU) $(MEV_NORM) $(CAM_THETA_DEG) data/sigma/sigma_nuN_CC_GBW.dat gaussian $(TORUS_RADIAL_POWER) 0.0 $(FUNNEL_THETA_DEG) 0.0 $(ENVELOPE_ALPHA) $(TORUS_R_MIN_RG) $(TORUS_R_MAX_RG) $(RHO_FLOOR) $(CAM_R_OBS_RG) $(USE_F3) $(OUTPUT_DIR)/rays/kerr_geodesics_small_cuda_cpu.bin $(SOURCE_ARGS)
	./compute_kerr_image_cuda --image-from-cache $(ENU) $(ASPIN) $(MBH_MSUN) $(TORUS_RHO0) $(TORUS_R0_RG) $(TORUS_SIGMA_RG) $(TORUS_H_OVER_R) $(SOURCE_R_RG) $(SOURCE_SIGMA_RG) $(SOURCE_THETA_DEG) $(SOURCE_POWERLAW) $(SOURCE_EMAX_GEV) $(SOURCE_NORM) $(MEV_ENU) $(MEV_NORM) $(CAM_THETA_DEG) $(GPU_CACHE_RAYS) gaussian $(TORUS_RADIAL_POWER) 0.0 $(FUNNEL_THETA_DEG) 0.0 $(ENVELOPE_ALPHA) $(TORUS_R_MIN_RG) $(TORUS_R_MAX_RG) $(RHO_FLOOR) $(SOURCE_ARGS)
	OMP_NUM_THREADS=$(NTHREADS) ./compute_kerr_image_from_cache $(ENU) $(ASPIN) $(MBH_MSUN) $(TORUS_RHO0) $(TORUS_R0_RG) $(TORUS_SIGMA_RG) $(TORUS_H_OVER_R) $(SOURCE_R_RG) $(SOURCE_SIGMA_RG) $(SOURCE_THETA_DEG) $(SOURCE_POWERLAW) $(SOURCE_EMAX_GEV) $(SOURCE_NORM) $(MEV_ENU) $(MEV_NORM) $(CAM_THETA_DEG) data/sigma/sigma_nuN_CC_GBW.dat powerlaw_funnel_envelope $(TORUS_RADIAL_POWER) 1.0 20.0 1.0e-4 $(ENVELOPE_ALPHA) $(TORUS_R_MIN_RG) $(TORUS_R_MAX_RG) $(RHO_FLOOR) $(CAM_R_OBS_RG) $(USE_F3) $(OUTPUT_DIR)/rays/kerr_geodesics_small_cuda_cpu.bin $(SOURCE_ARGS)
	./compute_kerr_image_cuda --image-from-cache $(ENU) $(ASPIN) $(MBH_MSUN) $(TORUS_RHO0) $(TORUS_R0_RG) $(TORUS_SIGMA_RG) $(TORUS_H_OVER_R) $(SOURCE_R_RG) $(SOURCE_SIGMA_RG) $(SOURCE_THETA_DEG) $(SOURCE_POWERLAW) $(SOURCE_EMAX_GEV) $(SOURCE_NORM) $(MEV_ENU) $(MEV_NORM) $(CAM_THETA_DEG) $(GPU_CACHE_RAYS) powerlaw_funnel_envelope $(TORUS_RADIAL_POWER) 1.0 20.0 1.0e-4 $(ENVELOPE_ALPHA) $(TORUS_R_MIN_RG) $(TORUS_R_MAX_RG) $(RHO_FLOOR) $(SOURCE_ARGS)
	$(PYTHON) $(SCRIPT_DIR)/compare_cuda_cpu_metrics.py $(OUTPUT_DIR)/validation/cuda_cpu_comparison.txt gaussian $(OUTPUT_DIR)/images/kerr_image_cuda_cache_GBW_rho0_torus_1e-2_Enu_1e_05_MeVEnu_1e_01_MeVNorm_1e_00_CamTheta_80p0.dat $(OUTPUT_DIR)/images/kerr_image_cuda_cache_rho0_torus_1e-2_Enu_1e_05_MeVEnu_1e_01_MeVNorm_1e_00_CamTheta_80p0.dat powerlaw_funnel_envelope $(OUTPUT_DIR)/images/kerr_image_cuda_cache_GBW_rho0_torus_1e-2_Profile_powerlaw_funnel_envelope_Enu_1e_05_MeVEnu_1e_01_MeVNorm_1e_00_CamTheta_80p0.dat $(OUTPUT_DIR)/images/kerr_image_cuda_cache_rho0_torus_1e-2_Profile_powerlaw_funnel_envelope_Enu_1e_05_MeVEnu_1e_01_MeVNorm_1e_00_CamTheta_80p0.dat

validate_medium_cuda: dirs build compute_kerr_image_cuda
	mkdir -p $(OUTPUT_DIR)/validation
	OMP_NUM_THREADS=$(NTHREADS) ./compute_kerr_geodesics $(ASPIN) $(CAM_R_OBS_RG) $(CAM_THETA_DEG) $(CAM_FOV_DEG) $(MEDIUM_NX) $(MEDIUM_NY) $(MEDIUM_R_MAX_RG) $(MEDIUM_CAM_STEP) $(MEDIUM_CACHE_PATH)
	./compute_kerr_image_cuda --geodesic-cache $(ASPIN) $(CAM_R_OBS_RG) $(CAM_THETA_DEG) $(CAM_FOV_DEG) $(MEDIUM_NX) $(MEDIUM_NY) $(MEDIUM_R_MAX_RG) $(MEDIUM_CAM_STEP) $(GPU_CACHE_RAYS) $(GPU_CACHE_MAX_POINTS) $(GPU_TOL)
	OMP_NUM_THREADS=$(NTHREADS) ./compute_kerr_image_from_cache $(ENU) $(ASPIN) $(MBH_MSUN) $(TORUS_RHO0) $(TORUS_R0_RG) $(TORUS_SIGMA_RG) $(TORUS_H_OVER_R) $(SOURCE_R_RG) $(SOURCE_SIGMA_RG) $(SOURCE_THETA_DEG) $(SOURCE_POWERLAW) $(SOURCE_EMAX_GEV) $(SOURCE_NORM) $(MEV_ENU) $(MEV_NORM) $(CAM_THETA_DEG) data/sigma/sigma_nuN_CC_GBW.dat gaussian $(TORUS_RADIAL_POWER) 0.0 $(FUNNEL_THETA_DEG) 0.0 $(ENVELOPE_ALPHA) $(TORUS_R_MIN_RG) $(TORUS_R_MAX_RG) $(RHO_FLOOR) $(CAM_R_OBS_RG) $(USE_F3) $(MEDIUM_CACHE_PATH) $(SOURCE_ARGS)
	./compute_kerr_image_cuda --image-from-cache $(ENU) $(ASPIN) $(MBH_MSUN) $(TORUS_RHO0) $(TORUS_R0_RG) $(TORUS_SIGMA_RG) $(TORUS_H_OVER_R) $(SOURCE_R_RG) $(SOURCE_SIGMA_RG) $(SOURCE_THETA_DEG) $(SOURCE_POWERLAW) $(SOURCE_EMAX_GEV) $(SOURCE_NORM) $(MEV_ENU) $(MEV_NORM) $(CAM_THETA_DEG) $(GPU_CACHE_RAYS) gaussian $(TORUS_RADIAL_POWER) 0.0 $(FUNNEL_THETA_DEG) 0.0 $(ENVELOPE_ALPHA) $(TORUS_R_MIN_RG) $(TORUS_R_MAX_RG) $(RHO_FLOOR) $(SOURCE_ARGS)
	OMP_NUM_THREADS=$(NTHREADS) ./compute_kerr_image_from_cache $(ENU) $(ASPIN) $(MBH_MSUN) $(TORUS_RHO0) $(TORUS_R0_RG) $(TORUS_SIGMA_RG) $(TORUS_H_OVER_R) $(SOURCE_R_RG) $(SOURCE_SIGMA_RG) $(SOURCE_THETA_DEG) $(SOURCE_POWERLAW) $(SOURCE_EMAX_GEV) $(SOURCE_NORM) $(MEV_ENU) $(MEV_NORM) $(CAM_THETA_DEG) data/sigma/sigma_nuN_CC_GBW.dat powerlaw_funnel_envelope $(TORUS_RADIAL_POWER) 1.0 20.0 1.0e-4 $(ENVELOPE_ALPHA) $(TORUS_R_MIN_RG) $(TORUS_R_MAX_RG) $(RHO_FLOOR) $(CAM_R_OBS_RG) $(USE_F3) $(MEDIUM_CACHE_PATH) $(SOURCE_ARGS)
	./compute_kerr_image_cuda --image-from-cache $(ENU) $(ASPIN) $(MBH_MSUN) $(TORUS_RHO0) $(TORUS_R0_RG) $(TORUS_SIGMA_RG) $(TORUS_H_OVER_R) $(SOURCE_R_RG) $(SOURCE_SIGMA_RG) $(SOURCE_THETA_DEG) $(SOURCE_POWERLAW) $(SOURCE_EMAX_GEV) $(SOURCE_NORM) $(MEV_ENU) $(MEV_NORM) $(CAM_THETA_DEG) $(GPU_CACHE_RAYS) powerlaw_funnel_envelope $(TORUS_RADIAL_POWER) 1.0 20.0 1.0e-4 $(ENVELOPE_ALPHA) $(TORUS_R_MIN_RG) $(TORUS_R_MAX_RG) $(RHO_FLOOR) $(SOURCE_ARGS)
	$(PYTHON) $(SCRIPT_DIR)/compare_cuda_cpu_metrics.py $(OUTPUT_DIR)/validation/cuda_cpu_comparison_medium.txt gaussian $(OUTPUT_DIR)/images/kerr_image_cuda_cache_GBW_rho0_torus_1e-2_Enu_1e_05_MeVEnu_1e_01_MeVNorm_1e_00_CamTheta_80p0.dat $(OUTPUT_DIR)/images/kerr_image_cuda_cache_rho0_torus_1e-2_Enu_1e_05_MeVEnu_1e_01_MeVNorm_1e_00_CamTheta_80p0.dat powerlaw_funnel_envelope $(OUTPUT_DIR)/images/kerr_image_cuda_cache_GBW_rho0_torus_1e-2_Profile_powerlaw_funnel_envelope_Enu_1e_05_MeVEnu_1e_01_MeVNorm_1e_00_CamTheta_80p0.dat $(OUTPUT_DIR)/images/kerr_image_cuda_cache_rho0_torus_1e-2_Profile_powerlaw_funnel_envelope_Enu_1e_05_MeVEnu_1e_01_MeVNorm_1e_00_CamTheta_80p0.dat

run_production: dirs compute_kerr_geodesics compute_kerr_image_from_cache
	@echo "Running explicit CPU production with NTHREADS=$(NTHREADS). Use make -j2 run_production to limit make jobs."
	$(MAKE) kerr-criar-cache NTHREADS=$(NTHREADS)

# =========================================================
# Kerr image directly on CUDA GPU
# =========================================================

compute_kerr_image_cuda: \
	$(BUILD_DIR)/compute_kerr_image_cuda.cuda.o
	$(NVCC) $(NVCCFLAGS) $^ -o $@

kerr-image-gpu: dirs compute_kerr_image_cuda
	./compute_kerr_image_cuda $(ENU) $(ASPIN) $(MBH_MSUN) $(TORUS_RHO0) $(TORUS_R0_RG) $(TORUS_SIGMA_RG) $(TORUS_H_OVER_R) $(SOURCE_R_RG) $(SOURCE_SIGMA_RG) $(SOURCE_THETA_DEG) $(SOURCE_POWERLAW) $(SOURCE_EMAX_GEV) $(SOURCE_NORM) $(MEV_ENU) $(MEV_NORM) $(CAM_THETA_DEG) $(CAM_R_OBS_RG) $(CAM_FOV_DEG) $(CAM_NX) $(CAM_NY) $(CAM_R_MAX_RG) $(CAM_STEP) $(GPU_TILE) $(GPU_STEPS_PER_LAUNCH) $(GPU_MAX_STEPS) $(GPU_TOL) $(DENSITY_PROFILE) $(TORUS_RADIAL_POWER) $(FUNNEL_DEPLETION) $(FUNNEL_THETA_DEG) $(ENVELOPE_RHO0) $(ENVELOPE_ALPHA) $(TORUS_R_MIN_RG) $(TORUS_R_MAX_RG) $(RHO_FLOOR) $(SOURCE_ARGS)

kerr-image-gpu-plot: kerr-image-gpu
	$(PYTHON) $(SCRIPT_DIR)/plot_kerr_image.py $(ENU) $(MEV_ENU) $(MEV_NORM) $(CAM_THETA_DEG) $(TORUS_RHO0)

kerr-geodesics-to-cache-gpu: dirs compute_kerr_image_cuda
	./compute_kerr_image_cuda --geodesic-cache $(ASPIN) $(CAM_R_OBS_RG) $(CAM_THETA_DEG) $(CAM_FOV_DEG) $(CAM_NX) $(CAM_NY) $(CAM_R_MAX_RG) $(CAM_STEP) $(GPU_CACHE_RAYS) $(GPU_CACHE_MAX_POINTS) $(GPU_TOL)

image-from-gpu-cache: dirs compute_kerr_image_cuda
	./compute_kerr_image_cuda --image-from-cache $(ENU) $(ASPIN) $(MBH_MSUN) $(TORUS_RHO0) $(TORUS_R0_RG) $(TORUS_SIGMA_RG) $(TORUS_H_OVER_R) $(SOURCE_R_RG) $(SOURCE_SIGMA_RG) $(SOURCE_THETA_DEG) $(SOURCE_POWERLAW) $(SOURCE_EMAX_GEV) $(SOURCE_NORM) $(MEV_ENU) $(MEV_NORM) $(CAM_THETA_DEG) $(GPU_CACHE_RAYS) $(DENSITY_PROFILE) $(TORUS_RADIAL_POWER) $(FUNNEL_DEPLETION) $(FUNNEL_THETA_DEG) $(ENVELOPE_RHO0) $(ENVELOPE_ALPHA) $(TORUS_R_MIN_RG) $(TORUS_R_MAX_RG) $(RHO_FLOOR) $(SOURCE_ARGS)

# =========================================================
# Automatic dependency tracking
# =========================================================

-include $(BUILD_DIR)/*.d
