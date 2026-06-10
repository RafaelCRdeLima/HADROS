# Recovery Log: 2026-06-09

## Freeze Event

The machine froze during simulation/validation work on 2026-06-09. Files being
written at that moment are treated as potentially incomplete or corrupted.

## Quarantined Files

Potentially unsafe outputs were moved to:

```text
output/quarantine_after_freeze/20260609_1252/
```

Quarantined files include:

```text
output/quarantine_after_freeze/20260609_1252/output_rays/kerr_geodesics.bin
output/quarantine_after_freeze/20260609_1252/output_images/kerr_image_cuda_cache_GBW_rho0_torus_1e-2_Enu_1e_12_MeVEnu_1e_01_MeVNorm_1e_00_CamTheta_80p0.dat
output/quarantine_after_freeze/20260609_1252/output_profiles/torus_synthetic_grid.dat
output/quarantine_after_freeze/20260609_1252/validation_plots/validation_density_backgrounds/
```

The quarantined image file listed above had size `0` bytes and must not be used.

## Caches Not To Use

Do not use the quarantined CPU ray cache:

```text
output/quarantine_after_freeze/20260609_1252/output_rays/kerr_geodesics.bin
```

Do not use any quarantined `.dat`, `.bin`, or `.png` file for final figures or
paper diagnostics.

## CUDA Cache Status

The file:

```text
output/rays/kerr_geodesics_cuda.bin
```

was not modified during the 2026-06-09 freeze/recovery session. It was later
moved out of the active cache path before CUDA validation so it could not be
used accidentally.

Update after the CUDA validation attempt: the legacy unvalidated CUDA cache was
preserved at:

```text
output/quarantine_after_freeze/20260609_cuda_validation/kerr_geodesics_cuda_unvalidated_legacy.bin
```

`nvcc` is available in the `dis` micromamba environment. A tiny 8x8 CPU/CUDA
smoke validation was completed after running with host GPU access. A medium
32x32 CPU/CUDA validation with nonzero image intensity was also completed.

Reports:

```text
output/validation/cuda_cpu_comparison.txt
output/validation/cuda_cpu_comparison_medium.txt
```

This validates validation-resolution opacity/survival and image-intensity
consistency for `gaussian` and `powerlaw_funnel_envelope`. It does not validate
the quarantined 89 GB legacy cache or replace a production-resolution CPU/CUDA
comparison.

Required comparisons:

```text
tau_CPU vs tau_CUDA
Psurv_CPU vs Psurv_CUDA
```

## Safe Small Validation

Run a cheap validation without touching production caches:

```bash
make validate_small NTHREADS=4
```

For an even smaller smoke test:

```bash
make validate_small NTHREADS=2 SMALL_NX=16 SMALL_NY=16 SMALL_R_MAX_RG=30 SMALL_CAM_STEP=0.02 SMALL_TORUS_NR=32 SMALL_TORUS_NTH=24
```

This writes a small separate cache:

```text
output/rays/kerr_geodesics_small.bin
```

## Production Runs With Limited Threads

Production must be launched explicitly. The default `make` target does not run
production simulations.

Use limited OpenMP threads:

```bash
make run_production NTHREADS=4
```

If using make parallelism, limit make jobs as well:

```bash
make -j2 run_production NTHREADS=4
```

For custom production parameters, pass them explicitly, for example:

```bash
make -j2 run_production NTHREADS=4 CAM_NX=300 CAM_NY=300 ENU=1e12 TORUS_RHO0=1e10
```

CUDA validation and CUDA production require `nvcc` plus an accessible NVIDIA
driver/device and should be run separately.
