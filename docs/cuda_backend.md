# CUDA Kerr backend

This backend computes the Kerr ray tracing and radiative transfer directly on
the GPU. It is intentionally separate from the CPU/cache workflow, so the CPU
implementation remains the reference path.

## Safe default run

```bash
make kerr-image-gpu CAM_NX=100 CAM_NY=100 GPU_TILE=32 GPU_STEPS_PER_LAUNCH=128
```

The CUDA executable writes:

```text
output/images/kerr_image_cuda_*.bin
output/images/kerr_image_cuda_*.dat
```

Both files contain the same columns as the CPU image file. The binary file is
intended for production plots and large grids; the text file is kept for quick
inspection and debugging.

## Parameters that control GPU load

- `GPU_TILE`: image tile width/height in pixels. Smaller values keep each batch
  lighter. Start with `16` or `32` on a desktop GPU that is also driving the
  display.
- `GPU_STEPS_PER_LAUNCH`: maximum adaptive geodesic steps per kernel launch.
  Smaller values reduce the time spent inside one CUDA kernel.
- `GPU_MAX_STEPS`: maximum number of accepted geodesic steps per ray.
- `GPU_TOL`: RKF45 local error tolerance. The default is `1.0e-8`, matching the
  CPU integrator default.

## Suggested workflow

Start with a tiny grid:

```bash
make kerr-image-gpu CAM_NX=8 CAM_NY=8 GPU_TILE=8 GPU_STEPS_PER_LAUNCH=32
```

To compute and plot in one command:

```bash
make kerr-image-gpu-plot CAM_NX=8 CAM_NY=8 GPU_TILE=8 GPU_STEPS_PER_LAUNCH=32
```

Then compare with the CPU/cache path before increasing resolution:

```bash
make kerr-criar-cache CAM_NX=8 CAM_NY=8
```

For publication runs, compare CPU and GPU results pixel by pixel for a small
grid before using the GPU backend for large parameter scans.

## Cached GPU geodesics

For scans over many neutrino energies with the same geometry, generate the
geodesic cache once:

```bash
make kerr-geodesics-to-cache-gpu CAM_NX=100 CAM_NY=100 GPU_CACHE_RAYS=16 GPU_CACHE_MAX_POINTS=50000
```

This writes:

```text
output/rays/kerr_geodesics_cuda.bin
```

Then reuse that cache for different energies:

```bash
make image-from-gpu-cache ENU=1e8 GPU_CACHE_RAYS=64
make image-from-gpu-cache ENU=1e9 GPU_CACHE_RAYS=64
make image-from-gpu-cache ENU=1e10 GPU_CACHE_RAYS=64
```

The cached-image path writes:

```text
output/images/kerr_image_cuda_cache_*.bin
output/images/kerr_image_cuda_cache_*.dat
```

The plotting scripts prefer these binary cache outputs when present.
