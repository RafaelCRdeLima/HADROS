#!/usr/bin/env bash
# Compute all IIM image files needed for the tau vs energy publication figure.
# Densities: 1e6 1e7 1e8 1e9 1e10 1e11   (skip 1e-2)
# Energies:  1e3 1e4 1e5 1e6 1e7 1e8 1e9 1e10 1e11 1e12
# Params: same as Makefile defaults
# Run sequentially (one at a time) to avoid competing for the 89 GB cache.

set -euo pipefail
cd "$(dirname "$0")"

export OMP_NUM_THREADS=6

EXEC=./compute_kerr_image_from_cache
SIGMA=data/sigma/sigma_nuN_CC_IIM.dat
ASPIN=0.0001
MBH=3.0
TR0=10.0  TSIG=5.0  THR=0.25
SR=3.5   SSIG=1.0  STHETA=15.0  SPOW=2.0  SEMAX=1e12  SNORM=1.0
MEV_E=10.0  MEV_N=1.0  CAM=80.0

DENSITIES=(1e6 1e7 1e8 1e9 1e10 1e11)
ENERGIES=(1e3 1e4 1e5 1e6 1e7 1e8 1e9 1e10 1e11 1e12)

total=0
done_count=0

for rho in "${DENSITIES[@]}"; do
    for enu in "${ENERGIES[@]}"; do
        total=$((total + 1))

        # Build expected output filename to skip if already exists
        # Naming convention from C++ code: IIM_rho0_torus_<compact>_Enu_<tag>_...
        # Use python to build the exact tag (matches the C++ snprintf logic)
        outfile=$(python3 -c "
import re
def sci_tag(v):
    s='%.0e'%v
    s=s.replace('+','_').replace('-','m').replace('.','p')
    return s
def compact(v):
    s='%.0e'%v
    m,e=s.split('e')
    return f'{m}e{int(e)}'
rho=$rho; enu=$enu
mev_e=10.0; mev_n=1.0; cam=80.0
rt=f'Enu_{sci_tag(enu)}_MeVEnu_{sci_tag(mev_e)}_MeVNorm_{sci_tag(mev_n)}_CamTheta_{cam:.1f}'.replace('.','p')
print(f'output/images/kerr_image_cuda_cache_IIM_rho0_torus_{compact(rho)}_{rt}.dat')
")

        if [[ -f "$outfile" ]] && [[ $(wc -l < "$outfile") -ge 10000 ]]; then
            echo "[SKIP] $outfile (already complete)"
            done_count=$((done_count + 1))
            continue
        fi

        echo "[$((done_count+1))/$total] rho=$rho  E=$enu  -> $outfile"
        $EXEC $enu $ASPIN $MBH $rho $TR0 $TSIG $THR $SR $SSIG $STHETA $SPOW $SEMAX $SNORM $MEV_E $MEV_N $CAM $SIGMA
        done_count=$((done_count + 1))
        echo "   done."
    done
done

echo "All $total combinations processed."
