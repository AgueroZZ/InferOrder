#!/bin/bash
set -euo pipefail
export LC_ALL=C LANG=C TZ=UTC
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1 BLIS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
export R_PROFILE=/dev/null R_PROFILE_USER=/dev/null
export R_ENVIRON=/dev/null R_ENVIRON_USER=/dev/null
experiment_dir=experiments/isomap_elbo_screen_v040
export R_LIBS="$PWD/$experiment_dir/library:$PWD/experiments/estimate_intrinsic_m_smooth_v032/library${R_LIBS:+:$R_LIBS}"
experiment_r=${MPCURVE_R_BIN:-/home/ziangzhang/R/r-4.3.3-fixed/bin/R}
mkdir -p "$experiment_dir/library" "$experiment_dir/logs"
if [[ ${1:-} == install ]]; then
  exec "$experiment_r" CMD INSTALL --no-multiarch --no-byte-compile \
    -l "$experiment_dir/library" ../MPCurver
fi
exec "$experiment_r" --vanilla --slave --file="$1" --args "${@:2}"
