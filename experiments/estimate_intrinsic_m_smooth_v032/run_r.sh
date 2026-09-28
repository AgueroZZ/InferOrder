#!/bin/bash
set -euo pipefail
export LC_ALL=C LANG=C
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1 BLIS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
export R_PROFILE=/dev/null R_PROFILE_USER=/dev/null
export R_ENVIRON=/dev/null R_ENVIRON_USER=/dev/null
STUDY_R=/home/ziangzhang/R/r-4.3.3-fixed/bin/R
exec "$STUDY_R" --vanilla --slave --file="$1" --args "${@:2}"
