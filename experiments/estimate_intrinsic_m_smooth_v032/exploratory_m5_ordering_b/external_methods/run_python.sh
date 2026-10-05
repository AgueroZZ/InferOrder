#!/bin/bash
set -euo pipefail
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export LC_ALL=C LANG=C MPLCONFIGDIR=/tmp/mpcurve-gpy-matplotlib
SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
exec "$SCRIPT_DIR/.venv/bin/python" "$SCRIPT_DIR/run_gpy.py" "$@"
