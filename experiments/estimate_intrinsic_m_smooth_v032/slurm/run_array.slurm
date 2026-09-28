#!/bin/bash
#SBATCH --job-name=mpcurve-smooth-M
#SBATCH --account=pi-zemmour
#SBATCH --partition=zemmour-hm
#SBATCH --qos=zemmour
#SBATCH --array=0-89%12
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=06:00:00
#SBATCH --output=/project/mstephens/ziangzhang/MPCurve/InferOrder/experiments/estimate_intrinsic_m_smooth_v032/slurm/%x_%A_%a.out
#SBATCH --error=/project/mstephens/ziangzhang/MPCurve/InferOrder/experiments/estimate_intrinsic_m_smooth_v032/slurm/%x_%A_%a.err

set -euo pipefail
export LC_ALL=C LANG=C
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1 BLIS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
export R_PROFILE=/dev/null R_PROFILE_USER=/dev/null
export R_ENVIRON=/dev/null R_ENVIRON_USER=/dev/null
export R_DEFAULT_PACKAGES=utils,grDevices,graphics,stats,methods
STUDY_ROOT=/project/mstephens/ziangzhang/MPCurve/InferOrder
STUDY_R_RUNNER=/project/mstephens/ziangzhang/MPCurve/InferOrder/experiments/estimate_intrinsic_m_smooth_v032/run_r.sh
STUDY_PYTHON=/software/python-anaconda-2022.05-el8-x86_64/bin/python3
cd "$STUDY_ROOT"
STUDY_ID=$("$STUDY_PYTHON" - "$SLURM_ARRAY_TASK_ID" <<'PY_ID'
import csv, sys
from pathlib import Path
rows = list(csv.DictReader(Path("experiments/estimate_intrinsic_m_smooth_v032/active_manifest.csv").open()))
index = int(sys.argv[1])
assert 0 <= index < len(rows) == 90
print(rows[index]["id"])
PY_ID
)
echo "Dataset: $STUDY_ID; job: $SLURM_ARRAY_JOB_ID/$SLURM_ARRAY_TASK_ID; host: $(hostname)"
bash "$STUDY_R_RUNNER" experiments/estimate_intrinsic_m_smooth_v032/run_dataset.R "$STUDY_ID"
