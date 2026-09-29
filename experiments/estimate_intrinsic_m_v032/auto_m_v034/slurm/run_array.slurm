#!/bin/bash
#SBATCH --job-name=mpcurve-auto-M-v034-monotone
#SBATCH --account=pi-zemmour
#SBATCH --partition=zemmour-hm
#SBATCH --qos=zemmour
#SBATCH --array=0-89%12
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=02:00:00
#SBATCH --output=/project/mstephens/ziangzhang/MPCurve/InferOrder/experiments/estimate_intrinsic_m_v032/auto_m_v034/slurm/%x_%A_%a.out
#SBATCH --error=/project/mstephens/ziangzhang/MPCurve/InferOrder/experiments/estimate_intrinsic_m_v032/auto_m_v034/slurm/%x_%A_%a.err

set -euo pipefail
export LC_ALL=C LANG=C
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1 BLIS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
export R_PROFILE=/dev/null R_PROFILE_USER=/dev/null
export R_ENVIRON=/dev/null R_ENVIRON_USER=/dev/null
export R_DEFAULT_PACKAGES=utils,grDevices,graphics,stats,methods

STUDY_ROOT=/project/mstephens/ziangzhang/MPCurve/InferOrder
R_RUNNER=$STUDY_ROOT/experiments/estimate_intrinsic_m_smooth_v032/run_r.sh
TASK_INDEX=$((SLURM_ARRAY_TASK_ID + 1))
cd "$STUDY_ROOT"
echo "Manifest index: $TASK_INDEX; job: $SLURM_ARRAY_JOB_ID/$SLURM_ARRAY_TASK_ID; host: $(hostname)"
bash "$R_RUNNER" \
  experiments/estimate_intrinsic_m_v032/auto_m_v034/run_dataset.R \
  "$TASK_INDEX"
