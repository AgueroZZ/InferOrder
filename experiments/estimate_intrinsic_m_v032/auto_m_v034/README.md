# Package automatic-M extension with MPCurver 0.3.4

This extension fits the 90 fixed inputs from the all-monotone study with the
current package-level automatic ordering-count initialization in MPCurver
0.3.4. The original adaptive EB and uniform-forward fits remain
unchanged and provide paired controls.

Each fit uses `intrinsic_dim = "auto"`, fixed-df natural-cubic-spline variance
explained (`spline_r2_df = 5`), single linkage, candidate cuts up to eight
orderings, and a minimum cluster size of two. The selected cluster-size
proportions initialize the adaptive global ordering probabilities. Subsequent
empirical-Bayes updates may reduce the effective ordering count.

All other fitting controls match the published study: 300 samples, 60
features, `K = 50`, RW2 smoothing, Isomap applied independently to every
selected feature group, normalized convergence tolerance `1e-6`, and at most
10,000 sweeps.
The seed is the original dataset seed plus 1,000,000.

The shared source archive at `experiments/mpcurve_v034_site/source/` is frozen
at MPCurver commit `15f2b0bbe5dfa61cd46da5160b2bc251e75a0475`;
its SHA-256 is recorded in the study design. Run one dataset from the
InferOrder root with:

```sh
bash experiments/estimate_intrinsic_m_smooth_v032/run_r.sh \
  experiments/estimate_intrinsic_m_v032/auto_m_v034/run_dataset.R 1
```

The Slurm array script maps tasks 0 through 89 to one-based manifest indices.
Completed compact results are retained under `results/`; checkpoints, full
replicate-one fits, installed packages, and scheduler logs are ignored.

## Completed run

Slurm array 59675733 completed all 90 tasks with exit code zero. All 90 fits
converged without warnings. Both the similarity cut and the final adaptive-EB
effective M recovered the generating count in all 90 datasets. Nine retained
replicate-one fits verify that every selected feature group requested and used
Isomap independently, with no PCA component.
