# Package automatic-M extension with MPCurver 0.3.3

This extension fits the 90 fixed inputs from the one-monotone-anchor study with
the package-level automatic ordering-count initialization released in
MPCurver 0.3.3. The original adaptive EB and uniform-forward fits remain
unchanged and provide paired controls.

Each fit uses `intrinsic_dim = "auto"`, fixed-df natural-cubic-spline variance
explained (`spline_r2_df = 5`), single linkage, candidate cuts up to eight
orderings, and a minimum cluster size of two. The selected cluster-size
proportions initialize the adaptive global ordering probabilities. Subsequent
empirical-Bayes updates may reduce the effective ordering count.

All other fitting controls match the published study: 300 samples, 60
features, `K = 50`, RW2 smoothing, Isomap for a possible one-ordering
fallback, normalized convergence tolerance `1e-6`, and at most 10,000 sweeps.
The seed is the original dataset seed plus 1,000,000.

The source archive under `source/` is frozen at MPCurver commit
`c905901424e43eab78b58bdcc0d1de367ec8fd73`; its SHA-256 is recorded beside
the archive. Run one dataset from the InferOrder root with:

```sh
bash experiments/estimate_intrinsic_m_smooth_v032/run_r.sh \
  experiments/estimate_intrinsic_m_smooth_v032/auto_m_v033/run_dataset.R 1
```

The Slurm array script maps tasks 0 through 89 to one-based manifest indices.
Completed compact results are retained under `results/`; checkpoints, full
replicate-one fits, installed packages, and scheduler logs are ignored.
