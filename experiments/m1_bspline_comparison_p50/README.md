# Fifty-feature extension of the single-ordering comparison

This experiment increases P from 12 to 50 while retaining N=200, M=1,
30 replicates, random cubic B-spline trajectories, noise SD 0.25 and 1, and a
shared PCA initialization across MPCurve, principal curves, GPLVM and Bayesian
GPLVM. It compares ordering recovery and fitting time with the saved P=12 runs.
The canonical report remains `analysis/simulation_m1_comparison.Rmd`.

## Results

| Method | Low-noise median recovery | High-noise median recovery | Low-noise median seconds | High-noise median seconds |
| --- | ---: | ---: | ---: | ---: |
| Initial PCA | 0.698 | 0.695 | -- | -- |
| MPCurve | 0.693 | 0.707 | 0.26 | 1.22 |
| Principal curve | 0.496 | 0.568 | 3.32 | 2.42 |
| GPLVM | 0.686 | 0.701 | 2.53 | 0.64 |
| Bayesian GPLVM | 0.699 | 0.680 | 15.84 | 2.05 |

MPCurve and the two GP methods improve their median recovery relative to P=12;
the initial PCA recovery improves as well. Principal curves remain variable,
with 9/30 low-noise and 8/30 high-noise endpoints reaching 0.95, but 16/30 losses
larger than 0.05 relative to PCA at each noise level. Convergence: MPCurve 60/60,
principal curves 52/60 (eight caps), GPLVM 60/60, Bayesian GPLVM 59/60 (low_r10
line-search error). All 240 endpoints are finite and retained.

MPCurve remains fastest at low noise. GPLVM is fastest at high noise, where its
median evaluation count falls from 430 at P=12 to 152 at P=50. Relative to P=12,
median paired runtime ratios for MPCurve are 2.03 and 2.66, for principal curves
3.18 and 2.40, for GPLVM 0.78 and 0.36, and for Bayesian GPLVM 1.18 and 0.87
(low/high noise). These paired ratios differ from ratios of marginal medians.

## Paired construction

For each replicate, reuse the P=12 sample positions, its first 12 spline
coefficient sets and its first 12 standard-normal noise columns. Add 38
independent eight-coefficient cubic spline features and noise columns, using
`extra_seed = parent_seed + 1000000`. All 50 trajectories are centered and scaled
together to average signal variance one on the same 2,001-point grid. This
retains the same variance SNR (16 or 1). The first 12 coefficient sets and error
draws are unchanged; their signal amplitudes change by the common normalization
factor. The observed P=12 matrix is therefore not an exact submatrix of P=50.

PCA is recomputed from each entire centered, unscaled-feature P=50 matrix, then
shared across all four methods. Increasing P can change both the initial ordering
and the subsequent fit. Within each replicate, low/high noise continue to share
all trajectories, positions, and standard-normal errors. No data are selected
based on initialization quality or final recovery.

## Methods and timing

All fitting settings and versions match the P=12 experiment: MPCurver 0.3.4,
princurve 2.1.6, GPy 1.13.2; one PCA start per dataset; MPCurve K=50 and RW2;
principal-curve smoother df=5; Bayesian GPLVM 50 inducing points. Stopping
criteria and budgets are unchanged. The copied method drivers preserve the
pilot implementations, with the R input/output directory updated.

Every method estimates noise from the data. All finite endpoints are included,
with unmet stopping rules flagged. Recovery uses absolute Spearman correlation,
with absolute Kendall as a secondary check. Runtime measures model setup and
fitting, excluding shared simulation/PCA, file loading and serialization. The R
measurements also include a zero-iteration fit to record the initial state.
The Python measurement includes its parameter archive write, as in P=12.

There are two R processes, each processing a disjoint set of 15 replicates, and
two GP workers. Each fit uses a single BLAS thread. P=12 used one R process and
three GP workers. Thus runtime comparisons are descriptive, not a controlled
hardware scaling benchmark. The runtime table reports medians separately for
P=12 and P=50; the saved dimension summary also reports median within-replicate
runtime ratios, which can differ from the ratio of those two medians.

## Reproduction and artifacts

Run from the InferOrder root using the same frozen R library and Python
environment as the parent experiment:

```bash
bash experiments/estimate_intrinsic_m_smooth_v032/run_r.sh experiments/m1_bspline_comparison_p50/prepare.R
bash experiments/estimate_intrinsic_m_smooth_v032/run_r.sh experiments/m1_bspline_comparison_p50/verify_generator.R
bash experiments/estimate_intrinsic_m_smooth_v032/run_r.sh experiments/m1_bspline_comparison_p50/run_r_methods.R
experiments/estimate_intrinsic_m_smooth_v032/exploratory_m5_ordering_b/external_methods/.venv/bin/python experiments/m1_bspline_comparison_p50/run_gpy.py --workers 2
bash experiments/estimate_intrinsic_m_smooth_v032/run_r.sh experiments/m1_bspline_comparison_p50/summarize.R
bash experiments/estimate_intrinsic_m_smooth_v032/run_r.sh experiments/m1_bspline_comparison_p50/compare_dimensions.R
bash experiments/estimate_intrinsic_m_smooth_v032/run_r.sh experiments/m1_bspline_comparison_p50/build_site.R
python3 scripts/check_current_site.py
```

Existing fit files are resumable checkpoints. `inputs/` saves truth, generator
coefficients, noise, common coordinates, and parent input hashes. `results/`
saves full R fits and GP parameters, coordinates, stopping status and warnings.
`summary.csv` and `paired_summary.csv` compare methods at P=50;
`dimension_summary.csv` and `dimension_pairs.csv` compare P=50 with P=12.
`figures/` contains PNG and PDF artifacts. Source archive, seeds, file hashes,
package versions and dependency versions preserve provenance. The original
P=12 inputs and fits remain the comparison reference.
