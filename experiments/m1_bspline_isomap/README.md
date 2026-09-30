# Common Isomap initialization for the replicated M1 comparison

This study replaces PCA with a common one-dimensional Isomap initialization
for MPCurve, principal curves, GPLVM and Bayesian GPLVM. It reuses the exact
observed matrices from the P=12 and P=50 random B-spline experiments: N=200,
M=1, noise SD 0.25 or 1 (SNR 16 or 1), and 30 replicates per condition.
There are 120 datasets and 480 fitted endpoints. The canonical report remains
`analysis/simulation_m1_comparison.Rmd`.

## Results

| Method | P=12 low | P=12 high | P=50 low | P=50 high |
| --- | ---: | ---: | ---: | ---: |
| Initial Isomap | 0.8687 | 0.6146 | 0.9997 | 0.8159 |
| MPCurve | 0.8811 | 0.6360 | 0.9997 | 0.8392 |
| Principal curve | 0.9649 | 0.5504 | 0.9997 | 0.8828 |
| GPLVM | 0.9095 | 0.6049 | 0.9997 | 0.8463 |
| Bayesian GPLVM | 0.9116 | 0.5433 | 0.9997 | 0.8572 |

Entries are median absolute Spearman recovery. At P=50 and low noise, counts
at or above 0.95 are 24/30, 29/30, 24/30 and 28/30 for the four fitted methods.
Principal curves' mean recovery is 0.9957 versus MPCurve's 0.9721 in that condition,
a paired difference of 0.0235 (Monte Carlo SE 0.0082). At P=50 and high noise,
principal curves have the highest median but a lower mean (0.8035) than MPCurve
(0.8340), reflecting substantial deterioration in several replicates.

All 480 endpoints are finite; convergence is 120/120 for MPCurve and GPLVM,
104/120 for principal curves (16 caps) and 115/120 for Bayesian GPLVM (five
line-search errors). One converged principal-curve run, p12_low_r22, emits 12
smoothing warnings and reports a numerical df=1 fallback despite requesting df=5.
It remains in summaries and appears as a triangle in the figures. All warnings
and actual initial projections are saved.

Including initialization, MPCurve median seconds are 0.2195/0.5200 at P=12 and
0.5105/1.4210 at P=50 (low/high noise). It is fastest in three conditions; GPLVM
is faster at P=50/high noise (0.7647 seconds). Median common embedding time is
0.04--0.08 seconds. Runtime tables retain fitting and total preparation-inclusive
times separately.

## Common upstream initialization

The frozen MPCurver 0.3.4 `isomap_ordering` implementation uses Euclidean distance
on the same column-centered observations, without feature variance scaling.
The neighborhood size is fixed at k=15, following the earlier pilot and package
default. Directed neighbor lists form an undirected union graph. All N=200
points are landmarks, and one-dimensional classical MDS uses the shortest-path
distance matrix. The package applies its usual inverse-distance landmark
extension (self-distance epsilon 1e-8), orients the coordinate by PC1, and
scales it to [0,1]. We then center it and normalize its population variance to
one before supplying it to the models. Orientation uses observations, not truth.

Every one of the 120 graphs is connected. Preparation explicitly requires this,
retains all 200 observations, and saves the complete geodesic matrix. There is
no neighborhood search, truth-based selection or disconnected-component removal.
The independent verification checks classical-MDS rank agreement against those
saved geodesics. The fixed k=15 study describes this choice, not optimal Isomap
performance over possible neighborhoods.

MPCurve maps the common Isomap ordering into 50 quantile bins, as for PCA.
Both GP methods receive the same continuous coordinates as their initial latent
locations (posterior means for Bayesian GPLVM). Bayesian variances and inducing
selection use the same seeds as the corresponding PCA runs.

Principal curves require a geometric starting curve. Following the pilot, each
feature is smoothed against the common Isomap coordinate using a df=5 spline,
evaluated in coordinate order. The package then projects observations onto that
curve before its first update. This projection can change sample ranks. Thus
all methods share an upstream Isomap start, but their internal representations
differ. `positions.csv` saves upstream, actual initial, and final coordinates;
`summary.csv` reports actual-initial recovery and upstream rank agreement.
The principal-curve initial projection is replayed independently during checks.

## Fitting and estimands

Package versions, kernels, smoothing settings, noise estimation, seeds, stopping
rules and iteration/evaluation budgets match the preceding PCA experiments.
MPCurver 0.3.4 uses K=50, RW2, learned feature-specific noise/precision, adaptive
position weights, relative ELBO tolerance 1e-6 and at most 2,000 sweeps.
princurve 2.1.6 uses df=5, stretch=2, tolerance 1e-6, at most 1,000 iterations.
GPy 1.13.2 uses one latent dimension: native RBF plus bias for GPLVM, native RBF
and 50 inducing points for Bayesian GPLVM. L-BFGS-B uses gtol=1e-5,
bfgs_factor=1e7, up to three blocks of 2,000 evaluations; only evaluation-limit
stops receive another block. No endpoint is selected using truth.

Absolute Spearman recovery allows a global reversal; absolute Kendall is a
secondary metric. Each final endpoint is paired with the same method's saved
PCA endpoint on identical observations. All finite endpoints are retained and
unmet stopping rules are flagged. Differences among methods are also paired
within the same Isomap dataset. Monte Carlo standard errors use replicate SD
of the difference divided by sqrt(30), and are descriptive for each condition.

Two R processes and two GP workers run with one BLAS thread per fit. Fit-time
measurement follows the earlier runs: R includes a zero-iteration fit used for
initialization auditing; Python includes its parameter archive write. Shared
Isomap preparation and principal-curve construction occur beforehand and are
recorded separately. `total_seconds` adds the common embedding time to each
method's fitting time, and adds smoothing construction for principal curves.
The reported medians of total time are calculated per replicate before taking
the median; they are not sums of marginal medians. Timings are descriptive
under concurrent processes and method-specific stopping criteria.

## Reproduction

From the InferOrder root, with the same frozen R libraries and Python environment
as the parent comparisons:

```bash
bash experiments/estimate_intrinsic_m_smooth_v032/run_r.sh experiments/m1_bspline_isomap/prepare.R
bash experiments/estimate_intrinsic_m_smooth_v032/run_r.sh experiments/m1_bspline_isomap/run_r_methods.R
experiments/estimate_intrinsic_m_smooth_v032/exploratory_m5_ordering_b/external_methods/.venv/bin/python experiments/m1_bspline_isomap/run_gpy.py --workers 2
bash experiments/estimate_intrinsic_m_smooth_v032/run_r.sh experiments/m1_bspline_isomap/summarize.R
bash experiments/estimate_intrinsic_m_smooth_v032/run_r.sh experiments/m1_bspline_isomap/build_site.R
python3 scripts/check_current_site.py
```

The R driver accepts replicate indices to split disjoint jobs; existing endpoints
are resumable checkpoints. `manifest.csv` records parent datasets, original input
hashes, graph components, and preparation time. `inputs/` retains Isomap/geodesics,
initial curves, the unchanged matrices and truth. `results/` retains full R fits,
GP parameters, coordinates, statuses and warnings. `metrics.csv`, `summary.csv`,
`method_comparison.csv`, `initialization_pairs.csv` and `positions.csv` retain
numeric comparisons. Source archives, driver hashes, dependency versions and
`verification.txt` document provenance and checks. Original PCA fits are retained.
