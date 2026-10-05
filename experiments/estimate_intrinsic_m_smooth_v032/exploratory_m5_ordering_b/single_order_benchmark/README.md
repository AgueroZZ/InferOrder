# B as a single-ordering benchmark with native initialization

This benchmark isolates the 300-sample, 12-feature B signal and compares MPCurve,
P-curve, classical GPLVM and Bayesian GPLVM, each estimating one ordering from its
native initialization. P-curve recovers the ordering at zero and half noise when
allowed to continue beyond its default stopping point. All four methods leave
substantial ordering errors at original noise.

## Fixed input and protocol

The input is B from the fixed dataset `main_M5_S4_r001`. Feature membership is
supplied, and no feature assignment is learned. The three matrices reuse the same
signal and noise realization at noise multipliers 0, 0.5 and 1; multiplier 1 is
the original SNR 4 case. These are diagnostic variants of one selected dataset,
not independent simulation replicates. Input matrices and latent truth are saved
in `../external_methods/inputs/`. Each result's input SHA-256 is checked against
that matrix, and `source_manifest.csv` identifies the exact fit files reused.

| Method | Implementation | Native initialization | Main optimization settings |
| --- | --- | --- | --- |
| MPCurve | MPCurver 0.3.4, public `fit_mpcurve` | `method="PCA"`, centered and unscaled; native quantile discretization | `intrinsic_dim=1`, CAVI, up to 2000 sweeps, default tol=1e-6 |
| P-curve | princurve 2.1.6, `principal_curve` | `start=NULL`: centered, unscaled PC1 | maxit=1000, thresh=1e-6; default smoothing spline df=5 and stretch=2 |
| GPLVM | GPy 1.13.2, `GPLVM` | Native PCA, with feature standardization inside the initializer | L-BFGS-B, gtol=1e-5, bfgs_factor=1e7, up to three 2000-evaluation blocks |
| Bayesian GPLVM | GPy 1.13.2, `BayesianGPLVM` | Native PCA; native random variances and inducing-point selection | Same optimizer protocol, default 10 inducing points |

Default initialization and default stopping are separate choices. The main
comparison retains native starts and increases optimization budgets. For the GP
models a new block is run only after a budget-limited block, up to approximately
6000 evaluations in total. Package-default stopping results for MPCurve (100
sweeps) and P-curve (10 iterations, threshold 0.001) are retained separately.
There is no truth-based choice of initialization, seed, hyperparameter or endpoint.
The fixed seed is 20260929. GP observations are column centered, without variance
scaling; feature standardization occurs only inside native GPy PCA. Classical
GPLVM uses its default RBF plus bias kernel and Bayesian GPLVM its default RBF.

MPCurve was rerun through the public single-ordering API with its native
initializer. Its native K is 50 here; RW2, adaptive smoothness and position weights,
and other model settings use package defaults. This replaces the previous
ordering-injection helper for this benchmark, while preserving those earlier
fits. Frozen package source commit: `15f2b0bbe5dfa61cd46da5160b2bc251e75a0475`.
The public API's native initial moment calculation yields slightly different
iteration counts and a negligible half-noise difference from that helper.

## Results

Entries below are initial to final absolute Spearman correlations with truth.
Global reversal is treated as equivalent.

| Method | No noise | Half original noise | Original noise |
| --- | ---: | ---: | ---: |
| MPCurve | 0.103 -> 0.123 | 0.123 -> 0.097 | 0.142 -> 0.085 |
| P-curve | 0.103 -> 0.9998 | 0.123 -> 0.9988 | 0.142 -> 0.479 |
| GPLVM | 0.103 -> 0.150* | 0.096 -> 0.161 | 0.080 -> 0.168 |
| Bayesian GPLVM | 0.103 -> 0.158* | 0.096 -> 0.179 | 0.080 -> 0.143 |

The two starred noiseless GP runs exhaust their evaluation budgets. Their
endpoints are not converged optima. All other main runs satisfy their own
convergence criteria. All noisy methods have poor initial rank recovery; GP
methods' native initial ranks differ from MPCurve and P-curve because of the
PCA standardization. At original noise P-curve changes the ordering substantially,
but its result still contains large folds and is not a successful recovery.

[Main initial/final figure](native_initialization_benchmark.png)
([PDF](native_initialization_benchmark.pdf)). Every panel overlays actual initial
and final sample ranks; for MPCurve the initial series is the continuous PCA
ordering before its native grid discretization. Blue open circles show initial
ranks, orange points final ranks, and gray segments pair samples. Initial direction
is aligned to truth and final direction to initial; coordinate spacing is removed
by average ranks.

Stopping sensitivity is substantial for P-curve. With its package-default
stopping, final correlations are 0.172, 0.197 and 0.248, compared with 0.9998,
0.9988 and 0.479 under the main protocol. Default MPCurve gives 0.123, 0.097 and
0.085; the original-noise run hits its 100-sweep cap and converges after 154 sweeps
when the budget is raised. Principal-curve recovery at low noise is therefore
an iterative improvement from its native start under tighter stopping.

[Package-stopping sensitivity figure](package_stopping_sensitivity.png)
([PDF](package_stopping_sensitivity.pdf)). GP endpoints are unchanged in that
figure; it is not a claim that all four optimization budgets are package defaults.

`benchmark_summary.csv` contains the 12 main fits; `all_stopping_results.csv`
contains all 18 endpoints. Additional descriptive metrics are pairwise ordering
error (ties count as half errors) and mean absolute normalized-rank error, with
each estimate oriented to positive Spearman correlation with truth. The pairwise
metric underscores the remaining folds: P-curve at original noise still has
about 47.8% discordant/tied pair error despite its higher Spearman correlation.
These metrics evaluate ordering, not latent coordinate spacing or feature-curve
prediction. Likelihoods and ELBOs are not compared across methods.

## Reproduction and checks

Run from the InferOrder root:

```bash
bash experiments/estimate_intrinsic_m_smooth_v032/run_r.sh experiments/estimate_intrinsic_m_smooth_v032/exploratory_m5_ordering_b/single_order_benchmark/fit_mpcurve.R
bash experiments/estimate_intrinsic_m_smooth_v032/run_r.sh experiments/estimate_intrinsic_m_smooth_v032/exploratory_m5_ordering_b/single_order_benchmark/render_benchmark.R
```

The renderer reuses already-computed P-curve native-PCA fits and GP jobs 1–3 and
8–10 from [the external-method experiment](../external_methods/README.md), where
input creation, installation, fitting commands, pinned dependencies and literature
references are documented. Both P-curve stopping settings have the same default
start. The GP results have no injected Isomap or externally computed PCA coordinates.

Checks cover the 12 main and 18 total records, 300 finite initial/final positions
per fit, 300-by-12 input dimensions, source hashes, public MPCurve native-init
metadata, NULL P-curve starts and native-PCA GP starts with 10 Bayesian inducing
points. Both PNG figures were visually inspected. Fit warnings, optimizer statuses,
full MPCurve fits, exact versions and source references are retained. The experiment
is local and leaves package and website sources unchanged.
