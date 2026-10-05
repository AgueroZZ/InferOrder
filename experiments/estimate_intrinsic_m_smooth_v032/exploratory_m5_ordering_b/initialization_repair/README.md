# Can MPCurve repair an imperfect initial ordering?

On this selected challenging dataset, MPCurver 0.3.4 can substantially repair
randomly perturbed positions, but it does not reliably correct coherent
rearrangements such as folds, block reversals, and shifts. Most recovery in the
successful perturbation experiments occurs during the two subset CAVI updates
inside initialization, before the structural multi-ordering iterations begin.
The structural iterations can refine or worsen true-order recovery, and are not
a general mechanism for recovering the correct global order from arbitrary starts.

## Fixed data and experimental design

All experiments use the original `main_M5_S4_r001` observations: 300 samples,
60 features, true M=5, variance SNR=4. No signal or measurement noise is changed.
Only the initial positions for B's original 12-feature cluster are replaced;
the other four initial orderings are unchanged. M remains fixed at the correctly
selected five groups. Subsequent feature assignments and both adaptive priors
are learned normally. Truth is used to construct artificial errors and evaluate
recovery; the fitting updates do not receive it.

The predeclared 19 starts are:

- True positions, their global reversal, observed-data Isomap k=10 and k=15,
  and observed-data PC1.
- True positions plus independent Gaussian errors with SD 0.1, 0.3, or 0.6
  in unit-interval position units; three fixed random seeds per SD.
- Exchange the second and third rank quarters; reverse the middle half;
  circularly shift normalized ranks by one quarter.
- Two independent completely random starts.

`common.R` specifies every transformation and perturbation seed. Perturbations
are discretized by the existing quantile rule; no clipping of the Gaussian
perturbations is applied. Each start runs the ordinary two subset CAVI updates,
then zero-sweep expansion to all 60 features, followed by the original joint
structural fit (K=50, RW2, adaptive noise/smoothing, T=5 to 1 over 25 sweeps,
normalized tolerance 1e-6). Initialization and fit seed is the original dataset
seed plus 1,000,000. At most 10,000 sweeps are allowed.

The frozen package is 0.3.4, source commit
`15f2b0bbe5dfa61cd46da5160b2bc251e75a0475`. Per-run RDS files save the input hash,
package identity, fitting and perturbation seeds, and driver/common script hashes.
They also retain coordinates, initialization stages, every-sweep metrics,
selected position snapshots, final position probabilities, trajectory means,
noise/smoothing estimates, convergence, and warnings. A process-local exit trace
observes structural sweeps without changing their inputs or outputs. Instrumented
Isomap k=10, k=15, and PCA controls reproduce their previous objectives and recovery.
True-position and globally reversed controls have matching objectives.

## Results

All 19 fits converge without fitting warnings and retain feature-partition ARI 1.
Ranges below are minimum to maximum across the three perturbation seeds on
this one fixed dataset; they are not confidence intervals or independent datasets.

| Initialization | Raw absolute Spearman | After two subset CAVI updates | Final absolute Spearman |
| --- | ---: | ---: | ---: |
| True positions | 1.000 | 0.997 | 0.997 |
| Isomap k=10 | 0.995 | 0.996 | 0.997 |
| Gaussian error SD 0.1 | 0.937--0.956 | 0.996--0.996 | 0.997--0.997 |
| Gaussian error SD 0.3 | 0.641--0.779 | 0.964--0.995 | 0.959--0.997 |
| Gaussian error SD 0.6 | 0.400--0.483 | 0.790--0.928 | 0.789--0.961 |
| Swap adjacent middle quarters | 0.812 | 0.821 | 0.828 |
| Reverse middle half | 0.750 | 0.748 | 0.748 |
| Circular shift | 0.125 | 0.149 | 0.114 |
| Isomap k=15 | 0.750 | 0.745 | 0.700 |
| PCA | 0.142 | 0.139 | 0.145 |
| Random starts | 0.058 / 0.073 | 0.170 / 0.527 | 0.074 / 0.374 |

For example, SD=0.3 replicate 2 improves from 0.640842 to 0.996486, reducing
pair-order error from 27.35% to 2.56%. Its two initialization updates already
reach 0.986555. In contrast, the middle-half reversal starts at 0.750008 and
ends at 0.747599. Initial Spearman correlation alone therefore does not measure
the difficulty of repairing an ordering. The spatial organization of errors
matters: dispersed errors around a coherent trend can be repaired, while a
coherent but incorrectly connected ordering can remain self-consistent.

The true-order initialization is a diagnostic reference, not an attainable
perfect estimator under noisy observations: fitting takes it to about 0.997.
The two random starts do not recover the correct order. The current algorithm
is not supported here as an initialization-independent global ordering estimator.

Pair-order error counts discordant sample pairs, with ties contributing one
half. Global orientation is chosen by Spearman correlation with truth, consistently
with the other plots. Consequently, this pair-error metric can exceed 0.5 for
some nonlinear orders whose Spearman and Kendall orientation preferences disagree.
Mean absolute normalized-rank error is also retained in `repair_summary.csv`.

## Does the stopping rule explain the failures?

From each converged Isomap k=10, k=15, and PCA solution, run another 500 ordinary
T=1 structural sweeps with stopping disabled. The actual package sweep function
is used directly; updates, priors, and controls are unchanged.

| Initialization | Before extra sweeps | After 500 extra sweeps | Objective increase |
| --- | ---: | ---: | ---: |
| Isomap k=10 | 0.996823 | 0.996781 | 1.784827 |
| Isomap k=15 | 0.700139 | 0.700133 | 1.135286 |
| PCA | 0.145150 | 0.145150 | 1.136875 |

The global objective continues to improve slightly, but the incorrect B
orderings do not unfold. The final normalized objective increments are around
1.7e-8 to 5.1e-8, well below the original 1e-6 threshold. This is evidence against
ordinary early stopping as the cause of the ordering failures; it is not a
proof of a mathematical local optimum or a guarantee about infinitely many sweeps.

## What the implemented updates do

`frozen_update_functions.R.txt` records the functions from the installed frozen
package. Conditional on current trajectories, the position update evaluates all
50 positions for every sample. With estimated per-feature variance it is
proportional to

```
q(C_i = k) = pi_k * exp[-0.5 * sum_j w_j *
              ((X_ij - mu_jk)^2 + posterior_var_jk) / sigma2_j]
```

followed by normalization (terms constant in k omitted). There is no distance
penalty to the initial position, neighborhood restriction, or multiplicative
factor of the preceding position probability. A numerical floor preserves
positive probabilities; for K=50 it is approximately 2.98e-10 before normalization.
The dependence on the previous ordering comes through the newly fitted trajectory
posterior, the adaptive parameters, and the position prior.

The annealing temperature divides feature-assignment logits in the q(Z) update.
It is not supplied to q(C), the sample-position update. Current feature-assignment
annealing therefore is not position annealing and does not directly flatten
sample-position likelihoods to encourage global rearrangement.

A conditional-update diagnostic further isolates the issue: keep a poor fit's
noise estimates, position prior, and feature-assignment weights, but replace its
trajectory posterior with that from the good Isomap k=10 fit. One ordinary
position update changes B recovery from 0.700139 to 0.996583 for Isomap k=15,
and from 0.145150 to 0.996628 for PCA. Respectively 19.7% and 65.0% of samples
move more than one quarter of normalized rank. This is an oracle intervention
using a known good solution, not a proposed estimator. It demonstrates that the
position-update implementation permits large moves; it does not explain how
the optimizer would discover the good trajectory from a poor start.

## How the curves accommodate the wrong order

The generating design contains one monotone anchor, but the fitted model does
not know which feature is the anchor or impose monotonicity on it. Under the
bad orderings, V20 is fitted as a curved, nonmonotone function of inferred
position instead of correcting the ordering. Its estimated RW2 precision is
about 69,727 under Isomap k=10, 61.5 under k=15, and 12.2 under PCA (smaller
precision allows more curvature). Its estimated noise variance remains near
the true 0.25 in all three fits: 0.2460, 0.2499, and 0.2463. Thus this anchor's
ordering error is being accommodated primarily by a more flexible trajectory,
not simply by labeling its variation as extra noise.

The bad fits can be confident conditional on their learned curves: mean maximum
sample-position probability is about 0.927 for Isomap k=15 and 0.945 for PCA,
versus 0.691 for k=10 in the self-update audit. These conditional variational
probabilities do not represent uncertainty across alternative global orderings.
The good k=10 fit has the highest objective among these 19 completed fits, so
there is evidence for an optimization failure on this dataset, rather than a
claim that the objective necessarily prefers the displayed poor solutions.

## Interpretation and next research question

The evidence supports a real ability to denoise imperfect initial positions,
including large random perturbations, together with a substantial weakness in
repairing coherent global ordering errors. Most successful repair here is
inside initialization; it should not be attributed exclusively to the later
joint structural iterations. Feature-group recovery and high local posterior
confidence do not establish correct sample-order recovery.

A practical hypothesis supported by this case is to compare diverse data-only
starts by the same fitted objective: among observed-data Isomap k=5/10/15/20/30
and PCA, k=10 has the best completed objective and recovers B well. Whether this
is reliable across independent datasets is untested here. Position annealing
or coordinated block-reordering proposals would require separate implementation
and validation; current feature-assignment annealing is not evidence that either
already occurs. No package default is changed by this investigation.

This is one deliberately selected failure case with correct initial feature
groups. Three perturbation seeds characterize variation in starts, not sampling
uncertainty over datasets. Broader robustness claims require independently
replicated simulations with predeclared initialization errors and recovery metrics.

## Reproduction and figures

From the InferOrder root:

```sh
bash experiments/estimate_intrinsic_m_smooth_v032/run_r.sh \
  experiments/estimate_intrinsic_m_smooth_v032/exploratory_m5_ordering_b/initialization_repair/run_repair.R \
  1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19
bash experiments/estimate_intrinsic_m_smooth_v032/run_r.sh \
  experiments/estimate_intrinsic_m_smooth_v032/exploratory_m5_ordering_b/initialization_repair/check_stopping.R \
  isomap10 isomap15 pca
bash experiments/estimate_intrinsic_m_smooth_v032/run_r.sh \
  experiments/estimate_intrinsic_m_smooth_v032/exploratory_m5_ordering_b/initialization_repair/audit_updates.R
bash experiments/estimate_intrinsic_m_smooth_v032/run_r.sh \
  experiments/estimate_intrinsic_m_smooth_v032/exploratory_m5_ordering_b/initialization_repair/render_report.R
```

Cached compact outputs are retained under `results/`; the larger reconstructible
full fits used for continuation remain locally ignored under `full_fits/`.
If reproducing continuation in a fresh checkout, regenerate those full fits by
rerunning the corresponding baseline indices after moving their cached compact
result files aside.

- `repair_by_error_type.png`: raw, warmup, annealing, and final ordering recovery.
- `ordering_iteration_traces.png`: every-sweep recovery for representative starts.
- `forced_continuation.png`: 500 extra sweeps without the stopping rule.
- `anchor_trajectory_adaptation.png`: the same monotone anchor fitted under
  different final orderings, with true-position colors.
- CSV files retain all reported numerical comparisons and mechanism diagnostics.

All plots are exploratory standalone artifacts. Website sources, generated
pages, formal study fits, and package code remain unchanged.
