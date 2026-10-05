# Can the first ELBO select the Isomap neighborhood?

The first joint structural update is not a reliable selector on this case.
At the original noise level it selects k=15, while the best converged candidate,
k=10, ranks last among k=5/10/15/20/30 and PCA. At sweep 2 the preferred candidate
is k=30. At sweeps 3, 5, 10, and 25 it is k=10.

All these comparisons use the identical M=5, SNR=4 observations and model.
Only B initialization changes, following the controlled fits in the parent
experiment. A sweep here is a joint structural update after the two subset
CAVI initialization updates and the zero-sweep full-data expansion; it is not
a raw Isomap score or the initial objective before any structural update.

## Scores on the original observations

Larger ELBO is better. These are ordinary T=1 ELBO evaluations of each saved
annealed state, including at the first update where the fitting temperature is 5.

| Initialization | Sweep 1 | Sweep 5 | Sweep 10 | Converged | Final B correlation |
| --- | ---: | ---: | ---: | ---: | ---: |
| k=5 | -28744.38 | -18629.96 | -18084.95 | -17839.41 | 0.996793 |
| k=10 | -28746.93 | -18626.80 | -18083.91 | -17839.19 | 0.996823 |
| k=15 | -28660.99 | -18748.23 | -18222.89 | -17989.57 | 0.700139 |
| k=20 | -28700.33 | -18749.94 | -18224.30 | -17997.36 | 0.688720 |
| k=30 | -28686.55 | -18720.63 | -18217.03 | -17982.78 | 0.676222 |
| PCA | -28726.95 | -18760.51 | -18244.98 | -17995.31 | 0.145150 |

Selecting at the first update loses 150.38 final ELBO units relative to k=10;
selecting after two updates loses 143.59. Even retaining the top two candidates
at the first update would discard k=10 and k=5. By five updates k=10 and k=5
are the top two candidates and both produce high final recovery.

## Temperature and exact replay checks

The original schedule is T=5 to 1 over 25 joint sweeps. The reported optimization
objective at T>1 includes T times the feature-assignment entropy. It is therefore
not the ordinary ELBO. For each state we recompute assignment information at T=1
and reevaluate the package objective without changing any variational state.
We also verify the equivalent correction

```
standard_ELBO = annealed_objective - (T - 1) * entropy(q(Z)).
```

This correction does not fix the first-step misranking. At the original noise
level the first-step corrections for k=10 and k=15 are about 30.23 and 36.54,
respectively, and k=15 still scores higher.

`score_early.R` replays only the original 25-sweep annealing phase, with no
post-annealing convergence iterations. It covers all 17 previously completed
candidates across B-noise multipliers 0, 0.5, and 1; the disconnected noiseless
k=5 candidate remains excluded. Every replay's 25 logged objective values match
the corresponding prefix of the saved complete fit within 1e-5. Modified input
hashes are checked, and seeds and frozen package commit are retained in the
compact per-candidate results. No fitting model, package source, or original
fit is modified.

## Interpretation

A short-run score is a useful screening hypothesis, but a one-step winner need
not be the converged winner: different starts have different transient parameter
and assignment adjustments and different remaining objective gains. Monotone
improvement along an individual optimization path would not imply preservation
of rankings between paths.

On the noiseless B signals all connected Isomap candidates are tied to numerical
precision. At half noise, the best Isomap candidates converge to nearly identical
scores: selecting k=30 instead of k=10 costs only 0.000195 ELBO units, and their
final B recovery is essentially identical. Rankings among these near-ties should
not be interpreted as meaningful differences in initialization quality.

A candidate practical strategy is equal-budget short runs, scoring each state
with the ordinary ELBO, retaining multiple competitive candidates, and continuing
them with their intended optimization schedule. Five to ten sweeps is motivated
by this example, not a validated universal budget or rejection threshold. Both
this case and its noise-scaled versions share the same data-generating realization;
independent datasets are required to validate selection accuracy and runtime.
Initialization and scoring overhead must be included in any claimed speedup.

## Reproduce

From the InferOrder root:

```sh
bash experiments/estimate_intrinsic_m_smooth_v032/run_r.sh \
  experiments/estimate_intrinsic_m_smooth_v032/exploratory_m5_ordering_b/early_selection/score_early.R \
  1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18
bash experiments/estimate_intrinsic_m_smooth_v032/run_r.sh \
  experiments/estimate_intrinsic_m_smooth_v032/exploratory_m5_ordering_b/early_selection/summarize.R
```

Indices 1--15 are the original neighborhood/noise grid; 16--18 are PCA at
noise multipliers 0, 0.5, and 1. Index 1 is skipped as disconnected.
`early_scores.csv` retains every evaluated state; `selection_by_budget.csv`
records the selected candidate and eventual score loss; `comparison.csv`
provides the displayed score table. `early_elbo_selection.png`/PDF plots the
ordinary-ELBO gap to the best candidate at the same sweep on the original data.
The plot compares scores within each sweep, rather than subtracting objectives
at different temperatures. Website files remain unchanged.

## Distinction from the user's intended within-group selector

The user subsequently clarified that the intended score is obtained by fitting
one single-ordering update to the group's own features, before joint fitting.
That experiment is in `../within_group_selection/`. It selects k=10 on the
original B data and succeeds at choosing a high-recovery candidate in all three
noise settings. The joint-update first-step failure documented here does not
refute that different local scoring proposal.
