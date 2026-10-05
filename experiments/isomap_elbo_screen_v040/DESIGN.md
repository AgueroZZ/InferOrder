# One-sweep Isomap screening and continuation

This internal experiment asks whether a single MPCurve CAVI sweep can select
an Isomap neighborhood that repairs the previously identified folded ordering,
while costing less than fitting every candidate to convergence.

The fixed input is ordering B from `main_M5_S4_r001` in
`../estimate_intrinsic_m_smooth_v032/data/`: 300 samples by 12 features,
variance SNR 4, original observations and noise. The group is supplied from the
saved truth for this isolated single-ordering experiment. No grouping or
multi-ordering optimization is performed. This is a selected failure case,
not a new independent dataset or a robustness study across sampling densities.

The primary candidate set is every integer `num_neighbors` from 5 through 30.
A prespecified cheaper comparison uses 5, 10, 15, 20, 30. Both include the
default-neighborhood comparator 15. Fix `num_bins = 50`, RW2, ridge zero,
initial precision one, adaptive precision/noise/position probabilities,
quantile discretization, no feature standardization, and no supplied noise SDs.
Use the public MPCurver 0.4.0 API from the adjacent source repository.

For each candidate, run `fit_mpcurve(initial_method = "isomap", max_iter = 1,
tol = 0)` and select the highest final ELBO in the two-entry trace. Exact ties
are resolved by the smallest neighborhood. Initialization includes the usual
initial trajectory-posterior calculation; X=1 means one CAVI sweep after that
initial state. Use all 300 landmarks, embedding component one, and strict
initializer failure handling. Disconnected or failed candidates are excluded
and recorded, with no PCA fallback or missing-position imputation.

Continue the selected saved fit with `do_mpcurve()` at tolerance 1e-6,
normalized by 300 times 12 observations, up to 10,000 total sweeps. Continue
every other connected candidate using the identical controls only to audit the
selection against the highest converged candidate ELBO. These extra fits do
not enter the deployed screen-and-continue cost or the selection decision.
Also fit neighborhoods 15 and the selected neighborhood uninterrupted to
check continuation against the ordinary public fitting pipeline.

Recovery is absolute Spearman correlation between true positions and posterior
mean positions. Truth is used only for defining this known feature group and
reporting recovery. The selector sees the observation matrix and early ELBOs.
Plots align global orientation and use normalized ranks to show reordering.

Measure complete in-memory, sequential, single-thread pipelines on this host:
fixed k=15, dense one-sweep screening plus selected continuation, sparse
screening plus selected continuation, and dense fitting of all candidates to
convergence. Exclude package loading, file I/O, truth metrics, and plots from
all timed pipelines. Include embedding and initial posterior calculation in
every fit. Warm each pipeline once, then time five repeats in randomized order
(seed 20261002). Retain every timing, stage durations, selected k, final ELBO,
recovery, iterations, convergence, and the ordering of timed runs. Repetition
quantifies runtime variability on this dataset, not scientific replication.

Save compact states and histories, input/package/script SHA256 hashes, package
and dependency versions, source commits, seeds, controls, and checks. Keep
reconstructible full fitting objects and installed libraries locally ignored.
Do not change package defaults, existing results, or public website pages.
