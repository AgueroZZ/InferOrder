# Estimate Intrinsic M Implementation Plan

**Goal:** Compare adaptive EB occupancy estimation and uniform-prior forward
selection on known intrinsic dimensions, and add a reproducible simulation
study to the InferOrder Simulation section.

**Architecture:** One independent experiment directory owns the pinned package
source, design manifest, generated data, resumable candidate fits, compact
results, validation records, and R Markdown report. Both methods share every
generated dataset. The workflowr page reads completed experiment summaries.

**Tech Stack:** R, MPCurver 0.3.1, clue, mclust, ggplot2, rmarkdown, workflowr.

**Execution:** Continue inline in the current task under the user's approved
design. The project-local `plan/` location follows the user's instructions.

## Confirmed main design

- True M = 3, 4, 5; n = 300; D = 60.
- Balanced groups: 20, 15, or 12 features per true ordering.
- Independent Uniform(0,1) latent positions, varied smooth monotone curves.
- Standardize each noiseless feature to sample variance 1, then add independent
  Gaussian noise with variance 1/SNR. SNR = 1, 4, 16.
- Ten independent datasets per condition: 9 conditions and 90 datasets.
  Reduced from 50 repeats at the user's request on September 26, 2026;
  original seed registry and completed pilot fits are retained.
- Same data for both methods; fixed M_max = 8 and K = 50.
- Similarity-only starts: absolute Spearman similarity, single linkage,
  within-group Isomap, with package singleton fallbacks recorded.
- Adaptive EB: one fit starting at M = 8, effective count at omega > 1e-12.
- Uniform forward: fit M = 1,2,...,8 and stop at the first non-improvement in
  converged final T = 1 ELBO, using the fixed uniform assignment prior.
- Adaptive position priors, RW2, ridge = 0, estimated feature noise and
  smoothness, default 5-to-1 annealing for partition fits.
- Relative convergence tolerance 1e-8. Continue the same fit in 1,500-sweep
  blocks, with at most 10,000 total sweeps per candidate. Unresolved fits
  remain reported failures; they are not silently accepted or replaced.
- Nine separate pilot datasets, one per condition, excluded from main totals.
  Check runtime and numerical behavior before starting the main manifest.

## Metrics and interpretation

Primary results use all 10 planned repeats as the denominator: exact-M recovery,
underestimation, overestimation, and unresolved/failed outcomes. Report the full
estimated-M distribution and Wilson intervals for exact recovery. For completed
fits report ARI for feature partitions, absolute-Spearman ordering recovery
after one-to-one matching, unmatched true orderings scored as zero, and runtime.
Report convergence/fit failures explicitly. EB threshold sensitivity at 1e-9,
1e-6, and 1e-3 is diagnostic and does not replace the declared primary count.

Methods are compared through truth recovery, not by comparing their ELBOs
against each other. Within the forward method, comparisons use identical K,
fixed-uniform priors, soft assignments, and final temperature 1.

## Files and steps

- [x] Create `experiments/estimate_intrinsic_m_v031/common.R` for pinned-library
  loading, design constants, data generation, and atomic checkpoint writes;
  separate `fit.R` and `evaluate.R` implement fitting and evaluation.
- [x] Create `prepare.R` to freeze distinct pilot/main seeds and save inputs.
  Validate dimensions, group sizes, unit noiseless variances, and noise scale.
- [x] Create `run_dataset.R` to run one manifest entry, save each candidate,
  continue unconverged fits, and write paired compact method results.
- [x] Create `run_study.R` as a resumable serial dispatcher with phase and
  optional row-range arguments for portable compute batches.
- [ ] Run and inspect all nine pilot datasets, including representative
  initialization metadata and convergence traces. Estimate main runtime.
- [ ] Freeze the final execution settings and complete all 90 main datasets.
- [x] Create `summarize.R`, `analysis_report.Rmd`, and `render_examples.R`.
  Use fixed representative datasets (replicate 1 at SNR 4) for recovery figures.
- [ ] Render and inspect the completed standalone HTML report.
- [ ] Add `analysis/estimate_intrinsic_m.Rmd`; update `analysis/_site.yml` and
  the Simulation links in `analysis/index.Rmd`.
- [ ] Build and inspect the workflowr page, key numerical tables, figures,
  math, and links. Record scientific results and deployment state in `log/`.

Commands run from the InferOrder root with `R_PROFILE_USER=/dev/null` and
`R_ENVIRON_USER=/dev/null`, using experiment-local package libraries. Each
worker uses one CPU thread. Pilot results determine the execution venue
and runtime estimate before the full study is dispatched.

Targeted stress scenarios (imbalanced groups, weak signal, nuisance features,
nonmonotone curves, correlated orderings) remain a subsequent design stage;
their numerical settings will be agreed after the main benchmark is reviewed.

## Runtime audit, September 26

All nine pilots completed (54 candidate fits, 18 method outcomes) in 9.30
single-worker fitting hours under tolerance 1e-8. The main benchmark has not
started. Investigate tolerance sensitivity and per-sweep hotspots before
launching the reduced ten-repeat study. Any revised fitting tolerance will
be recorded separately from the completed strict-tolerance pilots.
