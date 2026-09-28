# Estimating M with one monotone feature and smooth trajectories

This study measures how estimating the number of sample orderings changes when
each ordering has one monotone feature and its remaining features follow smooth,
nonmonotone trajectories. It uses the same 90 datasets' latent positions, feature
groups, feature order, and measurement-noise realizations as the completed
`estimate_intrinsic_m_v032` monotone benchmark.

There are 300 samples and 60 features, with true M = 3, 4, or 5; variance SNR =
1, 4, or 16; and ten repetitions per condition. Features form balanced groups of
20, 15, or 12. The first current feature column within each true group retains
its original monotone signal. Each other feature follows

\[
f_j(t)=\sum_{k=2}^{4}
\frac{a_{jk}\sin(k\pi t)+b_{jk}\cos(k\pi t)}{k^2},
\qquad a_{jk},b_{jk}\overset{\mathrm{iid}}{\sim}N(0,1).
\]

These mixtures provide varied peaks and troughs over one to two oscillation
cycles. Coefficients are redrawn until the analytic derivative has both positive
and negative values on a fixed 1,001-position grid, using a relative threshold
of 1e-6 times its largest absolute value. Each sampled signal is centered and
scaled to variance one, then receives its original residual noise. The retained
monotone feature makes the complete noiseless feature vector vary uniquely with
the true position, while the other features have turning points.

The shape seed is the original dataset seed plus 2,000,000. Each saved dataset
contains the coefficients, normalization constants, derivative certificates,
anchor indices and names, exact residual-noise matrix, and hashes of the paired
baseline. `source_provenance.json` records the frozen MPCurver 0.3.2 source
archive's SHA-256, source commit, and baseline design hash.

## Methods

Adaptive EB fits eight candidate ordering slots. Uniform + forward fixes equal
ordering-assignment priors and adds slots until the converged temperature-one
ELBO first fails to improve, up to eight. Both methods report effective M by
counting posterior mean feature-assignment probabilities greater than 1e-12.

Fitting settings match the monotone study: MPCurver 0.3.2, CAVI, K = 50,
second-order random-walk smoothing with ridge = 0, absolute Spearman similarity,
single-linkage grouping, and Isomap initialization. The convergence threshold is
absolute ELBO change divided by N times D below 1e-6, with a maximum of 10,000
sweeps per candidate. Each candidate uses one initialization. Comparing paired
outcomes therefore measures the effect of the generating trajectories on this
complete fitting procedure, including its initialization.

## Results

Adaptive EB recovers the true effective M in 60/90 datasets, compared with 87/90
in the all-monotone benchmark. Uniform + forward recovers 88/90, compared
with 90/90. EB makes 29 underestimates and one overestimate; forward makes two
underestimates, both at true M = 5 and SNR = 1. All 180 method outcomes and 570 candidate
fits converged, with no fitting warnings.

Uniform + forward retains mean feature-partition ARI 0.995, while its mean
ordering score decreases from 0.985 to 0.843. Adaptive EB's corresponding ARI
and ordering scores are 0.730 and 0.685. Paired results, condition tables,
figures, and the standalone report are saved under `main_summary/` and
`analysis_report.html`.

## Reproduction

Run from the InferOrder root. Install the baseline's frozen
`experiments/estimate_intrinsic_m_v032/source/MPCurver_0.3.2.tar.gz` into this
study's `library/` directory with the required dependencies available. Set
R_PROFILE_USER and R_ENVIRON_USER to `/dev/null` and BLAS/OpenMP thread limits
to one, then prepare and validate inputs:

```sh
bash experiments/estimate_intrinsic_m_smooth_v032/run_r.sh experiments/estimate_intrinsic_m_smooth_v032/prepare.R
bash experiments/estimate_intrinsic_m_smooth_v032/run_r.sh experiments/estimate_intrinsic_m_smooth_v032/validate_setup.R
```

`manifest.csv` contains the fixed 90 main rows, and `dispatch_manifest.csv`
assigns task indices 0 through 89 for a Slurm array of single-thread tasks,
with at most 12 tasks running concurrently. `data/` retains the paired inputs;
`input_hashes.csv` and
`setup_validation.json` record reproducibility checks. Validation independently
regenerates every input and verifies the current monotone benchmark totals of
87/90 exact recoveries for adaptive EB and 90/90 for uniform + forward.

The included `run_r.sh` selects the verified R 4.3.3 environment in this workspace
and applies the clean, single-thread settings. After preparing inputs, submit
`slurm/run_array.slurm` from the repository root. Each array task executes
`run_dataset.R` for its assigned input, reusing completed candidates or saved
checkpoints. Finalize the study with:

```sh
bash experiments/estimate_intrinsic_m_smooth_v032/run_r.sh experiments/estimate_intrinsic_m_smooth_v032/validate_results.R
bash experiments/estimate_intrinsic_m_smooth_v032/run_r.sh experiments/estimate_intrinsic_m_smooth_v032/summarize.R main
bash experiments/estimate_intrinsic_m_smooth_v032/run_r.sh experiments/estimate_intrinsic_m_smooth_v032/render_report.R
```

`main_validation.json` records complete-result checks. The paired comparison
reevaluates saved baseline fits using the same occupied-slot matching rule,
then saves aligned baseline metrics and pair/condition/overall comparison tables.
The original baseline artifacts retain their published values.
