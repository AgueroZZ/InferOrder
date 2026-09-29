# Estimating M with one monotone feature and smooth trajectories

This study measures how estimating the number of sample orderings changes when
each ordering has one monotone feature and its remaining features follow smooth,
nonmonotone trajectories. It uses the same 90 datasets' latent positions, feature
groups, feature order, and measurement-noise realizations as the completed
`estimate_intrinsic_m_v032` monotone benchmark.

The public workflowr report presents only the current MPCurver 0.3.4 automatic-M
fits under `auto_m_v034/`. Earlier package fits are retained here as experiment
provenance and are not displayed as current software performance.

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

## MPCurver 0.3.3 automatic-M extension

The `auto_m_v033/` extension applies the released package-level initializer to
the same 90 fixed input matrices. It uses fixed-df natural-cubic-spline
variance explained (`spline_r2_df = 5`), single linkage, cuts up to M = 8, and
a minimum cluster size of two. If the selected group sizes are
`d_1, ..., d_M`, adaptive EB is initialized with global assignment
probabilities `d_m / 60`.

The similarity cut estimates the true M in 88/90 datasets, with two
overestimates and no underestimates. After adaptive EB, the effective M is
correct in 88/90, with two underestimates and no overestimates. Compared with
the original adaptive M = 8 fits, the extension gains exact recovery in 28
paired datasets and loses it in none. Its mean feature-partition ARI is 0.994,
mean ordering recovery is 0.839, and median runtime is 32.5 seconds. Uniform +
forward is also exact in 88/90, with mean ARI 0.995, mean ordering recovery
0.843, and median runtime 132.4 seconds.

All 90 extension fits converged without warnings. `auto_m_v033/validation.json`
records strict input, source, result, and occupancy-threshold checks. The
extension freezes MPCurver 0.3.3 commit
`c905901424e43eab78b58bdcc0d1de367ec8fd73` and source archive SHA-256
`cb57e1f6e8d859b7ff9aeda17c97d538fae0686eb0ceb616f126023a26d645f9`.

## Exploratory ordering-B diagnostic

`exploratory_ordering_b/` diagnoses the ordering-B correlation of 0.756 in the
original adaptive fit for `main_M3_S4_r001`. The feature partition is exact,
but the old M = 8 absolute-Spearman initialization splits B into groups of 3,
16, and 1 features. The surviving 16-feature Isomap ordering already has
correlation 0.757, and adaptive fitting preserves that local solution. A
truth-aware sensitivity analysis identifies one excluded trajectory that
restores the Isomap geometry; the complete B group, uniform-forward fit, and
automatic-M fit all recover B at about 0.998. This diagnostic does not alter
the completed benchmark or propose a truth-dependent fitting rule.

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

To reproduce the automatic-M extension after installing its frozen source
archive into `auto_m_v033/library/`, run the 90 array tasks and then validate
and summarize:

```sh
sbatch experiments/estimate_intrinsic_m_smooth_v032/auto_m_v033/slurm/run_array.slurm
bash experiments/estimate_intrinsic_m_smooth_v032/run_r.sh \
  experiments/estimate_intrinsic_m_smooth_v032/auto_m_v033/validate_results.R
bash experiments/estimate_intrinsic_m_smooth_v032/run_r.sh \
  experiments/estimate_intrinsic_m_smooth_v032/auto_m_v033/summarize.R
```
