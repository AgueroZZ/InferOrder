# GP precision substitution on ordering B

This exploratory sanity check asks whether replacing MPCurve's RW2 precision with the precision of a GP kernel makes its inferred ordering resemble Bayesian GPLVM. The same 300 observations and 12 B features from `main_M5_S4_r001` are used at noise multipliers 0, 0.5, and 1. No website or package source was changed.

Half-noise results are strongly concordant after aligning additional model settings: GP-kernel CAVI and Bayesian GPLVM have ordering agreement 0.983 (absolute Spearman correlation), while their truth correlations remain only 0.169 and 0.158. Agreement is weaker at no noise and original noise (0.902 and 0.828). This supports related behavior in this example, not equivalence of the inference procedures or successful recovery.

## Kernel and comparison design

The installed GPy 1.13.2 BayesianGPLVM default is RBF with ARD, not Matern. Its latent prior is standard normal, its observation variance is shared across features, and its zero-mean GP kernel amplitude and lengthscale are shared. Initial latent means use GPy's PCA implementation; latent variances are initialized uniformly in (0,0.1). The exact installed source was inspected. The main comparator reuses the earlier saved native-PCA, seed-20260929 fits with 50 inducing points. The constructor default is 10 inducing points; 50 is an explicit approximation refinement, not a default. Earlier 10-point fits remain available in `../external_methods/gpy_results/`.

Two experiments separate precision substitution from broader model alignment:

1. **Q replacement:** retain the original 50-bin MPCurve PCA responsibilities, feature-specific adaptive noise and prior scales, and adaptive position probabilities. Set Q to the inverse of the RBF or Matern-3/2 covariance on 50 equally spaced coordinates in [-sqrt(3),sqrt(3)]. This fixes a unit-scale coordinate convention without using truth. Initial kernel variance and lengthscale are 1; subsequent feature amplitudes are 1/lambda_j, while lengthscale stays fixed. Initial noise and position probabilities are taken from the native fit. RW2's initial lambda values are not transferred because their scale has a different meaning. A relative diagonal nugget of 1e-6 ensures stable inversion. Proper GP priors use rank 50 and ordinary log determinant, rather than RW2's rank 48 and pseudo-determinant.
2. **Aligned grid GP:** use the same centered data, GPy initial latent means and variances, standard-normal position prior, zero-mean RBF covariance, shared amplitude/lengthscale, and shared noise variance as Bayesian GPLVM. Discretize positions onto 100 points in [-4,4], initializing responsibilities by integrating the native initial Gaussian distributions over grid cells. Use Gaussian curve updates and analytic softmax position updates. Every five sweeps (and initially), optimize the three shared log hyperparameters, conditional on current responsibilities and with the curve Gaussian optimized analytically. This small hyperparameter optimization does not replace the analytic update of the 30,000 position responsibilities. Kernel/noise bounds are amplitude [exp(-8),exp(5)], lengthscale [0.03,10], and noise variance [exp(-12),exp(4)]; these numerical bounds differ from GPy's positive-only constraints. Stop when a hyperparameter-update sweep improves ELBO by less than n*d*1e-6, with a 2000-sweep cap. Half-noise controls use 200 grid points and a Matern-3/2 kernel.

The aligned experiment still differs from Bayesian GPLVM: discrete versus Gaussian position factors, finite latent support, a mean-field q(positions)q(grid function values) versus GPy's inducing-variable variational construction, and different optimization trajectories. It is not obtained solely by replacing Q in the public API. No new Bayesian GPLVM Matern benchmark was fitted; the Matern results are kernel sensitivity controls against the saved RBF comparator.

## Results

All correlations allow a global ordering reversal. Recovery compares truth to posterior expected position; agreement compares the two methods' final positions.

| Noise multiplier | Native RW2 truth rho | RBF Q replacement truth rho | Matern Q replacement truth rho | Aligned RBF CAVI truth rho | Bayesian GPLVM truth rho | Aligned CAVI/BGPLVM agreement |
|---|---:|---:|---:|---:|---:|---:|
| 0 | 0.123 | 0.069 | 0.110 | 0.222 | 0.158 | 0.902 |
| 0.5 | 0.097 | 0.121 | 0.117 | 0.169 | 0.158 | 0.983 |
| 1 | 0.085 | 0.031 | 0.136 | 0.064 | 0.151 | 0.828 |

For half noise, RBF Q substitution alone has agreement 0.958 with Bayesian GPLVM, versus native RW2's 0.906. The aligned 200-grid result has truth rho 0.177 and agreement 0.978, supporting qualitative agreement but not numerical grid convergence. Aligned Matern-3/2 has truth rho 0.164 and agreement 0.978 with the RBF Bayesian GPLVM. All primary fits satisfy their specified stopping rule; all saved objective traces are monotone to numerical precision. These are local optimizer results. Similarity is partly explained by shared PCA initialization and the same folded ordering basin; high mutual agreement alone does not prove model equivalence.

The saved timings are exploratory single runs, not a controlled speed benchmark. Sparse precisions are represented and solved densely here to isolate approximation accuracy, so these runs do not test sparse-solver speed gains. Do not compare ELBO values across the intrinsic RW2, proper GP, and distinct variational models as if they were the same objective.

## Sparse SPD precision

A sparse precision is constructed from ordered nearest-predecessor GP conditionals. If B is the unit lower-triangular matrix of conditional regression coefficients and D contains conditional variances, Q_sparse = B' D^-1 B is SPD by construction. This is a Vecchia-style covariance approximation, not arbitrary thresholding of Q and not an SPDE finite-element implementation. It retains the same analytic CAVI updates.

For 50 nodes and five predecessors, the precision has 20.8% nonzero entries. Relative covariance Frobenius errors are 0.288% for Matern-3/2 and 25.018% for RBF. The Matern ordering results remain essentially unchanged across all noise settings. The five-neighbor RBF result is therefore a coarse approximation and is not used as evidence of equivalence. Increasing the RBF neighborhood to 20 reduces covariance error to 0.393% (65.2% nonzero precision entries); half-noise truth rho is 0.12136 versus exact 0.12125. This illustrates a kernel-dependent accuracy/sparsity tradeoff.

## Verification

`verification.json` checks RBF and Matern kernel values against installed GPy, and checks Gaussian posterior mean/covariance and conditional evidence against exact GP regression when positions are known grid nodes. Maximum errors are below 5e-14. This verifies the mathematical GP-prior substitution independently of latent-position recovery.

`audit.R` creates a local copy of the frozen MPCurver 0.3.4 CAVI function and overrides only its precision constructor in a private function environment. The package's own metadata routine recognizes the proper GP precision. It runs all 15 exact/sparse Q replacements using identical starts and saves their R fits. `package_audit.csv` shows maximum responsibility discrepancy below 3e-9 and maximum ELBO discrepancy below 1e-7 against the independent implementation. Very small noiseless rank differences can result from floating-point near ties even when probabilities agree to this precision. The prior comparison's native RW2 CAVI validation is retained under `../collapsed_comparison/`.

## Artifacts and reproduction

- `run.py`: kernels, sparse precision construction, Q replacement fits, aligned-grid fits, and GP checks.
- `export_audit.py`, `audit.R`: independent frozen-package validation and saved `_package.rds` fits.
- `report.py`, `summary.csv`, `fit_checks.json`: recovery, ordering agreement, and convergence checks.
- `ordering_comparison.png`: gray CAVI initialization, blue Bayesian GPLVM, orange GP-kernel CAVI; top row replaces Q, bottom row additionally aligns model settings. Bayesian GPLVM is oriented positively toward truth; other series are oriented toward Bayesian GPLVM to remove arbitrary reversal. Titles report absolute truth correlation.
- `ordering_agreement.png`: final aligned CAVI versus Bayesian GPLVM ranks; color denotes true latent position. Both figures were visually inspected.
- Per-fit JSON/NPZ files contain settings, histories, final states, and script hashes. `provenance.json` includes input hashes, software versions, package source commit, and comparator paths.

Run from the InferOrder root using the existing `experiments/estimate_intrinsic_m_smooth_v032/exploratory_m5_ordering_b/external_methods/.venv/bin/python` interpreter. Execute `gp_kernel_sanity/run.py`, then invoke `run.swap(case, 'rbf', 20)` for each of B_noise0, B_noise05, and B_noise1 to reproduce the larger-neighborhood follow-up. Run `report.py` and `export_audit.py`. Execute `audit.R` with `bash experiments/estimate_intrinsic_m_smooth_v032/run_r.sh` followed by its experiment-relative path. Do not import NumPy before `run.py` sets single-thread environment variables if reproducing timing.
