# Standard principal-curve and GPLVM comparisons for the challenging M = 5 case

This exploratory comparison asks whether established one-dimensional curve methods
repair the imperfect initialization seen for ordering B in Simulation 3. Classical
principal curves can unfold the poor PCA initialization at zero and half noise,
where MPCurve remains folded. At original noise, GPLVM improves the same problematic
Isomap start more than either single-ordering or joint MPCurve, but does not fully
recover B. All methods recover B well from the better Isomap k = 10 start.

## Methods and literature

The comparators are existing package implementations of established methods:

- **Hastie–Stuetzle principal curves**, R **princurve 2.1.6**, using
  `principal_curve`. [Hastie and Stuetzle (1989), JASA](https://doi.org/10.2307/2289936);
  [maintainer documentation](https://rcannood.github.io/princurve/reference/principal_curve.html).
  The native start is centered, unscaled PC1. The default smoother uses a smoothing
  spline with df = 5; updates smooth each feature along the current arc length and
  reproject observations onto the resulting curve. We retain package defaults
  (`maxit=10`, `thresh=0.001`, `stretch=2`) and separately use `maxit=1000`,
  `thresh=1e-6`, keeping the smoother unchanged.
- **Classical GPLVM**, **GPy 1.13.2**, `GPy.models.GPLVM`, one latent dimension,
  native RBF plus bias kernel and Gaussian likelihood. See
  [Lawrence (2005), JMLR](https://www.jmlr.org/papers/v6/lawrence05a.html) and
  [the Sheffield GPy implementation](https://github.com/SheffieldML/GPy).
- **Bayesian GPLVM**, the same GPy release, `GPy.models.BayesianGPLVM`, one latent
  dimension and native RBF kernel. See
  [Titsias and Lawrence (2010), AISTATS](https://proceedings.mlr.press/v9/titsias10a.html).
  We check the native 10 inducing points and an explicit 50-point setting; matched
  initialization comparisons use 50. GPy initializes latent variances randomly
  and chooses inducing points from initial latent locations. Seeds and parameters
  are saved; two additional seeds check the original-noise native PCA case.

GPy is a standard implementation from the GPLVM research community; this choice
is not a claim that it is the only implementation or the latest available release.
The pinned version has a binary compatible with the available Python 3.9 environment.

Both GPy models use PCA by default, but GPy's PCA **standardizes each feature**
internally. This is visible in its [PCA implementation](https://gpy.readthedocs.io/en/deploy/_modules/GPy/util/pca.html)
and was checked against the installed 1.13.2 source saved here. Consequently,
native GPy PCA and the earlier R PCA are separate controls. Jobs 20–25 explicitly
supply the earlier centered, unscaled PCA coordinates to GPy, normalized to latent
mean zero and variance one. Their initial ranks agree exactly with the R controls.

## Design and estimand

The fixed input is `main_M5_S4_r001`: 300 samples, 60 features, five true feature
groups, with 12 features in B. Input columns are selected by `true_assign`, not by
physical column ranges. Comparators receive only the same 300-by-12 B matrix and
estimate one ordering. This tests ordering recovery conditional on the correct
feature group; it does not test discovery of multiple feature groups. Truth is
used for group selection and evaluation, not to supply latent positions or select
a favorable optimization endpoint.

Noise scales 0, 0.5 and 1 use the same signal and the same noise realization as the
previous analysis. Scale 1 is the original SNR 4 case. GP observations are column
centered without variance scaling; principal curves and MPCurve receive the raw
matrix. All comparisons use absolute Spearman correlation, treating a global
reversal as equivalent. Figures show average ranks to remove coordinate stretching.

The earlier joint MPCurve fits retain all 60 features and M = 5. We additionally
fit frozen **MPCurver 0.3.4**, source commit
`15f2b0bbe5dfa61cd46da5160b2bc251e75a0475`, to B alone using its existing
`.cavi_build_from_ordering` initialization: K = 50, RW2, quantile bins, initial
lambda = 1, adaptive position weights and smoothness, tolerance 1e-6, maximum
2000 sweeps. All nine single-ordering fits report convergence without fit warnings.
This separates within-group ordering estimation from joint feature assignment.

GPy uses L-BFGS-B with `gtol=1e-5`, `bfgs_factor=1e7`, and a budget of 2000 per
optimization block. Runs ending at the evaluation limit receive up to two more
blocks; actual evaluation counts are saved (SciPy can exceed the nominal limit
slightly). Neither objective values nor training errors are compared numerically
between model families: likelihoods, ELBOs and curve projection distances differ.

Principal curves require a geometric starting curve, rather than latent
coordinates alone. For Isomap starts, each feature is smoothed against the
Isomap coordinates with df = 5, and the package projects observations onto this
curve before its first update. The figure records this **actual initial
projection**, not the upstream Isomap positions. This extra step changes original
k = 15 correlation from 0.749564 to 0.803014. These are common upstream starts,
not identical internal states across algorithms.

## Results at original noise

Final absolute Spearman correlations with true B, using common upstream starts:

| Method | Centered, unscaled PCA | Isomap k = 15 | Isomap k = 10 |
| --- | ---: | ---: | ---: |
| Initial upstream ordering | 0.141566 | 0.749564 | 0.995146 |
| MPCurve, five orderings | 0.145150 | 0.700139 | 0.996823 |
| MPCurve, B only | 0.084701 | 0.758614 | 0.996760 |
| Principal curve, package defaults | 0.248186 | 0.870330 | 0.996211 |
| Principal curve, tighter stopping | 0.479039 | 0.221100 | 0.996214 |
| Classical GPLVM | 0.178035 | 0.806974 | 0.996780 |
| Bayesian GPLVM, 50 inducing points | 0.118477 | 0.804596 | 0.996698 |

All final fits in this table satisfy their own stopping rules. This is not a
certificate of a global optimum or correct ordering. From native GPy PCA instead,
classical GPLVM gives 0.167867 and Bayesian GPLVM gives 0.143159 (10 inducing) or
0.151014 (50 inducing); the two extra 50-point seeds give 0.132559 and 0.151015.
Their native initial correlation is only 0.080309. Thus native PCA does not solve
this case, and the matched-PCA controls do not change that conclusion.

The principal-curve k = 15 result is strongly stopping-dependent. The default
stops at iteration 5 with correlation 0.870330. Tighter stopping reaches iteration
48 and correlation 0.221100. Continuing reduces the final projection error but
worsens ordering recovery; the path is not monotonically improving in either
quantity. Selecting the iteration with best truth correlation would be an oracle
procedure and is not used here.

[Initial/final ordering panels at original noise](original_noise_matched_starts.png)
([PDF](original_noise_matched_starts.pdf)).

## Repair from exactly the same PCA ordering

| Noise multiplier | Initial PCA | Joint MPCurve | B-only MPCurve | Principal curve, tighter stopping |
| --- | ---: | ---: | ---: | ---: |
| 0 | 0.103159 | 0.122720 | 0.123470 | 0.999777 |
| 0.5 | 0.123075 | 0.153552 | 0.097074 | 0.998820 |
| 1 | 0.141566 | 0.145150 | 0.084701 | 0.479039 |

The low-noise principal-curve success occurs after initialization, mostly around
iterations 20–30. At half noise the default convergence threshold stops at
iteration 9 with correlation 0.197190; tighter stopping reaches 0.998820 at
iteration 49. At zero noise the default hits its 10-iteration cap with correlation
0.171665; the tighter run converges at iteration 47. This is direct evidence that
another established update scheme can repair this folded PCA ordering, rather
than evidence of a superior initial embedding.

The classical and Bayesian GPLVM PCA controls remain folded. Some zero-noise
GP fits do not converge: classical native and matched PCA and Bayesian 10-point
native PCA exhaust the budget; Bayesian 50-point matched PCA terminates with a
line-search error. These outcomes are labeled in the figures and CSVs, and should
not be interpreted as converged local optima. All original-noise GP fits converge.

[Matched PCA at three noise levels](matched_pca_three_noise.png)
([PDF](matched_pca_three_noise.pdf));
[principal-curve iteration paths](princurve_iteration_recovery.png)
([PDF](princurve_iteration_recovery.pdf)).

## Interpretation

MPCurve's initialization sensitivity is not unique: GPLVM also remains close to
poor PCA embeddings here. However, the low-noise principal-curve recovery is a
specific counterexample to the idea that this initialization cannot be repaired.
Its smoothing followed by orthogonal projection is worth investigating. The
fixed df = 5 smoother, curve parameterization, projection update and endpoint
stretching all differ from MPCurve; attributing the difference to one of them
requires an ablation. The current data support an algorithmic difference, not a
causal claim about which component produces it.

The Isomap k = 15 comparison similarly shows more repair by GPLVM than MPCurve,
while k = 10 gives near-identical good recovery across methods. Joint fitting
contributes to the original MPCurve deterioration: B-only correlation is 0.759
rather than 0.700, although that still falls below GPLVM's 0.807. This selected
case and correlated noise variants are a diagnostic, not a replicated benchmark
or evidence of general method superiority. A reviewer-facing study should compare
both native starts and common starts, retain convergence information, and separate
known-group recovery from learning multiple orderings.

Supplementary fits of A, C, D and E at original noise using native principal-curve
and classical GPLVM PCA starts are retained in `all_comparisons.csv`. They are not
used to select settings for B. Principal curve with tighter stopping on A reaches
its 1000-iteration limit; default runs on C, D and E also hit their 10-step limit.

## Reproduction and verification

Run from the InferOrder root. R uses the existing study wrapper and frozen package
library. Install `princurve` 2.1.6 into that library if missing. Python packages are
pinned in `requirements.txt`; `.venv/` is ignored and `python_freeze.txt` records
all installed dependencies.

```bash
bash experiments/estimate_intrinsic_m_smooth_v032/run_r.sh experiments/estimate_intrinsic_m_smooth_v032/exploratory_m5_ordering_b/external_methods/prepare_inputs.R
bash experiments/estimate_intrinsic_m_smooth_v032/run_r.sh experiments/estimate_intrinsic_m_smooth_v032/exploratory_m5_ordering_b/external_methods/run_princurve.R
python3.9 -m venv experiments/estimate_intrinsic_m_smooth_v032/exploratory_m5_ordering_b/external_methods/.venv
experiments/estimate_intrinsic_m_smooth_v032/exploratory_m5_ordering_b/external_methods/.venv/bin/pip install -r experiments/estimate_intrinsic_m_smooth_v032/exploratory_m5_ordering_b/external_methods/requirements.txt
bash experiments/estimate_intrinsic_m_smooth_v032/exploratory_m5_ordering_b/external_methods/run_python.sh {1..25}
bash experiments/estimate_intrinsic_m_smooth_v032/run_r.sh experiments/estimate_intrinsic_m_smooth_v032/exploratory_m5_ordering_b/external_methods/run_controls.R
bash experiments/estimate_intrinsic_m_smooth_v032/run_r.sh experiments/estimate_intrinsic_m_smooth_v032/exploratory_m5_ordering_b/external_methods/render_report.R
```

The 69 comparison records comprise 26 principal-curve, 25 GP, nine single-ordering
MPCurve and nine previously saved joint MPCurve fits. Inputs, hashes, actual
initial coordinates, final coordinates, status and parameters are retained.
The first 19 GP results record the initial driver hash, whose source is preserved
as `run_gpy_initial_source.py.txt`; jobs 20–25 use its extension for matched PCA.
GPy emitted NumPy/paramz deprecation warnings, saved in each JSON. R fit warnings
are empty; session reporting emitted an environment timezone warning.

Checks verify 300 finite coordinates per fit, all 25 GP results, matched PCA rank
correlations, and the expected 69 records. All three PNG figures were visually
inspected for panel mappings, labels and convergence flags. Principal-curve
iteration diagnostics recompute squared projection distance directly because
its package-returned initial PCA distance uses a different scaling; both fields
are saved. Original fits, package sources and website sources are preserved.
