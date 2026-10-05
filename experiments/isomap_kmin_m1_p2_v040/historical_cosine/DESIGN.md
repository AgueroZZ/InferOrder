# Paired Isomap neighborhood comparison: M=1, P=2

Prespecified on 2026-10-02, before generating or fitting these datasets.

## Question and comparison

Compare MPCurve position recovery from automatic Isomap `k_min` with the former
fixed k=15 default. Fixed k=10 is a secondary reference. All three settings use
the same frozen MPCurver 0.4.0.9000 implementation, observations, fit seed,
position grid, model, and convergence rule within each replicate. Only the
neighbor-count argument differs. No ELBO screening or truth-based initialization
selection is performed. M=1 denotes one latent ordering; P=2 denotes two
observed signal features, with no noise-only features or feature grouping.

## Generator and independent repetitions

- Twenty independent replicates, N=200 samples each. Sample positions are
  independent Uniform(0,1), stored in their original random row order.
- Two random cubic B-spline trajectories share those positions. Each trajectory
  has eight independent N(0,1) coefficients; boundary knots are 0 and 1 and
  interior knots are 0.2, 0.4, 0.6, 0.8. This adapts the existing
  `m1_bspline_comparison` generator from twelve features to two.
- Both trajectories must be nonmonotone: on a 2,001-point grid, each has at
  least one increment >1e-8 and one <-1e-8 after scaling. Coefficient pairs
  that fail this generating-condition check are redrawn before adding noise
  or fitting; attempts are saved. Curve intersections and difficult shapes
  are retained. Fitting outcomes never determine which datasets enter.
- Center each trajectory using its dense-grid mean, then use a single common
  scale so the average dense-grid centered signal variance across features is
  one. Relative feature amplitudes are preserved. Independent Gaussian errors
  have SD 0.25, corresponding to average variance SNR 16 under this definition.
- Replicate generating seeds are 202610020 plus replicate index. Fit seeds are
  202620020 plus replicate index. A separate seed 202630020 randomizes the
  order of the three fits within each replicate. No observed-variance feature
  standardization is performed.

## Fitting and failure accounting

Use the public `fit_mpcurve()` and `do_mpcurve()` API. All fits have one
ordering, 50 position bins, quantile initialization, RW2, ridge zero, initial
precision one, estimated feature noise/precision, and adaptive position weights.
Use normalized ELBO tolerance 1e-6. The first call allows 2,000 sweeps; continue
in blocks of up to 2,000, with a total budget of 10,000. Preserve the actual
state when continuing. The default PCA fallback policy and default Isomap
largest-component handling are retained and recorded. All landmarks are used
because N=200 is below the default landmark cap of 1,000.

Save every fit's status, warnings, actual initializer, realized k, positions,
ELBO history, iteration count, and input/source hashes. Retain finite
budget-limited endpoints with their convergence flags. Report failed fits and
missing metrics explicitly; do not replace replicates or silently restrict to
successful pairs. Timing is descriptive fitting time, not a speed benchmark.

## Position metric and requested figure

The provisional primary metric is ordinary cosine similarity on the true and
posterior-mean estimated positions in [0,1], allowing a global reflection:

`max(cosine(t, q), cosine(t, 1-q))`, where
`cosine(a, b) = sum(a*b) / sqrt(sum(a^2)*sum(b^2))`.

This retains the user's requested cosine metric. A pending clarification asks
whether centered cosine is preferred; that changes evaluation/reporting only,
not simulation or fitting. Save ordinary cosine, absolute centered cosine
(equivalent to absolute Pearson correlation), and absolute Spearman correlation
for every fit. Truth is used only in data generation and evaluation. No
rank transform or monotone remapping is applied to the primary metric.
Also save the constant-position cosine reference to identify the high baseline
of uncentered cosine; it is not an additional fitted method.

The primary boxplot shows the twenty final endpoint scores per setting, all
individual replicates, and the pairing across methods. Boxes show quartiles and
median; whiskers extend to the most extreme values within 1.5 IQR. Report paired
auto-minus-fixed differences, including their Monte Carlo SE (replicate SD
divided by sqrt(number of complete pairs)). No population-wide superiority
claim or multiplicity-adjusted inference is intended from this small study.

## Provenance and scope

Archive the current package source before installation into an ignored local
library. Record package version, base commits plus runtime source hashes,
archive SHA256, generator settings, seeds, complete inputs, script hashes,
and R session information. Preserve compact fits and plotted data in normal
experiment paths; full reconstructible fits and operational logs are ignored.
This is an internal experiment. Existing package sources, historical results,
and public workflowr pages are not edited or published.
