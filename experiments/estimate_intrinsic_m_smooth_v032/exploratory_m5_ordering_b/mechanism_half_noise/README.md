# Why principal-curve iteration repairs the half-noise PCA fold

The half-noise B case separates the smoothness penalty from the complete
ordering algorithm. Replacing principal curves' cubic smoothing spline with a
50-node RW2 smoother, calibrated to five effective degrees of freedom, retains
near-perfect recovery (absolute Spearman correlation 0.998817). Increasing the
cubic spline's degrees of freedom to 10 destroys recovery (0.118739). Yet simply
making MPCurve smoother does not repair its ordering. Thus the penalty family
alone does not explain the gap: smoothness, latent-coordinate parameterization,
and the curve/position update scheme interact.

## Fixed design and source

All diagnostic fits use the same 300-by-12 matrix `B_noise05_X.csv` from
`../external_methods/inputs/`, the same default unscaled PCA ordering, and one
ordering. This is one fixed signal with half the original noise realization,
with nominal noise variance 0.0625. Ground truth is used for evaluation only,
except the explicitly labeled known-noise oracle control. No experiment selects
an endpoint using truth. These settings were investigated after seeing the
benchmark; they are mechanism diagnostics, not a new prespecified benchmark.

MPCurver 0.3.4 is frozen at source commit
`15f2b0bbe5dfa61cd46da5160b2bc251e75a0475`. The installed CAVI and initialization
functions are saved in `frozen_cavi_source.R.txt`. Principal curves use princurve
2.1.6 with default PCA, `maxit=1000`, `thresh=1e-6`; smoother df or other settings
change only when specified. MPCurve uses the public single-ordering PCA API,
2000 maximum sweeps and the default 1e-6 convergence tolerance. Seeds are fixed
at 20260929 for the MPCurve controls. Both baseline endpoints reproduce the
previous benchmark.

## Similar smoothing penalties do not imply the same iteration

An RW2 prior penalizes squared second differences on a regular grid. Up to grid
scaling, this approximates the integrated squared second derivative used by
cubic smoothing splines; boundary treatment and discrete representation differ.
See [Lindgren and Rue's RW2 note](https://www.math.ntnu.no/preprint/statistics/2005/S6-2005.pdf)
for the second-difference and continuous-process connection. The precise
comparison is conditional on the same input coordinates and smoothing strength.
Ordinary PCA restricts the trajectory to a straight line; MPCurve is more
naturally described as a probabilistic nonlinear principal-curve model.

For fixed responsibilities and noise, MPCurve's feature-j posterior mean is

```
m_j = [diag(N_k / sigma_j^2) + lambda_j Q]^{-1} R' x_j / sigma_j^2.
```

This is a penalized smoother. Its conditional smoothing degrees of freedom are
`trace([diag(N_k / sigma_j^2) + lambda_j Q]^{-1} diag(N_k / sigma_j^2))`.
This trace measures the quadratic expected-complete-data smoothing subproblem
for weighted bin pseudo-observations. It conditions on responsibilities and fitted
hyperparameters; it is not the observed-data hat-matrix trace with soft assignments,
nor the total degrees of freedom of the nonlinear estimator learning positions.
In this case the initial conditional df are 35.51–47.14 and the converged values
16.63–22.04. Principal curves hold each feature's smoother near df=5. Identical
PCA ranks therefore do not imply similarly regularized initial trajectories.

Principal curves smooth feature values against current projected arc length,
then project each observation onto a continuous polygonal curve in Euclidean
feature space. Projected arc length supplies the next iteration's coordinates.
This is the alternating construction in
[Hastie and Stuetzle (1989)](https://hastie.su.domains/Papers/Principal_Curves.pdf),
and the installed implementation is saved in the parent external-method report.

MPCurve instead computes discrete variational position probabilities,

```
log r_ik = constant_i + log pi_k
           - 0.5 sum_j [(x_ij - m_jk)^2 + Var(U_jk)] / sigma_j^2.
```

It then updates position weights, per-feature noise, per-feature smoothness and
trajectory posteriors. Its coordinate is a fixed latent grid, with no arc-length
reparameterization step. Continuous projection, feature weighting, uncertainty,
position weights and coordinate geometry therefore differ even with the same
penalty family. The responsibilities have support over every grid point; there
is no restriction to adjacent moves or explicit penalty for leaving the previous
sample position. A poor curve/position configuration can nevertheless reinforce
itself under these coupled updates.

## Controlled results

| Diagnostic | Final absolute Spearman correlation | What changes |
| --- | ---: | --- |
| Native MPCurve | 0.097072 | Baseline |
| Native P-curve, cubic df=5 | 0.998820 | Baseline |
| P-curve, cubic df=3 | 0.280611 | Less flexible smoother |
| P-curve, cubic df=10 | 0.118739 | More flexible smoother |
| P-curve, cubic df=20 | 0.113660* | More flexible smoother |
| P-curve, RW2 df=5 | 0.998817 | Replace cubic smoother with RW2; retain projection/arc-length iteration |
| P-curve, rank-spaced cubic df=5 | 0.448493 | Smooth against average ranks instead of arc-length spacing |
| P-curve, no endpoint stretching | 0.998820 | Set stretch=0 |
| MPCurve, initial conditional df=5, then adaptive lambda | 0.070582 | Smoother initial trajectories |
| MPCurve, initial df=5 calibration, lambda subsequently fixed | 0.124982 | Final conditional df=3.96–6.77 |
| MPCurve, fixed lambda=10000 | 0.060886 | Final conditional df=5.05–5.91 |
| MPCurve, fixed lambda=100000 | 0.348918 | Final conditional df=2.88–3.50 |
| MPCurve, known true noise variance | 0.187227 | Supply S=0.25, removing noise estimation and imposing equal noise weights |
| MPCurve, fixed uniform position prior | 0.077229 | Remove adaptive location weights |
| MPCurve, equal-width initial bins | 0.091316 | Change quantile bin initialization only |

The starred df=20 run reaches its 1000-iteration limit; every other diagnostic
in the table reports convergence. All fit warning lists are empty. The additional
fixed lambda=100 and 1000 controls are retained in `ablation_summary.csv` and
also remain poor (0.085044 and 0.097649). Greater smoothness is not uniformly
better: the principal-curve df=3 control also fails to unfold correctly.

The RW2 hybrid is an experimental principal-curve algorithm, not an MPCurve fit.
At each smoothing step, observations are linearly interpolated on 50 regular
latent nodes with the exact Q from the frozen MPCurve fit. Its penalized
least-squares operator is calibrated to trace=5 using current coordinates only.
The rest of `principal_curve` is unchanged. Verification checks df=5 to 1e-6 and
reproduction of constant/linear functions to 1e-7 (observed errors below 1.2e-12).
This shows that a 50-node RW2 smoothing representation can support recovery;
it does not isolate continuous projection from all other differences with CAVI.

The rank-spacing hybrid changes the coordinates supplied to the spline smoother
at every iteration. Projection still returns arc lengths. Its result supports
sensitivity to parameterization, but is not an exact reproduction of MPCurve's
fixed-grid latent model and does not by itself identify that model's failure cause.

[Initial/final diagnostic figure](mechanism_comparison.png)
([PDF](mechanism_comparison.pdf)). All six displayed fits converge; curves share
the initial PCA ordering and use average normalized ranks for comparison.

## What the evidence establishes

The intuition that RW2 and cubic splines are closely related is supported by
both the smoothing equations and the successful RW2 principal-curve hybrid.
The native implementations differ substantially in their effective smoothing,
parameterization and position updates. The principal-curve df=10 failure shows
that its success is sensitive to regularization. The failed strongly smoothed
MPCurve controls show that regularization alone does not reproduce that success.

At the native MPCurve endpoint, mean maximum position probability is 0.976 and
23 of 50 bins are occupied by MAP positions. Posterior averaging alone is thus
not a sufficient explanation. Estimated noise variances are approximately
0.055–0.078, around the true 0.0625; known-noise fitting does not repair the fold.
Neither fixing location weights nor removing endpoint stretching explains the
observed method gap by itself. These controls narrow the explanation but do not
prove a unique causal component or demonstrate a general fix for MPCurve.

A useful working hypothesis is that sufficiently smooth curve updates with
continuous projection and arc-length reparameterization can escape this folded
PCA configuration, whereas the fitted discrete variational model settles into a
configuration that already explains the observations while preserving the wrong
global ordering. Establishing which individual change repairs CAVI requires
further controlled modifications of its position and coordinate updates.

## Reproduction

From the InferOrder root, use the existing study R wrapper for each script in
this directory, in order: `inspect.R`, `run_ablations.R`, `run_hybrids.R`,
`render_report.R`, `verify_hybrid.R`. The first script saves native fit prefixes
and feature-specific noise/smoothness diagnostics. The next two save all fits,
convergence flags and warnings. Verification checks baseline reproduction,
nominal noise variance and the RW2 smoother operator. Source hashes are in
`source_hashes.csv`. The rendered figure was visually inspected. Package and
website sources remain unchanged.
