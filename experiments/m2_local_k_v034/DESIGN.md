# Fixed-M=2 one-step neighborhood-selection study

Prespecified before fitting: 36 datasets, from two smooth random Fourier
families crossed with SNR 1/4/16 and six replicates. N=300, D=24, two balanced
12-feature groups, independent uniform latent positions. Every feature is a
random sine/cosine mixture; there are no imposed monotone anchors, monotonicity
filters, or outcome-based rejection/resampling. The broad family uses pi times
frequencies 2:4 with k^-2 attenuation; rich uses 2:6 with k^-1.5 attenuation.
Each signal feature has sample variance one. Independent Gaussian residuals
have variance 1/SNR. Signal, noise realization, feature permutation, and latent
positions are paired across SNR within each family/replicate (12 independent
base realizations, not 36 independent realizations).

Primary regime: only M=2 is known. Infer a two-group feature partition with the
frozen package's spline-R-squared df=5 similarity and single linkage. Both
methods receive this same partition. Diagnostic regime: initialize with true
feature groups, to isolate position-initialization effects. Subsequent joint
fits in both regimes learn feature assignments and adaptive priors normally.
True groups are not imposed during fitting.

Baseline uses Isomap k=15 and the package's ordinary two-update local warmup.
The selector considers k=5,10,15,20,30 separately within each initialized group,
scores each connected full-sample candidate after exactly one single-ordering
CAVI update, and selects maximum local ELBO (smallest k breaks exact ties).
Then run the same two-update warmup from the chosen raw positions and the same
full joint fitting controls as baseline. This preserves the ordinary warmup
and isolates the selection decision. The capital K grid size is fixed at 50.
No truth, recovery metric, or converged candidate fit is used for selection.

Disconnected candidates are excluded; the package's usual Isomap/PCA fallback
is retained if default k=15 cannot initialize all samples or no candidate is
eligible. Singleton groups use the package's direct feature-rank initialization
and have no meaningful neighborhood selection. Log these cases separately.

All generated datasets are retained. Report default failures (matched |rho|<0.9),
rescues to >=0.95, changes in ordering recovery, regressions >0.05, group ARI,
convergence, warning/fallback counts, and all timings. Also report candidate
coverage using truth only as a post-hoc diagnostic: whether any candidate's
one-step ordering reaches >=0.95. Counts are exploratory, conditional on the
specified families; no outcome-dependent enrichment or retry of seeds.

Timing includes similarity/clustering in the inferred regime, embedding each
required k, local scoring, the retained warmup and full-matrix expansion, and
joint fitting. Shared candidate computations may be reused in the diagnostic
execution; report pipeline time by summing separately measured required
components rather than treating the entire experimental driver time as either
method's deployment time. Give screening overhead and joint fitting time
separately. Wall times from concurrent single-thread workers are descriptive,
not a hardware-independent benchmark. Compare against the actual converged
M=2 outcomes, not only the initial ELBO.

Freeze MPCurver 0.3.4 at commit 15f2b0bbe5dfa61cd46da5160b2bc251e75a0475.
RW2, adaptive per-feature noise and smoothing, temperature 5->1 over 25 sweeps,
normalized stopping tolerance 1e-6, maximum 10,000 sweeps. Inputs are samples by
features. Keep source/package/input hashes, all generator coefficients, seeds,
candidate scores, selected k, fitted coordinates, assignment weights, and
objective histories. No package or website changes are included in this study.
