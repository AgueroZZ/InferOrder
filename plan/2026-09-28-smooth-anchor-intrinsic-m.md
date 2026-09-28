# Paired ordering-count experiment with smooth trajectories

## Scientific question

Measure how estimating the number of sample orderings changes when each
ordering retains one monotone trajectory and its remaining trajectories have
smooth peaks and troughs.

## Design

Use the original 90 datasets: N = 300, D = 60, true M = 3, 4, 5, variance
SNR = 1, 4, 16, and ten repeats per condition. Preserve latent positions,
feature assignments and order, realized measurement noise, and the first
feature in each ordering as a monotone anchor. Replace other trajectories
with independent random sine/cosine mixtures at frequencies k*pi, k = 2, 3, 4,
with coefficient weights k^-2. Center and standardize sampled signals to
variance one. Store shape seeds, coefficients, normalization constants,
derivative certificates, input hashes, and paired baseline provenance.

## Fitting and evaluation

Freeze the same MPCurver 0.3.2 archive and retain the baseline priors,
initialization, annealing, normalized stopping tolerance of 1e-6, maximum
M = 8, and 10,000-sweep candidate cap. Compare adaptive EB and fixed-uniform
forward selection on each input. Keep incomplete candidates resumable and
retain errors or nonconvergence as scientific outcomes.

For both methods, count posterior assignment columns with mean above 1e-12
as effective M. Evaluate feature groups with ARI and occupied sample orderings
with one-to-one absolute Spearman matching; assign zero to unmatched true
orderings. Reevaluate saved baseline fits with this same matching rule.

## Reporting

Use a separate workflowr page to compare paired count recovery and structural
recovery. Show saved replicate-one trajectories for all nine conditions and
preselected SNR = 4 examples. Include a descriptive initialization-similarity
diagnostic and report the modified study's runtime. Preserve the original
benchmark and its version provenance.
