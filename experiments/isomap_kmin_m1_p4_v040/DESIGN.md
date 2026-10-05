# Paired addition of trajectories: P=2 versus P=4

Prespecified on 2026-10-03 before generating added trajectories or fitting.

## Question and pairing

Test whether adding two independently generated nonmonotone signal features
improves position recovery when a two-feature curve has intersections or close
approaches. Preserve all twenty inputs from `../isomap_kmin_m1_p2_v040/`,
including their true Uniform(0,1) positions, first two trajectories, and exact
observation noise. Append two new features; compare P=4 with the saved P=2
fits within each replicate. No parent dataset or fit is regenerated or replaced.

Every dataset has N=200 and M=1. Fit automatic kmin, fixed k=15, and fixed k=10
using the same frozen MPCurver 0.4.0.9000 source and controls as the P=2 study.
The automatic count is recomputed using all four observed features. There is
no early ELBO selection or truth-based choice of initialization.

## Added features and scaling

For each replicate, use seed 202640020 plus replicate index to draw two cubic
B-spline coefficient vectors with eight independent N(0,1) coefficients each.
Use the same basis and 2,001-point grid as the parent study. Both added features
must have positive and negative increments exceeding 1e-8 in magnitude after
scaling; redraw only if this generator condition fails, before adding noise or
fitting. Retain all intersections and difficult configurations.

Center the two new trajectories on the dense grid and use a common scale to
set their average dense-grid variance to one. The first two features retain
their original scale. Thus average variance across all four features is also
one. Preserve relative amplitudes within each pair. Add independent Gaussian
noise with SD 0.25 to the new features, retaining average variance SNR=16.
Do not standardize observed feature variances. Coefficients, pair-specific
scales, signals, standard-normal errors, and final observations are saved.

## Fits and endpoint metric

Use the same fitting seeds (202620020 plus replicate) and randomized method
order seed (202630020) as the parent study. Controls are 50 position bins,
quantile initialization, RW2, ridge zero, initial precision one, adaptive
position weights, estimated noise/precision, and normalized ELBO tolerance
1e-6. Allow 2,000 sweeps per block and continue the actual state up to 10,000
sweeps. Retain default component/PCA fallback policies and record their use.

Evaluate absolute Spearman correlation between true positions and posterior
means, allowing global ordering reversal and average ranks for ties. Also
record raw Isomap Spearman before discretization or optimization. Truth is
used only in generation and evaluation. Retain all finite endpoints, record
failures and convergence flags, and include all twenty replicates in summaries.

Primary comparisons are P=4 minus P=2 final scores within each method; report
paired means, medians, and Monte Carlo SE (SD/sqrt(number of complete pairs)).
Within P=4, report auto-minus-fixed comparisons. These are descriptive findings
for this generator, with no superiority or equivalence claim across settings.

## Figures and provenance

Show endpoint boxplots for P=2 and P=4 with paired points for each method.
Inspect the previously selected replicate 6 without selecting a new case:
hold feature 1 and feature 2 on the axes while coloring by truth or raw
Isomap positions from all two/four-feature fits. Show all six feature-pair
projections for P=4, colored by truth. Two-dimensional projection crossings
need not correspond to the same pair of positions across all features; the
relevant ambiguity concerns the joint four-feature observations.

Archive the exact package source and record source/script/input/parent hashes,
commits, seeds, and sessions. Preserve compact fits and all plotted numeric
data. Libraries, full reconstructible fits, and operational logs are ignored.
This remains internal: no package default or public workflowr page changes.
