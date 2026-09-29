# Estimate-M simulation structure and MPCurver 0.3.4 refresh

- Agent: Codex
- Date: 2026-09-28 UTC
- Request: organize the workflowr simulation section into three study classes,
  update both estimate-M designs with the current package, and publish the
  refreshed results.

## Site structure

The homepage now presents three simulation classes: fixed M = 1, fixed M = 2,
and estimating M from the data. The estimate-M class contains two explicit
substudies: the simpler all-monotone design and the harder design with one
monotone anchor per ordering and nonmonotone trajectories for all other
features. The Simulation navigation labels distinguish these two estimate-M
pages as well.

## Current automatic-M fits

Both estimate-M extensions use the frozen MPCurver 0.3.4 archive from commit
`15f2b0bbe5dfa61cd46da5160b2bc251e75a0475`, SHA-256
`58d0c99c6170994c82eedba190fb1c28a63eb47530c575705323c3cf3992b7c6`.
They fit the original fixed 90-dataset inputs with `intrinsic_dim = "auto"`,
five-df natural-cubic-spline variance explained, single linkage, maximum M =
8, minimum cluster size two, cluster-size initial adaptive probabilities, and
Isomap initialized independently within each selected feature group.

- All-monotone Slurm array 59675733 completed 90/90 tasks with exit code zero.
  All fits converged without warnings. The initial similarity cut and final
  adaptive-EB effective M were both exact in 90/90 datasets. Mean final ARI
  was 1.0000 and mean ordering recovery was 0.9849. Median runtime was 95.0
  seconds, compared with 150.6 seconds for the original adaptive M = 8 fit and
  222.7 seconds for uniform + forward. The paired controls recovered M in
  87/90 and 90/90 datasets, respectively.
- One-anchor/nonmonotone Slurm array 59675746 completed 90/90 tasks with exit
  code zero. All fits converged without warnings. The initial similarity cut
  and final effective M were both exact in 88/90 datasets. Mean final ARI was
  0.9945 and mean ordering recovery was 0.8390. Median runtime was 34.6
  seconds, compared with 66.7 seconds for the original adaptive M = 8 fit and
  132.4 seconds for uniform + forward. The paired controls recovered M in
  60/90 and 88/90 datasets, respectively.

Strict validation reloaded every compact fit, checked fixed-input hashes and
package provenance, required positive finite noise estimates and normalized
assignment weights, confirmed convergence and zero warnings, and verified
effective-count stability at occupancy thresholds from 1e-12 to 1e-3. Nine
retained full fits per design confirmed that every block requested and used
Isomap independently and did not use a PCA component.

## Reporting and verification

The all-monotone and one-anchor report sources now use the current 0.3.4
automatic-M summaries, figures, runtime comparisons, and validation records.
The original MPCurver 0.3.2 adaptive and forward fits remain explicitly
labeled paired controls. The all-monotone page title now identifies its
trajectory design.

All eight workflowr pages were rebuilt. The curated site checker validates the
three-category homepage hierarchy, page titles and headings, navigation and
local links, embedded-image counts, package commit and archive checksum,
90-dataset completion, reported exact-recovery counts, minimum cluster size,
independent per-group Isomap initialization, and saved-result claims. The
checker and `git diff --check` pass. The principal recovery, count-distribution,
and structural-recovery figures were visually inspected for labels, clipping,
and consistency with the saved tables.

InferOrder commit `99deed26bbd9cead2b96e036d3c79648c511656f` was pushed to
`origin/master`. GitHub Pages deployment
[36512103474](https://github.com/AgueroZZ/InferOrder/actions/runs/36512103474)
completed successfully for that exact commit. The live homepage and both
estimate-M pages are byte-identical to their committed HTML files.
