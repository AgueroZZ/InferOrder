# MPCurver 0.3.3 automatic-M simulation update

## Question and design

The 90 fixed one-monotone-anchor datasets were refitted with the released
MPCurver 0.3.3 automatic ordering-count initializer. The extension changes
only the initialization: fixed-df spline variance explained with df = 5,
single linkage, the maximum-mean-silhouette eligible cut up to M = 8, minimum
cluster size two, and adaptive initial global probabilities equal to cluster
size divided by 60. The original adaptive M = 8 and uniform-forward fits are
unchanged paired controls.

The extension freezes package commit
`c905901424e43eab78b58bdcc0d1de367ec8fd73` and source archive SHA-256
`cb57e1f6e8d859b7ff9aeda17c97d538fae0686eb0ceb616f126023a26d645f9`.
Slurm array 59668185 completed all 90 one-thread tasks with exit status zero;
no stderr file was nonempty.

## Results

- The similarity cut estimated the true M in 88/90 datasets, with two
  overestimates and no underestimates. The mean initial partition ARI was
  0.9986, and every initial cluster contained at least two features.
- Automatic-M adaptive EB finished at the true effective M in 88/90 datasets,
  with two underestimates and no overestimates. Relative to original adaptive
  M = 8, it gained 28 exact results and lost none.
- Mean final feature-partition ARI increased from 0.7300 to 0.9945; mean
  matched ordering recovery increased from 0.6847 to 0.8390.
- Uniform + forward was also exact in 88/90, with mean ARI 0.9949 and ordering
  recovery 0.8432. Median runtime was 32.5 seconds for automatic-M adaptive,
  66.7 seconds for original adaptive, and 132.4 seconds for forward selection.
- The two initial and two final errors are not identical. Adaptive EB corrected
  one initial M = 5 overestimate when true M = 4, collapsed one correct M = 5
  initialization to four, and reduced one initial M = 7 overestimate to four
  when true M = 5.

## Validation and website

Strict validation passed for all 90 extension results: every fit converged,
there were no warnings, package/source/input hashes matched, and effective-M
counts were unchanged across occupancy thresholds from 1e-12 through 1e-3.
The experiment retains compact per-dataset results, source archive, pilot,
summary tables and figures, scheduler metadata, and reproduction scripts under
`experiments/estimate_intrinsic_m_smooth_v032/auto_m_v033/`.

The workflowr page and homepage were rebuilt with the repository's curated
builder. The audit passed all eight public pages, local links, 13 embedded
figures on the updated page, terminology, and saved-result claims; no absolute
local image paths remain. The key automatic-M accuracy, count-distribution,
and structural-recovery figures were visually inspected. Publication commit
and Pages verification are appended after push.

## Publication verification

The simulation and website were committed as
`a6003b57eb25997da26248696f055424ae5015fb` and pushed to `origin/master`.
GitHub Pages deployment
[36492621622](https://github.com/AgueroZZ/InferOrder/actions/runs/36492621622)
completed successfully for that commit. The live result page and homepage are
byte-identical to their committed HTML files; they contain the MPCurver 0.3.3
method, 88/90 and 60/90 comparison, and runtime results.
