# Ordering B in the MPCurver 0.3.4 challenging example

This exploration investigates `main_M5_S4_r001`: 300 samples, 60 features,
five true orderings, variance SNR 4, replicate 1. The saved automatic-M fit
selects five groups and recovers the feature partition exactly (ARI 1), but
recovers B poorly. All analysis outputs stay in this directory.

## Findings

The 12 B features, including monotone anchor V20, are grouped correctly both
initially and finally. True B corresponds to fitted slot A (slot 1); software
slot labels must not be interpreted as ground-truth labels. Absolute Spearman
correlation with true B positions is 0.749564 for raw Isomap, 0.744531 after
the two initialization sweeps, and 0.700139 in the final fit. C is also imperfect
(final 0.847023); A, D, and E are above 0.997. Mean realized B noise variance
is 0.249697, comparable to the other groups (0.248020--0.265553).

With the default 15-nearest-neighbor graph, three of 4,500 directed neighbor
entries connect samples whose true positions differ by more than 0.25. All
three actually span 0.937--0.965 of the unit interval, linking the two ends of
the trajectory. The 0.25 threshold is a truth-aware diagnostic, not an
estimation rule. These edges are listed in `long_range_edges.csv`.

On exactly the same observations, reducing the neighborhood to 10 eliminates
these end-to-end edges and improves raw Isomap recovery to 0.995146; five
neighbors gives 0.994322. At 15 neighbors, halving the same realized noise gives
0.998680, and noiseless signals give 1.000000. Graph-geodesic versus true-distance
Spearman correlation rises from 0.824845 at k=15 to 0.977011 at k=10 on the
observed data. The shape and noise together create shortcut edges that distort
the global one-dimensional embedding. The presence of a monotone anchor does
not guarantee that multivariate Euclidean nearest neighbors respect its order.

## Controlled full-data restarts

Reconstruct the original initial fits with the frozen package and original
seed. Fix M at the already correctly selected five groups. For the intervention,
replace only B's initialization with k=10 Isomap, applying the same two subset
CAVI sweeps and zero-sweep expansion to all features. Other initial fits and
all subsequent fitting controls remain the same. Balanced groups imply the
same initial global weights of 0.2 under this fixed-M restart.

| B initialization | Final B correlation | Final objective | Sweeps | Final ARI |
| --- | ---: | ---: | ---: | ---: |
| Original k=15 | 0.700139 | -17989.571542 | 161 | 1 |
| B-only k=10 | 0.996823 | -17839.191284 | 173 | 1 |

Both fits converge. The original restart reproduces the saved B recovery;
the replacement reaches a better objective by 150.380258. This supports an
initialization-induced poor local solution, rather than insufficient signal
to recover B in this dataset. Choosing B and comparing neighborhoods here is
an exploratory diagnosis of one selected example, not validation of a new
default. In particular, noiseless k=5 produces three disconnected graph
components; only 244 samples receive coordinates, and its full-sample recovery
is recorded as NA.

## Reproduce and inspect

From the InferOrder root:

```sh
bash experiments/estimate_intrinsic_m_smooth_v032/run_r.sh \
  experiments/estimate_intrinsic_m_smooth_v032/exploratory_m5_ordering_b/diagnose.R
```

Use `--diagnostics-only` to omit the two full-data restarts. The frozen 0.3.4
package is commit `15f2b0bbe5dfa61cd46da5160b2bc251e75a0475`; input hash,
package archive provenance, seed, and R session are saved in `provenance.rds`.
The script validates the fixed input hash and correspondence of reconstructed
initial clusters to saved clusters. Its inputs are the existing local dataset,
compact result, and retained full fit. The script emits an expected warning
for the disconnected noiseless k=5 graph.

- `initial_vs_final.png`: all five orderings before and after fitting.
- `geometry_sensitivity.png`: B positions across neighborhood and noise settings.
- `b_features.png`: all B signals and observations, with the anchor labeled.
- CSV files retain the plotted recovery metrics, edge diagnostics, and restarts.

No package code, saved study fit, formal summary, or website source/output is
modified. No commit or publication is part of this exploration.

## Within-panel initialization versus final ordering

`fit_sensitivity.R` extends the neighborhood/noise grid with full-data MPCurve
fits. It changes only the B observations to signal plus the specified multiple
of their original noise. M is fixed at the correctly selected five groups;
the original feature groups initialize the fit. Only B's Isomap neighborhood
size changes. All groups' initial responsibilities are expanded to the modified
full matrix with zero CAVI sweeps before joint structural fitting. Subsequent
feature assignments, trajectories, noise, and adaptive priors are learned.
The original annealing, tolerance, and 10,000-sweep study limit are retained.
Thus this is a controlled initialization study, not a new automatic-M benchmark.

Each setting is cached in `sensitivity_fits/setting_XX.rds`, retaining the input
hashes, seed, package provenance, initial/final coordinates, assignments,
objective trace, convergence status, and warnings. Settings are indexed with
noise multipliers 0, 0.5, 1 varying fastest within k=5, 10, 15, 20, 30.
The disconnected noiseless k=5 setting is explicitly left without a full-data
fit, rather than assigning invented coordinates to the missing samples.

```sh
bash experiments/estimate_intrinsic_m_smooth_v032/run_r.sh \
  experiments/estimate_intrinsic_m_smooth_v032/exploratory_m5_ordering_b/fit_sensitivity.R \
  1 2 3 4 5 6 7 8 9 10 11 12 13 14 15
bash experiments/estimate_intrinsic_m_smooth_v032/run_r.sh \
  experiments/estimate_intrinsic_m_smooth_v032/exploratory_m5_ordering_b/render_sensitivity_comparison.R
```

The resulting `geometry_initial_final_comparison.png` and PDF overlay raw
Isomap and final posterior-mean sample ranks in each panel. Both use normalized
average ranks, so coordinate stretching is not mistaken for order changes.
Isomap is oriented to truth and the final fit to Isomap; a global reversal
is equivalent. Gray segments connect the same sample at both stages. Panel
annotations give absolute correlation to truth at both stages and the
initial-to-final correlation. The earlier raw-coordinate sensitivity figure
is retained separately. The renderer verifies that observed-data k=10 and k=15
reproduce the earlier controlled restart objectives and B recovery.

All 14 connected settings converged without warnings and retained ARI 1.
At original noise, k=15/20/30 have final truth correlations
0.700139/0.688720/0.676222, with initial-final correlations
0.970361/0.972178/0.963130. The iterations adjust these folded solutions but
preserve the main ordering error. At k=5/10 final recovery is
0.996793/0.996823. The original-noise k=10/15 control checks pass.

## PCA comparison at the same three noise levels

The PCA extension changes only B's initialization to the package's PC1 ordering
on its 12 features, using the default centering and no variance scaling.
It retains the same two subset CAVI sweeps, full-data expansion, five initial
groups, and subsequent structural fitting controls. Other groups keep their
original Isomap initialization. This isolates the B initialization method;
it is not a fit initialized by PCA in every group.

```sh
bash experiments/estimate_intrinsic_m_smooth_v032/run_r.sh \
  experiments/estimate_intrinsic_m_smooth_v032/exploratory_m5_ordering_b/fit_sensitivity.R --pca 1 2 3
bash experiments/estimate_intrinsic_m_smooth_v032/run_r.sh \
  experiments/estimate_intrinsic_m_smooth_v032/exploratory_m5_ordering_b/render_pca_comparison.R
```

All three fits converge without warnings and retain feature-partition ARI 1.
Absolute Spearman recovery of B is:

| B noise multiplier | Initial PCA | Final MPCurve | Initial-final correlation |
| --- | ---: | ---: | ---: |
| 0 | 0.103159 | 0.122720 | 0.982568 |
| 0.5 | 0.123075 | 0.153552 | 0.984501 |
| 1 | 0.141566 | 0.145150 | 0.968001 |

PC1 is nonmonotone in the true B position even without noise: it folds the
trajectory into a U-shaped ordering. The subsequent iterations largely retain
that fold. On these same inputs, Isomap k=10 gives final recovery
0.999867, 0.999089, and 0.996823. The PCA objective is lower than the matched
Isomap k=10 objective at all three noise levels. Objectives should only be
compared within a noise level, since changing noise changes the observations.

`pca_initial_final_comparison.png`/PDF shows the three PCA panels;
`pca_vs_isomap_comparison.png`/PDF adds matched Isomap k=10 and k=15 rows.
`pca_vs_isomap_summary.csv` retains recovery, convergence, objectives, and ARI.
The renderer verifies identical modified-input hashes across PCA and Isomap
within each noise level. Cached PCA records are under `pca_fits/`.
