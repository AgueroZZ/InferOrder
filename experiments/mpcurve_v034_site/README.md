# MPCurver 0.3.4 results for the InferOrder website

This update regenerates the introductory simulations, pancreas analysis, and
fitness summary with MPCurver 0.3.4 from package commit
`15f2b0bbe5dfa61cd46da5160b2bc251e75a0475`. The frozen source archive is
`source/MPCurver_0.3.4.tar.gz`, with SHA-256
`58d0c99c6170994c82eedba190fb1c28a63eb47530c575705323c3cf3992b7c6`.

Run all commands from the InferOrder repository root with that archive installed:

```sh
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 \
  Rscript --vanilla experiments/mpcurve_v034_site/run_simulations.R

OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 \
  Rscript --vanilla experiments/mpcurve_v034_site/run_pancreas.R

python3 experiments/mpcurve_v034_site/extract_fitness.py \
  --mpcurver-root ../MPCurver
```

`run_simulations.R` reproduces the single- and two-ordering examples from the
public MPCurver introduction and feature-partitioning vignettes. The generator
parameters, random seeds, fitting calls, saved fits, tabular summaries, figures,
and R session are recorded under `results/simulations/`.
The two-ordering simulation uses Fiedler rather than PCA initialization because
the checked PCA fit separated the feature groups but did not recover the second
sample ordering on this data set. The published site reports the specified
Fiedler fit, allows up to 500 temperature-one refinement updates, and
identifies the result as one illustrative simulation.

`run_pancreas.R` uses the tracked `data/loading_order/pancreas_factors.RData`
and `data/loading_order/pancreas.RData` files. It reproduces the original
semi-NMF loading selection: up to 500 cells per cell type with seed 1, factor
columns 3, 8, 9, 12, 17, 18, 20, and 21, followed by the `inDrop3` filter.
The same finite 865-by-8 matrix is passed to fixed-M single- and two-ordering
CAVI fits. The M = 1 fit uses PCA initialization, `K = 50`, RW2, `ridge = 0`,
and up to 150 iterations. The M = 2 fit uses Fiedler initialization within
similarity-initialized feature groups, `K = 50`, RW2, `ridge = 0`, 25 annealing
sweeps, and up to 500 temperature-one refinement sweeps. The reported M = 2
fit assigns seven factors to one ordering and factor 20 to the other. No
automatic dimension selection is performed. If a fit completed but plotting or table
generation was interrupted, `run_pancreas.R --resume` reuses saved fits only
after checking that their stored input matrices equal the reconstructed input.

## Regenerated results

Both simulation fits converged under the normalized per-observation rule. The
single-ordering absolute Spearman correlation was 0.9997. The two-ordering fit
exactly recovered all 20 feature assignments, with ordering correlations
0.9975 and 0.9989; it used 369 recorded temperature-one objective values.

Both pancreas fits also converged. The two-ordering fit used 162 recorded
temperature-one objective values, retained the seven-versus-one feature split,
and placed factor 20 alone on ordering B. The absolute correlations between the
single-ordering positions and orderings A and B were 0.9679 and 0.0576.

The extracted fitness results assign 27 environments to ordering A and 18 to
B. The two mutant orderings have correlation 0.857, and the two single-ordering
fits under different smoothness settings have correlation 0.952.

The website pages load these saved results and do not rerun statistical fits
during a workflowr build. The fitness page is a concise derivative of the
MPCurver 0.3.4 public fitness vignette at the package commit recorded above.
Its UTF-8-corrected generated article comes from documentation-only commit
`eeaf65bcf4c8b021d590c88ad411506359202f1f`.
`extract_fitness.py` copies three figures and extracts the 45-environment
assignment table and both reported ordering correlations from that article. The
SHA-256 source manifest is in `data/site_assets/fitness/source_manifest.csv`.
The complete vignette remains the canonical analysis.

The exact article and figure hashes are regenerated in
`data/site_assets/fitness/source_manifest.csv` by `extract_fitness.py`.
