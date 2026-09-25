# MPCurver 0.3.0 results for the InferOrder website

Run all commands from the InferOrder repository root with MPCurver 0.3.0 installed:

```sh
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 \
  Rscript --vanilla experiments/mpcurve_v030_site/run_simulations.R

OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 \
  Rscript --vanilla experiments/mpcurve_v030_site/run_pancreas.R

python3 experiments/mpcurve_v030_site/extract_fitness.py \
  --mpcurver-root ../MPCurver
```

`run_simulations.R` reproduces the single- and two-ordering examples from the
public MPCurver introduction and feature-partitioning vignettes. The generator
parameters, random seeds, fitting calls, saved fits, tabular summaries, figures,
and R session are recorded under `results/simulations/`.
The two-ordering simulation uses Fiedler rather than PCA initialization because
the checked PCA fit separated the feature groups but did not recover the second
sample ordering on this data set. The published site reports the specified
Fiedler fit and identifies the result as one illustrative simulation.

`run_pancreas.R` uses the tracked `data/loading_order/pancreas_factors.RData`
and `data/loading_order/pancreas.RData` files. It reproduces the original
semi-NMF loading selection: up to 500 cells per cell type with seed 1, factor
columns 3, 8, 9, 12, 17, 18, 20, and 21, followed by the `inDrop3` filter.
The same finite 865-by-8 matrix is passed to fixed-M single- and two-ordering
CAVI fits. The M = 1 fit uses PCA initialization, `K = 50`, RW2, `ridge = 0`,
and up to 150 iterations. The M = 2 fit uses Fiedler initialization within
similarity-initialized feature groups, `K = 50`, RW2, `ridge = 0`, 25 annealing
sweeps, and up to 100 temperature-one refinement sweeps. The reported M = 2
fit assigns seven factors to one ordering and factor 20 to the other. No
automatic dimension selection is performed. If a fit completed but plotting or table
generation was interrupted, `run_pancreas.R --resume` reuses saved fits only
after checking that their stored input matrices equal the reconstructed input.

The website pages load these saved results and do not rerun statistical fits
during a workflowr build. The fitness page is a concise derivative of the
MPCurver 0.3.0 public fitness vignette at repository commit `4a40eaf`.
`extract_fitness.py` copies three figures and extracts the 45-environment
assignment table and both reported ordering correlations from that article. The
SHA-256 source manifest is in `data/site_assets/fitness/source_manifest.csv`.
The complete vignette remains the canonical analysis.

| Source path in MPCurver | SHA-256 |
| --- | --- |
| `docs/articles/fitness.html` | `5fdaf67630e037ce0bccab25573ebf06ac7eab05dfa6fc6451756b03df25ea9e` |
| `docs/articles/fitness_files/figure-html/compare-single-orderings-1.png` | `a5f2fd5172fc6a98f42d713a91db0e16dad8f87823443dca5146ce7021747861` |
| `docs/articles/fitness_files/figure-html/partition-plot-1.png` | `f75cd8e1ded4ab597398073cba650a117a62aef6079f66a93f72a1763f584fd5` |
| `docs/articles/fitness_files/figure-html/known-sd-partition-trajectories-1.png` | `75f57048ad9b4fce774503ea26aed87850cbb9bbd3171222656e42b357faed13` |
