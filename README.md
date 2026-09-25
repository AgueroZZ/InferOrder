# InferOrder

InferOrder is a [workflowr](https://github.com/workflowr/workflowr) site for
selected MPCurve simulation and analysis results. The current pages use
MPCurver 0.3.0 and cover the method, one- and two-ordering simulations,
mutant fitness, and pancreatic semi-NMF loadings.

The public pages are in `analysis/`; generated HTML is in `docs/`. Run the
targeted site build from the repository root after the result artifacts are
available:

```sh
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 \
  Rscript --vanilla scripts/build_current_site.R
python3 scripts/check_current_site.py
```

The site build reads saved results and does not fit models. Reproduction
commands, inputs, fitting settings, and provenance for the displayed results
are in [`experiments/mpcurve_v030_site/README.md`](experiments/mpcurve_v030_site/README.md).
Legacy workflowr sources are retained in `archive/legacy-workflowr/analysis/`
so they cannot be included by a broad build of the current site.
