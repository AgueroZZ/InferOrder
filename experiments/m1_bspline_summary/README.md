# Concise three-figure simulation summary

The canonical page is `analysis/simulation_summary.Rmd`. It summarizes the
existing MPCurve/principal-curve/GPLVM/Bayesian-GPLVM experiments with N=200,
P=12/50, two noise levels, 30 replicates, and common PCA or Isomap starts.
No new model fits or data-generating scenarios are introduced.

- Figure 1: fixed P=12 replicate 1, features 1-3; true spline functions and
  generated observations against true t. The same curves and noise draws appear
  at SD .25 on the left and SD 1 on the right. These illustrate easier/harder
  noise conditions, not different spline-complexity regimes. Values precede
  observation-column centering. Selection does not use fitting outcomes.
- Figure 2: 2-by-2 PCA boxplots, with P=12/50 as rows and low/high noise as
  columns. PCA is the shared centered, unscaled-feature PC1 from the original
  comparison, explicitly supplied to the GP implementations.
- Figure 3: matched 2-by-2 Isomap plots, k=15. Same input data and fixed settings;
  the principal-curve geometric-start adaptation is disclosed in the caption.

All 1,200 model/baseline endpoints are retained, with nonconvergence and the
numerical smoother warning marked. Plot data and source SHA256 hashes are
saved beside the generator. `figures/` contains PNG and PDF versions. The
longer `analysis/simulation_m1_comparison.Rmd` retains settings and full results.

From the InferOrder root:

```bash
bash experiments/estimate_intrinsic_m_smooth_v032/run_r.sh experiments/m1_bspline_summary/build_figures.R
bash experiments/estimate_intrinsic_m_smooth_v032/run_r.sh experiments/m1_bspline_summary/build_site.R
python3 scripts/check_current_site.py
```
