# Current-package-only smooth estimate-M report

The public one-anchor/nonmonotone estimate-M page previously combined current
MPCurver 0.3.4 automatic-M results with older MPCurver 0.3.2 comparison rows and
representative fits. This made an old imperfect ordering-B fit appear to describe
the current package.

The page now reads only the 90 validated MPCurver 0.3.4 automatic-M results in
`experiments/estimate_intrinsic_m_smooth_v032/auto_m_v034/`. New reporting code
generates current-only recovery, structural-recovery, runtime, and representative
fit figures. The old fits remain in the experiment archive for provenance but no
longer appear on the public page.

The current package recovers the final effective M in 88 of 90 datasets, with
two underestimates and no overestimates. Mean feature-partition ARI is 0.9945,
mean ordering recovery is 0.8390, and all 90 fits converged without warnings.
For `main_M3_S4_r001`, all three matched orderings, including ordering B, have
absolute Spearman correlation about 0.998.

The complete eight-page site build succeeded. `scripts/check_current_site.py`
passed all page, image, terminology, provenance, and saved-result assertions;
the checker now also rejects the old comparison labels on this page. The new
current-package figures were visually inspected, and `git diff --check` passed.
