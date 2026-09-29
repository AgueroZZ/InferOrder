# Ordering-B diagnostic for `main_M3_S4_r001`

This exploratory diagnostic investigates why the original adaptive-EB fit has
an absolute Spearman correlation of about 0.76 for true ordering B, although
it exactly recovers all three feature groups. It does not change the completed
simulation fits or the public workflowr page.

The script reconstructs the MPCurver 0.3.2 absolute-Spearman, single-linkage,
Isomap initialization and compares it with the saved adaptive, uniform-forward,
and automatic-M fits. It also tests how the ordering-B Isomap initialization
changes when features excluded from its surviving 16-feature initialization
cluster are added back one at a time. This subset experiment uses truth only as
a diagnostic and is not a proposed estimator.

Run from the InferOrder repository root:

```sh
bash experiments/estimate_intrinsic_m_smooth_v032/run_r.sh \
  experiments/estimate_intrinsic_m_smooth_v032/exploratory_ordering_b/diagnose_ordering_b.R
```

The script writes the diagnostic figure and CSV tables to this directory.

## Findings

The ordering is identifiable from the complete B group: the original
uniform-forward fit and the later automatic-M adaptive fit both recover B with
absolute Spearman correlation 0.998. The original adaptive fit also assigns all
features to the correct true group (ARI 1), so feature misclassification is not
the cause.

Instead, the M = 8 absolute-Spearman tree splits the 20 B features into initial
clusters of sizes 3, 16, and 1. The 16-feature component becomes the surviving
ordering, but its initial Isomap position has absolute Spearman correlation
0.757 and its final adaptive-EB position remains at 0.756. The smaller
three-feature component contains the monotone anchor and has a better initial
ordering (0.932), but that component collapses during adaptive fitting.

This is not explained by larger noise in B: mean realized noise variance is
0.253, 0.246, and 0.251 for A, B, and C. It is a geometric initialization
failure. Adding the excluded nonmonotone feature V28 to the 16-feature core
raises the Isomap ordering correlation from 0.757 to 0.998. The correlation
between its 15-nearest-neighbor graph distances and true latent distances rises
from 0.769 to 0.986, nearly matching the complete B group (0.988). V28 is
selected with truth only to localize the failure; it is not a fitting rule.
