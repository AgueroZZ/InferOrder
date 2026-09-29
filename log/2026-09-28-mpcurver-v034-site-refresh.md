# MPCurver 0.3.4 public-analysis refresh

- Agent: Codex
- Date: 2026-09-28 UTC
- Request: regenerate the current InferOrder simulation and observed-data
  pages with MPCurver 0.3.4, then publish the updated site.

## Scope and provenance

The introductory one- and two-ordering simulations, fitness summary, and
pancreas analysis were regenerated with MPCurver 0.3.4 from package commit
`15f2b0bbe5dfa61cd46da5160b2bc251e75a0475`. The frozen source archive is
`experiments/mpcurve_v034_site/source/MPCurver_0.3.4.tar.gz`, with SHA-256
`58d0c99c6170994c82eedba190fb1c28a63eb47530c575705323c3cf3992b7c6`.
Scripts, fit objects, inputs, settings, summaries, figures, session
information, and package provenance are retained under
`experiments/mpcurve_v034_site/`.

The formal intrinsic-M comparison pages remain frozen at their documented
MPCurver 0.3.2 and 0.3.3 versions. Their sources, saved results, and rendered
HTML were not changed by this refresh.

## Results

- The one-ordering simulation converged in 55 recorded iterations and had
  absolute Spearman correlation 0.9997 with the generating ordering.
- The two-ordering simulation converged with 369 recorded temperature-one
  objective values, recovered all 20 feature assignments, and had absolute
  ordering correlations 0.9975 and 0.9989.
- Both pancreas fits converged. The two-ordering fit retained the seven-versus-
  one feature split, with factor 20 assigned to ordering B, and used 162
  recorded temperature-one objective values. Correlations of the one-ordering
  positions with orderings A and B were 0.9679 and 0.0576.
- The refreshed fitness summary assigns 27 environments to ordering A and 18
  to B. Its two inferred orderings have correlation 0.857, and the two
  single-ordering smoothness fits have correlation 0.952.

The package fitness article was rebuilt in a valid UTF-8 locale before its
results were extracted. This preserves the `μg/ml` and `μM` units in the
environment labels. The article and copied-figure hashes are recorded in
`data/site_assets/fitness/source_manifest.csv`; the documentation-only source
commit is `eeaf65bcf4c8b021d590c88ad411506359202f1f`.

## Verification

All eight workflowr pages were rebuilt from their R Markdown sources. The
site checker validated page structure, titles, navigation, local links,
embedded-image counts, version labels, package provenance, archive checksum,
convergence flags, assignments, and reported values. `git diff --check`
passed. The principal simulation and pancreas figures were visually inspected
for readable labels, unclipped panels, and agreement with the saved summaries.

InferOrder commit `036298f4c3b78ca73a1445298706e192a10e0ce5` was pushed
to `origin/master`. GitHub Pages deployment
[36501337473](https://github.com/AgueroZZ/InferOrder/actions/runs/36501337473)
completed successfully for that exact commit. The live homepage, method page,
two simulation pages, fitness page, and pancreas page are byte-identical to
their committed HTML files. Live content checks found MPCurver 0.3.4 and the
reported recovery, assignment, and correlation values on the corresponding
pages.
