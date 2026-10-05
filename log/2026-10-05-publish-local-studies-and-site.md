# Publish accumulated experiments and research-site edits

With user authorization, prepare the accumulated local experiment code,
saved results, diagnostics, and research-site edits for publication on master.
The comparison page uses concise design specifications, shows PC1 against true
positions in Figure 1, and presents paired noise conditions in seven tables.
Shared navigation follows the five studies selected by the homepage.

Retain the exploratory ordering-B analyses, local-neighborhood experiments,
Isomap screening, twenty-replicate two- and four-feature comparisons, and
replicate-6 geometry/noise diagnostics in their existing experiment locations.
Their inputs, fitting versions, seeds, frozen source archives, plotted data,
and historical verification remain unchanged. Historical experiments using
MPCurver 0.4.0.9000 keep that actual fitting version. Generated internal pages
remain outside the main navigation. Python bytecode caches are now ignored;
R session-information files are retained as scientific provenance.

## Current publication checks

- Rebuilt all ten saved-result pages with `scripts/build_current_site.R` using
  the existing workflowr configuration; refreshed both retained exploratory
  pages' navigation. The final build completes without resource-fetch warnings.
- `scripts/check_current_site.py` passes its eleven-page result/link/image
  audit and checks navigation across all twelve public HTML pages.
- Independently checked all seven comparison tables against committed CSV
  summary inputs, including method labels, paired rows, precision, and headers.
- Browser inspection of the homepage and five main studies finds no image
  loading failures, math-rendering errors, or desktop page-width overflow.
  The edited comparison page also passes a 420-pixel viewport check; its
  Figure 1, caption, and grouped tables were visually inspected.
- Page sources, builders, and public HTML pass whitespace checks; frozen
  experiment snapshots and session-information outputs retain their original
  whitespace and CSV line endings to preserve provenance. A credential-pattern
  scan finds no matches. Inspected
  added paths for local-only artifacts; caches, fitting libraries, large full
  fits, and temporary logs follow their existing ignore rules.

Temporary checks and screenshots are under `/tmp/inferorder-v041-*` and
`/tmp/check_grouped_m1.py`. No experiment is refitted during publication.
