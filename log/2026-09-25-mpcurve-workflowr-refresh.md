# MPCurve workflowr refresh — 2026-09-25

## Scope

Replaced the public Smooth-EM presentation with six MPCurve pages: home,
method, one-ordering simulation, two-ordering simulation, fitness, and
pancreatic cell loadings. All current scientific fits use MPCurver 0.3.0 and
CAVI. The homepage and navigation now group the pages under Method,
Simulation, and Analysis.

The legacy analysis sources were moved to `archive/legacy-workflowr/analysis/`
so a broad build cannot republish them. Twelve tracked legacy HTML pages were
removed from `docs/`; two untracked legacy HTML pages were preserved under
`archive/legacy-workflowr/untracked-html/`. The public `docs/` directory now
contains only the six current top-level HTML pages. Unrelated code, gplvm,
data, and output changes in the working tree were left untouched.

## Results and provenance

- The one-ordering spiral simulation uses 1,000 samples, two observed
  features, noise standard deviation 0.12, and seed 1. Its fitted ordering
  has absolute Spearman correlation 0.9998 with generating position. This is
  one illustrative realization.
- The two-ordering simulation uses 500 samples with two ten-feature groups,
  noise standard deviation 0.05, and seed 1. The Fiedler-initialized fixed
  M = 2 fit has label-invariant feature partition accuracy 1.000 and absolute
  ordering correlations 0.9975 and 0.9989. A checked PCA initialization
  recovered the feature split but gave poor recovery for the second sample
  ordering; the displayed Fiedler result is specified explicitly.
- The pancreatic analysis reproduces the existing semi-NMF loading selection
  from `data/loading_order/pancreas_factors.RData` and `pancreas.RData`:
  at most 500 cells per cell type with seed 1, columns 3, 8, 9, 12, 17, 18,
  20, and 21, then the inDrop3 subset. The same finite 865-by-8 matrix is
  used for fixed M = 1 and M = 2 fits. The M = 2 fit assigns seven factors
  to ordering A and factor 20 to ordering B. Absolute Spearman correlations
  between M = 1 positions and M = 2 positions A and B are 0.967 and 0.177.
  These fits do not select the number of orderings.
- The fitness page summarizes the MPCurver 0.3.0 vignette at repository
  commit `4a40eaf`; it does not refit the data. The copied analysis has 421
  mutants and 45 environments, with 27 environments assigned to ordering A
  and 18 to B by largest assignment probability. The two fitted mutant
  position vectors have correlation 0.866. The vignette's Student-t4 to
  Gaussian-score calibration of measurement errors is stated as an analysis
  assumption. The source HTML and three copied figures have SHA-256 hashes
  in `data/site_assets/fitness/source_manifest.csv`.

The exact generators, seeds, fitting calls, serialized fits, numerical
summaries, figures, input hashes, and R session information are under
`experiments/mpcurve_v030_site/`. The pancreatic `--resume` path was
corrected to compare the saved fit's data matrix with the reconstructed input.

## Build and validation

- `Rscript --vanilla experiments/mpcurve_v030_site/run_simulations.R`
  completed for both simulations.
- `Rscript --vanilla experiments/mpcurve_v030_site/run_pancreas.R` completed
  for both fixed-M fits; `--resume` also completed after checking input
  identity without refitting.
- `python3 experiments/mpcurve_v030_site/extract_fitness.py --mpcurver-root ../MPCurver`
  completed, and every manifest hash matched the current source files.
- `Rscript --vanilla scripts/build_current_site.R` rendered all six pages
  using `workflowr::wflow_html` through `rmarkdown::render` and read only
  saved result artifacts. `python3 scripts/check_current_site.py` passed:
  titles, navigation, headings, local links, embedded images, terminology,
  exact public page set, and key saved-result claims.
- The pancreatic result figures and copied fitness figures were inspected
  as local images. `git diff --check` passed.

Two direct `workflowr::wflow_build()` attempts failed before rendering when
its subprocess log disappeared from the temporary directory. The targeted
render script is the reproducible site build path in this environment.
Pandoc still emits workflowr-template warnings about implicit div closure
and cannot fetch a workflowr badge image from GitHub while network access is
restricted; all six pages render and pass static checks. Browser visual QA
of local files was blocked by the browser security policy, so only local
image inspection and static HTML checks were performed.

At the initial review stage, the updated site was local and unpublished.

## Reader-facing publication pass

Before publication, the six pages were edited for direct interpretation of
their figures and results. The fitness summary now reports the 0.952 Pearson
correlation between its two single-ordering fits, extracted from the same
MPCurver vignette as the figures and assignment table. Repetitive caveats,
build details, and code for figure rendering were removed from the visible
pages; the underlying scripts and saved results remain linked.

The workflowr report and session-information panel were suppressed in the
public HTML through `_workflowr.yml`. A standard per-page figure path and a
clean knit environment removed visible workflowr warnings. The six pages
were rebuilt from saved results and passed `scripts/check_current_site.py`
and `git diff --check` before staging. A workflowr favicon fetch warning
remains in the local build log because the sandbox cannot reach GitHub, but
it is not displayed in the page body.

## Publication

The site refresh was committed as `40fa47e` and pushed to `origin/master` on
2026-09-25. GitHub Pages serves the updated site at
<https://aguerozz.github.io/InferOrder/>.

Public verification after deployment found:

- The canonical site root and all six current HTML pages return HTTP 200.
  Downloaded HTML for each current page matches its committed local file
  byte-for-byte by SHA-256.
- Fourteen checked historical HTML URLs, including the twelve removed
  tracked pages and the old `about.html` and `theory.html` paths, return
  HTTP 404 without query parameters.
- All nine distinct external links in the current pages to MPCurver methods,
  the InferOrder repository, reproduction scripts, and saved results return
  HTTP 200.
- A few immediate requests briefly returned old content from a CDN cache;
  subsequent requests to the canonical URLs returned the current pages and
  404 responses for the retired pages.

The unrelated modified and untracked files in the shared checkout were not
included in the publication commit. The locally retained untracked legacy
files and two exploratory pancreatic candidate fits also remain unpublished.
