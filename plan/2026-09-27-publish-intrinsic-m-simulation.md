# Publish the intrinsic-M simulation

**Goal:** Add the reviewed 90-dataset simulation to the InferOrder Simulation section and publish it.

**Architecture:** A workflowr page includes the experiment's shared report source and reads its saved summaries and figures. The existing curated site builder refreshes navigation on all public pages. Commit only this study and its publication changes.

**Tech Stack:** R Markdown, workflowr, saved MPCurver 0.3.2 results, GitHub Pages.

## Steps

- [x] Add `analysis/estimate_intrinsic_m.Rmd`, a Simulation menu item, and an index link.
- [x] Summarize the reviewed conclusion: both approaches perform well for balanced groups of independent, mostly mildly curved monotone trajectories; effective-M recovery is 87/90 for adaptive EB and 90/90 for uniform greedy forward.
- [x] Preserve the common posterior-occupancy definition and all nine setting illustrations. Keep fitting provenance distinct from the corrected reporting metric.
- [x] Make the experiment's preparation script reproducible without the unpublished historical experiment; include a checked pilot hash reference.
- [x] Extend `scripts/build_current_site.R` and `scripts/check_current_site.py` for seven pages and verify saved-result claims, navigation, and assets.
- [x] Rebuild the standalone report and workflowr site from saved results, then inspect the rendered page.
- [x] Commit the simulation bundle, relevant plans/logs, page source, and generated public pages; preserve unrelated edits and exclude installed libraries and full fitting states.
- [x] Push to origin/master, verify GitHub Pages deployment and the live page, and report the public URL.
