# Publish the intrinsic-M simulation

- Agent: Codex
- Update date: 2026-09-27
- Authorization: publish the reviewed simulation results and push the website.
- Plan: `plan/2026-09-27-publish-intrinsic-m-simulation.md`.

## Public content

- Added `analysis/estimate_intrinsic_m.Rmd`, linked from the Simulation menu and homepage.
- The page shares the standalone report source, including simulated trajectories for all nine settings and the corrected effective-M results.
- Summary: both methods perform well for these balanced, independent, mostly mildly curved monotone settings; uniform greedy forward recovers effective M in 90/90 datasets and adaptive EB in 87/90.
- Removed the homepage's site-wide 0.3.0 label; each study documents its own fitting version, preserving the older studies' provenance.

## Reproducibility and validation

- The committed study bundle includes code, fixed seeds, synthetic inputs, compact candidate/result artifacts, summary figures, reports, and the pinned 0.3.2 source archive. Installed libraries and large full fitting states remain excluded.
- Replaced a preparation-time dependency on the unpublished historical experiment with verified reference hashes for all nine pilot inputs and latent positions.
- All fits and metrics are reused; no simulations were rerun for publication.
- Extended the curated builder and site audit from six to seven pages.
- Site audit passed for all seven pages, navigation, images, terminology, and numerical claims. The new page includes ten embedded figures.
- Local Chrome inspection confirmed the page layout, rendered equations, figure display, effective-M explanation, and 87/90 versus 90/90 results.
- Fixed the child-report figure-path configuration so workflowr no longer inserts visible warnings. A network-restricted favicon fetch warning remains confined to the build log.
- `git diff --check` passed before staging. Unrelated existing edits are excluded from this publication.

## Deployment

Published commit `6d1fb4e719146e8854afe14aafcadf1e7c1ba7d0` to origin/master.
GitHub Pages run `36366495048` completed successfully for that commit.
A live Chrome visit verified the new simulation page, its figures and
conclusion, including 90/90 versus 87/90 effective-M recovery:
https://aguerozz.github.io/InferOrder/estimate_intrinsic_m.html

The simulation is accessible from the homepage and the Simulation menu.
The subsequent documentation-only commit records these completed checks.
