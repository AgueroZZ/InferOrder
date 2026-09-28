# Correct the reported intrinsic ordering count

- Agent: Codex
- Update date: 2026-09-27
- Request: an unused ordering must not count toward the estimated intrinsic M.

## Correction

The previous report used EB occupancy for the adaptive method but selected
model dimension for uniform forward. The user clarified that the target is
effective M. Both methods now use the number of posterior assignment columns
with mean probability greater than 1e-12, matching the existing EB occupancy
threshold. Selected model dimension remains a separate diagnostic.

Recomputed occupancy for all 180 saved method outcomes. Only
`main_M3_S16_r006_forward` changes: selected model dimension 4, effective M 3.
Its posterior masses are (20, 0, 20, 20). For all 180 results, the posterior
occupancy count also equals the number of nonempty MAP groups.

- Adaptive EB exact effective-M recovery: 87/90.
- Uniform forward exact effective-M recovery: 90/90.

## Artifacts and provenance

- Added `reporting_metrics.R`; updated summary/example rendering to use it.
- Updated Methods, recovery figures/tables, and the standalone report.
- Removed the misleading ARI paragraph that described the unused slot as an
  overestimate of the target effective M.
- Preserved original RDS files and the prior summary tables. Current runs.csv
  includes original_reported_M, selected_model_M, and effective_M separately.
- No fitting, dimension-search decisions, ARIs, or ordering-recovery values
  changed. This supersedes the reporting interpretation in the earlier
  `2026-09-27-intrinsic-m-ari-empty-ordering.md` note.
