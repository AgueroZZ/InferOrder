# Intrinsic-M main simulation execution plan

**Goal:** Run the approved 90 paired datasets using MPCurver 0.3.2 and two total CPU threads.

**Architecture:** A separate `experiments/estimate_intrinsic_m_v032/` directory reuses the approved seed registry, generator, fitting protocol, and evaluation scripts. Two single-thread workers process disjoint alternating rows. A detached supervisor records process status and generates local summaries after both workers finish.

**Tech Stack:** R, MPCurver 0.3.2, Python subprocess management, ggplot2, R Markdown.

**Authorization:** The user approved two threads and starting the study. Execute inline; no additional planning approval is required.

## Fixed design

- True M = 3, 4, 5; SNR = 1, 4, 16; ten repeats per condition.
- N = 300, D = 60, K = 50, maximum fitted M = 8.
- Same original main-study seeds and same input matrix for both methods.
- Similarity initialization: absolute Spearman, single linkage, within-group Isomap.
- Adaptive EB occupancy and fixed-uniform-prior forward selection retain their existing definitions.
- Normalized absolute ELBO change / (N * D), tolerance 1e-6, with package convergence safeguards.
- Keep the existing 1,500-sweep continuation blocks and 10,000-sweep cap.

## Execution steps

- [x] Copy the existing runner, scoring, and reporting scripts into the new experiment directory; update only the package pin, stopping rule, paths, and worker dispatch.
- [x] Build and install a frozen 0.3.2 source archive in the experiment-local library; record source commit and checksums.
- [x] Prepare all 90 datasets from the original seed registry. Validate dimensions, group balance, noise scales, and agreement with the nine previously saved pilot inputs.
- [x] Validate disjoint worker assignments: 45 datasets per worker, five replicates per condition, no overlap and complete coverage.
- [x] Start a detached supervisor with two single-thread R workers and saved process identifiers; verify actual fitting activity and resource settings.
- [ ] After both workers finish, validate saved candidate/result hashes, convergence, normalized stopping settings, and the 180 planned method outcomes, then generate summary figures and a standalone local report.
- [ ] Review scientific results and figures before publishing the Simulation page in a subsequent review step.

## Commands

Run from the InferOrder root, with all BLAS/OpenMP thread limits set to one:

```sh
Rscript --vanilla experiments/estimate_intrinsic_m_v032/prepare.R
Rscript --vanilla experiments/estimate_intrinsic_m_v032/validate_setup.R
Rscript --vanilla experiments/estimate_intrinsic_m_v032/run_study.R --phase=main --worker=1 --workers=2
Rscript --vanilla experiments/estimate_intrinsic_m_v032/run_study.R --phase=main --worker=2 --workers=2
```

The supervisor owns both worker commands and runs `validate_results.R`, `summarize.R main`, `render_examples.R`, and `render_report.R` after successful execution. Completed scientific failures remain reported outcomes; process or validation failures stop finalization and are recorded in `run_status.json`.
