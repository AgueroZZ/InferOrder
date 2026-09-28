# Intrinsic-M main simulation launch

- Agent: Codex
- Update date: 2026-09-26
- Authorization: the user approved two total CPU threads and starting the previously agreed study.
- Plan: `plan/2026-09-26-intrinsic-m-v032-main.md`.

## Configuration and provenance

- Separate experiment: `experiments/estimate_intrinsic_m_v032/`; prior 0.3.1 artifacts are preserved.
- 90 main datasets, 180 method outcomes: true M = 3, 4, 5; SNR = 1, 4, 16; ten repeats per cell.
- N = 300, D = 60, K = 50, M_max = 8; same generator, seeds, similarity/Isomap initialization, EB occupancy definition, and uniform-forward comparison protocol.
- MPCurver 0.3.2 built from clean tracked source at `dd1a903c9c316717419b18bb2b5f325c3a4c9185` and installed into an isolated library.
- Explicit normalized ELBO stopping at 1e-6; original continuation blocks and 10,000-sweep cap retained.
- Two disjoint workers, each with 45 datasets and five repeats per condition; BLAS/OpenMP and package thread counts are one per worker.

## Validation before launch

- Original seed registry verified column by column after CSV type normalization.
- All 90 main inputs saved; all nine separate pilot matrices and latent positions reproduce the old artifacts exactly.
- Existing data-generation and label/reversal-invariant scoring checks passed.
- All R scripts parsed; pinned package path, version, and normalized fit options verified.
- Supervisor syntax validated. A file lock prevents duplicate supervisors; candidates and fitting blocks retain atomic checkpoint writes.

## Launch and follow-up

- Started at 2026-09-26 23:33:27 America/New_York (2026-09-27 03:33:27 UTC).
- Detached supervisor PID: 75994; worker dispatcher PIDs: 75995 and 75996.
- Verified two active fitting processes at approximately 98% CPU each, running `main_M3_S1_r001` and `main_M3_S1_r002`; dispatchers and supervisor were idle.
- `run_status.json`, worker logs, and dataset logs provide progress. No completed main outcomes were claimed at launch.
- After both workers finish, the supervisor validates results and generates local tables, figures, and a standalone report. Process or validation failures are recorded and prevent successful finalization.
- Scientific review and workflowr publication remain subsequent steps.
