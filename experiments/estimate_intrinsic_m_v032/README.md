# Intrinsic-M main simulation with MPCurver 0.3.2

This study runs the approved 90 paired datasets: true M = 3, 4, 5;
variance SNR = 1, 4, 16; ten repeats per condition; N = 300 and D = 60.
Both methods use K = 50, similarity/Isomap initialization, and M <= 8.
Adaptive EB and fixed-uniform-prior forward selection use the same input
matrix for each dataset. The stopping rule is absolute ELBO change / (N * D)
below 1e-6, with MPCurver's annealing and numerical nondecrease safeguards.

## Reported dimension

Following the user's September 27 clarification, both methods report effective
M as the number of columns of the posterior assignment matrix whose mean
exceeds 1e-12. `reporting_metrics.R` implements this shared definition.
`main_summary/runs.csv` retains the selected model dimension separately in
`selected_model_M`, alongside `effective_M` and `original_reported_M`.
The original fit/result RDS files retain their execution-time fields; use the
reporting helper when interpreting their estimated counts. Pre-correction
tables are preserved under `reporting_audit/before_effective_m_correction/`.
This reporting correction changes one uniform-forward outcome from 4 model
slots to 3 occupied orderings. It does not alter the forward search or any fit.

The original 0.3.1 experiment is preserved separately. `manifest.csv` retains
its complete seed registry; `active_manifest.csv` selects the first ten
main repeats per condition. The legacy `repetitions = 50` field describes
the seed registry; `planned_repetitions = 10` specifies this execution.
All nine pilot input matrices and latent positions were regenerated exactly
as a data-provenance check. Those pilots are excluded from the main study
and are not refitted here. The package's completed 0.3.2 regression study
provides the fitting validation preceding this run.

## Reproduction

Run from the InferOrder root with R_PROFILE_USER and R_ENVIRON_USER set to
/dev/null and all BLAS/OpenMP thread limits set to one. Install the frozen
`source/MPCurver_0.3.2.tar.gz` in the experiment-local `library/` directory.
The source Git commit is recorded in `source/source_commit.txt`.

```sh
Rscript --vanilla experiments/estimate_intrinsic_m_v032/prepare.R
Rscript --vanilla experiments/estimate_intrinsic_m_v032/validate_setup.R
python3 experiments/estimate_intrinsic_m_v032/supervise.py
```

The supervisor runs two single-thread workers. Alternating manifest rows
give each worker 45 datasets and five replicates per condition. A file lock
prevents simultaneous supervisors. Logs append on restart; completed
candidates are reused and unfinished candidates resume from checkpoints.
Scientific failures remain reported outcomes rather than being silently
retried. The existing cap is 10,000 sweeps per candidate.

## Monitoring and outputs

- `run_status.json`: supervisor/worker process IDs, phase, and completed count.
- `worker_1.log`, `worker_2.log`, and `logs/`: dispatch and candidate progress.
- `dispatch_manifest.csv`: the fixed worker assignments.
- `execution_checksums.json`: code, configuration, and source archive hashes.
- `data/` and `input_hashes.csv`: all paired input matrices and fingerprints.
- `candidates/`, `results/`, `checkpoints/`, `full_fits/`: resumable fit artifacts.
- `main_validation.json`: complete-study validation after both workers finish.
- `main_summary/`: recovery distributions, accuracy, structure, and runtime figures.
- `analysis_report.html`: standalone local report, generated after validation.
- `finalize.log`: validation, summarization, and rendering output.

The computer must remain awake for continuous execution. Closing this chat
does not terminate the detached supervisor. The local report and scientific
results should be reviewed before publishing the workflowr Simulation page.
