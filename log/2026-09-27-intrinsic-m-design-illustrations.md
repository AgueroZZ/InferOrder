# Intrinsic-M simulation design illustrations

- Agent: Codex
- Update date: 2026-09-27
- Request: illustrate simulated trajectories under every setting in the Simulation design section.

## Changes

- Added `experiments/estimate_intrinsic_m_v032/render_design.R` to draw the saved replicate-one inputs for all nine true-M/SNR conditions.
- Added three figures, one for each true M, with ordering rows and SNR columns. Each panel shows three features, their true signals, and all 300 noisy observations on common axes.
- Recorded the exact selected features and input hashes in `main_summary/design_illustration_features.csv`; selection uses the first three features per group without consulting fitted results.
- Inserted the figures and a concise explanation into `report_body.Rmd` under Simulation design. The caption explains that conditions use independent datasets.
- Updated `render_report.R` to regenerate these illustrations before rendering the standalone report; PNG and PDF exports are saved in `main_summary/`.

## Verification

- Render completed successfully; inspected all three PNGs.
- Confirmed all three figures are embedded in the HTML before Methods, covering nine datasets and 108 selected features.
- Simulation data, fits, and numerical results were not modified or rerun. The original execution checksums remain the provenance for the completed fitting run.
- Updated the local report; no website publication in this update.
