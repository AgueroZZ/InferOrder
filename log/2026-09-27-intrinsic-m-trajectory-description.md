# Clarify the simulated trajectory shapes

- Agent: Codex
- Update date: 2026-09-27
- Request: reconcile the trajectory-family description with the nearly linear shapes in the illustrations.

## Findings

Inspected `simulate_intrinsic_trajectories()` and `.simulate_dual_signal_block()`
from the experiment-local MPCurver 0.3.2 installation. The monotone generator
forms each curve from a positive weighted sum of five bases: two power bases,
a logistic basis, an exponential basis, and a Hill-type basis. Weights are
independent exponential draws normalized to sum to one; direction is randomly
reversed. These are components within each curve, not separately sampled
trajectory categories.

Across all 90 saved main datasets (5,400 noiseless feature trajectories),
linear-fit R-squared, computed as squared Pearson correlation between each
signal and its assigned true latent position, has median 0.9879811.
The fraction with R-squared >= 0.95 is 0.9057407; the minimum is 0.8109475.
The near-linear appearance is therefore characteristic of this benchmark,
not solely a consequence of showing three features per group.

## Update

Replaced the list of function families in `report_body.Rmd` with a description
of random monotone basis combinations, random direction, and predominantly
mild curvature. Re-rendered the local standalone report. Simulation inputs,
fits, and recovery results are unchanged; no new fits were run.
