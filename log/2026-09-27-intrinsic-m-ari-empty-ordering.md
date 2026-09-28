# ARI and an unused ordering in uniform forward selection

- Agent: Codex
- Update date: 2026-09-27
- Request: explain why uniform forward has ARI 1 despite one overestimated M.

## Verified result

For `main_M3_S16_r006_forward`, the selected model has M = 4 but its posterior
assignment masses are (20, 0, 20, 20). The MAP feature counts are also
(20, 0, 20, 20). True groups A, B, C map exactly to slots 3, 1, 4,
respectively. Slot 2 has zero assignment probability for every feature.
The fixed prior remains uniform at (0.25, 0.25, 0.25, 0.25).

The stored estimated-M metric for uniform forward is the selected model
dimension, whereas ARI compares the true labels with the per-feature maximum
posterior assignment. Its three nonempty groups exactly recover the truth;
the unused fourth ordering contributes no cluster disagreement. Recomputed
all 180 ARIs from saved posterior assignment matrices and true labels;
every value matches the stored result within 1e-12.

The forward history accepted M = 4 over M = 3 with an ELBO gain of 1.960623,
then rejected M = 5 with a change of -15.455938. These are observed fitted
objectives; no additional optimization was performed in this audit.

## Documentation update

Clarified the MAP-label definition of ARI and the specific empty-ordering
example in the report's structure-recovery section. Re-rendered the local
report. No metrics, selection rules, data, or fits were changed.
