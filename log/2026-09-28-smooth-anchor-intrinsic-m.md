# Smooth trajectories with one monotone anchor per ordering

## Design and provenance

Prepared 90 paired datasets with N = 300, D = 60, true M = 3, 4, 5, variance
SNR = 1, 4, 16, and ten repeats per condition. Original latent positions,
feature assignments/order, noise realizations, and one monotone anchor per
ordering are retained. Other features use smooth random Fourier mixtures.
All 5,040 replacement curves have both analytic derivative signs and sampled
turning points. Every input independently regenerated identically.

Design hash: `bc315aa0089d96ff1dae6f572118a9e8704620450ca034d4205ff9f38e163777`.
The same frozen MPCurver 0.3.2 archive and original fitting settings were used.

## Execution and validation

The experiment-local library uses the existing fixed R 4.3.3 environment;
princurve 2.1.6 is installed from its retained source archive. The default
R configuration had mismatched startup/compiler files, so `run_r.sh` selects
the verified environment explicitly. Shared caslake submission lacked SU
allocation; the study used the available zemmour-hm private partition.

Slurm array 59628523 ran 90 one-thread tasks, with at most 12 simultaneous
workers. All tasks completed with exit status zero. Aggregate task time was
21,297 seconds; runtime details are recorded in
`execution_validation.json` and `slurm/accounting.txt`. Frozen fitting code
and the software source archive remained unchanged during execution.

Complete validation passed for all 180 method outcomes and 570 candidates:
all converged, zero fitting warnings/errors, no objective-guard violations,
and consistent forward-selection decisions. Paired baseline data and common
occupied-ordering evaluation were independently verified.

## Findings

- Adaptive EB: exact effective M in 60/90 datasets versus baseline 87/90;
  29 underestimates and one overestimate. At SNR = 1, recovery is 0/10 for
  true M = 4 and M = 5; at SNR = 4, true M = 5 recovery is 3/10.
- Uniform + forward: 88/90 versus baseline 90/90. Both errors estimate four
  occupied orderings when true M = 5 and SNR = 1.
- Mean feature-partition ARI: adaptive 0.730 versus 0.996 baseline;
  forward 0.995 versus 1.000.
- Mean matched ordering recovery: adaptive 0.685 versus 0.985 baseline;
  forward 0.843 versus 0.985. Forward ordering recovery averages 0.608 at
  SNR = 1 and 0.999 at SNR = 16.
- Median within-ordering absolute Spearman similarity falls from 0.800 to
  0.398; between-ordering similarity remains approximately 0.039. This is a
  descriptive diagnostic of the similarities used for initialization.

## Website preparation

Added `analysis/estimate_intrinsic_m_smooth.Rmd` and its shared report source,
plus navigation and curated-builder/audit entries. All eight pages passed
links, image, terminology, and saved-result checks. Inspected actual plots,
wrapped clipped example titles, and made trajectory-role colors consistent.
Desktop (1280 by 1000) and mobile (390 by 844) browser inspection passed: all
10 images and captions loaded, equations rendered, reported counts matched saved
results, and the page had no JavaScript errors or mobile horizontal overflow.
Wide tables scroll within the page; the curve formula and coefficient definition
are separated for readability. Evidence is in `visual_validation.json`.

The result page is prepared for publication at
https://aguerozz.github.io/InferOrder/estimate_intrinsic_m_smooth.html.
Live deployment and the publication commit are recorded in the workspace
`PROGRESS.md` after push.
