# One-step ELBO selection within a feature group

The user's intended selector fits each candidate ordering to its own feature
group for one single-ordering MPCurve update, then compares that group's ELBO.
This differs from the first joint structural update studied in `../early_selection/`.
The earlier joint-update counterexample does not test this within-group selector.

On the original B observations this local selector chooses Isomap k=10, the
candidate with the best previously converged joint-model objective and B recovery.

| Candidate | B-only ELBO after one update | Previously converged joint B recovery |
| --- | ---: | ---: |
| k=5 | -3884.95 | 0.996793 |
| k=10 | -3875.18 | 0.996823 |
| k=15 | -4351.64 | 0.700139 |
| k=20 | -4260.16 | 0.688720 |
| k=30 | -4549.98 | 0.676222 |
| PCA | -4500.57 | 0.145150 |

## Exact scoring protocol

Use only the 300-by-12 B matrix. For each Isomap k or PC1, obtain raw positions,
quantile-bin them into K=50 grid points, and use the existing frozen 0.3.4
`.cavi_build_from_ordering()` helper with `max_iter=1`. All candidates use RW2,
ridge zero, initial smoothing precision 1, and the same initialization rule for
noise and position probabilities. As in the package initializer, initial noise
estimates are computed from each candidate's bins, rather than held numerically
identical across candidates. This is the candidate's fully evaluated variational
objective including its fitted parameters, not a raw embedding score.

Each result is checked to have exactly one CAVI iteration and an ELBO trace
of length two: the initialized state followed by one update. There is no
feature-assignment uncertainty across the other four orderings and no q(Z)
annealing term. The first-step scores are ordinary single-ordering ELBO values.
The same data, group size, model, grid, and priors make comparisons across
candidate initializations within each noise setting meaningful. They should
not be compared numerically with the five-ordering joint-model ELBOs.

All 17 connected candidates in the three existing B-noise settings were scored:

| B noise multiplier | Selected candidate | Prior joint final B recovery | Loss in prior joint final ELBO versus its best candidate |
| --- | --- | ---: | ---: |
| 0 | k=10 (connected Isomap candidates tied) | 0.999867 | 0 |
| 0.5 | k=15 | 0.999089 | 0.001443 |
| 1 | k=10 | 0.996823 | 0 |

The disconnected noiseless k=5 candidate remains excluded. All 17 one-step fits
completed without warnings. Local initialization plus one update took about
0.010--0.054 seconds per candidate in this run, excluding embedding computation,
R startup, and result handling; this is not an end-to-end speedup benchmark.

## Scope of the conclusion

This is positive evidence for the user's within-group screening proposal on
this selected example. A group's initial candidates are judged directly by how
well their own trajectories explain that group's features, before joint feature
reassignment or initialization of other groups enters the optimization.

The downstream outcomes in the tables are the previously saved fits, which
used two subset CAVI updates before joint fitting. This experiment validates
candidate selection against those outcomes; it does not separately validate a
modified end-to-end pipeline that passes only the one-update fit to the joint
model. One could retain the selected candidate and complete the existing warmup.

No claim of reliability across independent datasets or incorrectly partitioned
feature groups follows from these three noise transformations of one dataset.
The package initializer remains unchanged. A broader test would apply this
within-group score independently to each automatically initialized group and
compare selected starts and resulting complete fits on independent datasets.

## Reproduce

From the InferOrder root:

```sh
bash experiments/estimate_intrinsic_m_smooth_v032/run_r.sh \
  experiments/estimate_intrinsic_m_smooth_v032/exploratory_m5_ordering_b/within_group_selection/score_one_step.R
```

`one_step_scores.csv` stores local scores and matched downstream outcomes;
`selection.csv` stores the within-noise winners. `one_step_results.rds` retains
initial/final positions, score traces, fitted parameters, feature identities,
input hashes, seed, warnings, and the frozen package commit. All package source,
original fits, and website files remain unchanged.
