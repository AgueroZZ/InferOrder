# Graph-based bounds for the Isomap neighborhood search

Read Samko, Marshall, and Rosin (2006), *Selection of the optimal parameter
value for the Isomap algorithm*, Pattern Recognition Letters 27:968–979,
[DOI](https://doi.org/10.1016/j.patrec.2005.11.017), particularly Section 3 and
Table 1 on page 970. The [author-uploaded full text](https://www.researchgate.net/profile/Paul-Rosin/publication/223824802_Selection_of_the_optimal_parameter_value_for_the_Isomap_algorithm/links/59f784c1aca272607e2d8632/Selection-of-the-optimal-parameter-value-for-the-Isomap-algorithm.pdf)
was accessible through the web reader; PDF screenshot retrieval failed.

The paper defines the lower bound by graph connectivity and gives a heuristic
upper-bound condition based on average degree: `2 * edge_count / n <= k + 2`.
Within its interval it finds reconstruction-cost minima, then selects among
them by input/output distance correlation. The authors report problems when
sparse sampling or noise causes shortcuts already at the connectivity threshold.

## Implications for our proposed ELBO selector

The following are our deductions and possible adaptations, not claims that the
paper establishes them for MPCurve. These were recorded at the reading stage,
before the graph evaluation reported below.

- Connectivity supplies a well-defined lower bound for fitting all samples
  with one connected Isomap graph. This does not guarantee correct geometry.
- Interpreting the upper bound literally as the largest qualifying k across
  all integers 1 through n-1 makes it n-1 for every dataset: at k=n-1 the
  simple graph is complete, its average degree is n-1, and the condition holds.
  The inequality is not a monotone exclusion rule. A small-k stopping convention,
  such as the first violation, would be an additional definition rather than
  an explicitly stated instruction from the paper.
- For the undirected union of directed k-nearest-neighbor lists, let A denote
  the number of directed neighbor entries whose reverse entry is absent.
  There are n*k directed entries and (n*k+A)/2 distinct undirected edges.
  Hence average degree is k+A/n, and the proposed condition is equivalent to
  A/n <= 2. It measures neighbor nonreciprocity rather than directly measuring
  the length or topological validity of shortcut edges.
- Count each unordered edge once. The current package constructs an undirected
  igraph from directed neighbor entries, retaining parallel edges for reciprocal
  pairs. Counting those entries directly would give average degree 2*k and
  incorrectly reject every k>2. A degree diagnostic must use a separate simple
  graph or deduplicated unordered pairs; existing shortest-path fits need no
  change merely to calculate this diagnostic.
- A plausible adaptation uses connectivity to define the lower bound, validates
  a local degree-based upper-bound convention, and lets the one-sweep MPCurve
  ELBO select within the resulting range. Before treating the upper bound as a
  hard cutoff, check whether it retains successful starts such as k=10 in the
  current case. An interval cannot guarantee the best downstream ELBO.

This reading leaves the upper-bound convention unresolved. The earlier
5:30 experiment, its selection, and its timing measurements are unchanged.

## Fixed-case graph evaluation: 2026-10-02

Evaluated all k=1 through 299 on the exact original 300-by-12 B-group matrix.
The graph is the simple undirected union of Euclidean k-nearest-neighbor lists;
reciprocal neighbor pairs count as one edge. No observations were scaled or
filtered. The neighbor query uses RANN's exact kd-tree search with eps=0.
Bound calculations use observations only; no true positions or ELBO enter them.

Table 1. Connectivity and the paper's degree condition at selected neighborhood
sizes. Average degree is twice the number of unique undirected edges divided
by 300; the condition passes when average degree <= k+2. Values are rounded
only for display; the condition is checked using integer edge counts.

| k | Connected components | Unique edges | Average degree | k+2 | Condition |
| --- | ---: | ---: | ---: | ---: | --- |
| 1 | 56 | 244 | 1.627 | 3 | Pass |
| 2 | 1 | 466 | 3.107 | 4 | Pass |
| 3 | 1 | 679 | 4.527 | 5 | Pass |
| 4 | 1 | 885 | 5.900 | 6 | Pass |
| 5 | 1 | 1103 | 7.353 | 7 | Fail |
| 10 | 1 | 2075 | 13.833 | 12 | Fail |
| 15 | 1 | 3000 | 20.000 | 17 | Fail |

Table 1 gives k_min=2. The complete set of k satisfying the degree condition
is {1,2,3,4,297,298,299}. Consequently:

- The paper's literal global-maximum definition gives [2,299]. This interval
  also contains the k values at which the degree condition fails.
- Our explicitly added first-violation convention stops at k=5 and gives [2,4].
  It excludes every k=5 through 11, which the preceding experiment already
  showed to recover above 0.9967, including its one-sweep choice k=10.

This is a concrete coverage limitation of that local upper-bound convention.
At this graph-only stage, recovery of k=2,3,4 had not been evaluated.
The condition becomes true again near the complete graph, confirming that it
is not a monotone exclusion rule on this dataset.

`graph_bounds.R` generates `graph_bounds.csv` for all 299 graph sizes,
`graph_bounds_summary.csv` for both interval conventions, and
`graph_rule_passing_runs.csv` for the two passing runs. `graph_bounds.rds`
retains those tables, provenance, and six independent graph checks. The latter
compare truncation of the all-neighbor query with separate queries at k=2,4,5,
10,15,299 and verify agreement with the package-style graph after removing
parallel edges. Connectivity is monotone and the final graph has 44,850 edges,
as required for 300 nodes. The matrix SHA256 matches the prior fitted experiment.

Reproduce from the InferOrder root:

```sh
BASH_ENV=/dev/null bash experiments/isomap_elbo_screen_v040/run_r.sh experiments/isomap_elbo_screen_v040/graph_bounds.R
```

Only neighbor-graph diagnostics were computed; existing fits and timings were
preserved. Package source, defaults, and public pages were not changed.

## Subsequent fits of [2,4]: 2026-10-02

The user requested fitting the three smaller neighborhoods. All three recover
well: final absolute Spearman correlation is 0.996796, 0.996766, and 0.996788
for k=2,3,4. One-sweep ELBO selects k=4. This establishes successful candidates
inside [2,4] on the same dataset, despite exclusion of the earlier good starts
k=5 through 11. See [the follow-up report](narrow_range.md) for selection,
controls, complete endpoints, checks, and Figure 1. Generalization to other
datasets and interpretation of the paper's upper-bound convention remain open.
