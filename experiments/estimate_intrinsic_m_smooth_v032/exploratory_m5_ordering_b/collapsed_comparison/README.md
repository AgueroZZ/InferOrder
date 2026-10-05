# Profiling the MPCurve curve block on ordering B

This exploratory experiment asks whether analytically profiling the Gaussian curve variational factor improves recovery or optimization speed for the PCA-initialized ordering B example. It uses the 300 observations and 12 B features from `main_M5_S4_r001`, fitting only one ordering. The three inputs retain the same true signal and realized noise, multiplied by 0, 0.5, and 1. Input files and their truth/PCA coordinates are in `../external_methods/inputs/`.

The profiled optimizer improves the attained ELBO in all three primary comparisons, but does not unfold the incorrect PCA ordering. Its terminal runtime is longer; time to a common objective level varies by noise setting. These results concern this implementation and initialization, not a general comparison of all collapsed optimizers.

## Objective and controls

The reference is MPCurver 0.3.4, source commit `15f2b0bbe5dfa61cd46da5160b2bc251e75a0475`, native PCA initialization, K=50, RW2, adaptive feature-specific noise variances and precision parameters, and adaptive grid probabilities. Frozen package and input provenance are exported in `*_reference.json`.

For responsibilities R and parameters theta, the alternative optimizes Lbar(R,theta) = max_q(U) L(R,q(U),theta). The Gaussian posterior precision is lambda_j Q + diag(colSums(R))/sigma_j^2. The curve mean and full covariance are recomputed analytically at every objective evaluation. Grid probabilities are also profiled as colMeans(R). Log noise variances and log precisions are optimized jointly with R, with the package's bounds of 1e-10 to 1e10. The intrinsic RW2 rank and pseudo-determinant convention match the package. This is the same variational objective with selected variables optimized out, not a different prior or a claim of exact marginal likelihood inference.

Native PCA responsibilities are hard assignments. Primary paired comparisons give both algorithms R0 = 0.99 R_PCA + 0.01/K to permit optimization in the simplex interior; this preserves the initial ordering but can alter the convergence basin. Original public-API default fits are reported separately. Half-noise sensitivity runs use epsilon=0.0001 and 0.1. A further half-noise comparison holds all noise, precision, and grid-probability parameters fixed at their initial values.

CAVI is reproduced in Python with the package update order, a 2000-sweep cap, and normalized ELBO increment threshold 1e-6. The main profiled method uses R_ik = a_ik^2 / sum_l a_il^2 with L-BFGS-B. This parameterization avoids the severe saturation observed with softmax logits. Optimization uses negative ELBO divided by n*d, ftol=1e-10, gtol=1e-8, maxiter=3000, maxfun=6000, maxls=40, maxcor=15. Signed square coordinates are bounded to [-40,40]. RW2 eigen-coordinates preserve the two null modes and improve numerical stability. Invalid Cholesky trial points are rejected and counted.

Both Python methods use the same posterior linear algebra and one BLAS thread. Primary runtimes are medians of three sequential runs, excluding validation and serialization. Iteration counts are not equivalent work units: one quasi-Newton iteration can evaluate the posterior multiple times. Stopping rules differ, so `summary.csv` also reports time to the paired CAVI final ELBO minus 0.0036 (n*d*1e-6), measured from the last repeat's recorded trace. This common-target timing is not a three-repeat median.

## Results

Absolute Spearman correlation uses truth versus posterior expected grid position, allowing reversal.

| Noise multiplier | Original package default rho | Paired CAVI rho | Profiled rho | CAVI ELBO | Profiled ELBO | CAVI median seconds | Profiled median seconds |
|---|---:|---:|---:|---:|---:|---:|---:|
| 0 | 0.1235 | 0.1262 | 0.0564 | 2247.829 | 3386.914 | 0.126 | 9.522 |
| 0.5 | 0.0971 | 0.0957 | 0.1002 | -1662.563 | -1466.001 | 0.146 | 4.461 |
| 1 | 0.0847 | 0.0868 | 0.1376 | -3622.659 | -3570.149 | 0.443 | 5.466 |

Half-noise CAVI takes 47 sweeps; the profiled optimizer takes 817 iterations and 844 objective/gradient evaluations. At the common CAVI objective target, recorded times are 0.142 versus 0.252 seconds. Across no noise and original noise, the respective common-target times are 0.123 versus 0.403 and 0.440 versus 0.333 seconds. Thus the original-noise profiled optimizer reaches that common target faster, despite taking longer to reach its own terminal result.

Original native R default median runtimes are 0.131, 0.191, and 0.599 seconds; their ELBO values are 2166.499, -1700.889, and -3615.323. These are reference results, not the matched-backend timing comparison.

At half noise, epsilon=0.0001 gives CAVI/profiled rho 0.0971/0.0980 and epsilon=0.1 gives 0.1407/0.1073. The failure to recover ordering is robust to these particular interiorization choices. Holding hyperparameters fixed gives rho 0.0695/0.0738; profiling does not repair that ordering either.

The softmax-logit implementation is retained as a numerical diagnostic. It stalls at substantially inferior ELBO values; in two cases L-BFGS reports success even though one CAVI sweep improves ELBO by over 800. It should not be interpreted as a successfully converged collapsed fit. The square-coordinate primary endpoints have one-sweep CAVI gains below 0.00007, smaller than the native 0.0036 stopping scale. This check supports local convergence but does not establish a global optimum.

Higher ELBO with persistently poor ordering shows that simply changing the optimization route is insufficient for this example. It does not establish that the global objective optimum prefers incorrect ordering, nor isolate the cause among prior, discretization, noise adaptation, and initialization.

## Verification and figures

`reference_checks.json` compares paired Python and R-package CAVI endpoints: maximum responsibility error < 1.1e-12 and ELBO error < 1.1e-9. Validation in each result JSON also checks the initial objective, first three native CAVI sweeps, an independent determinant expression, and directional gradients for all parameter blocks. `endpoint_audit.csv` reevaluates all six primary endpoints directly in the frozen R package; maximum objective discrepancy is 1.1e-9.

`ordering_comparison.png` shows true latent position versus estimated rank/n, with gray initial PCA and blue final positions, oriented toward positive truth correlation. Rows distinguish native default, paired CAVI, and profiled optimization. `convergence_half_noise.png` shows elapsed-time objective traces; its right panel restricts ELBO to [-1800,-1400]. Both images were visually inspected.

## Reproduction

From the InferOrder root, using the existing experiment environments:

```bash
bash experiments/estimate_intrinsic_m_smooth_v032/run_r.sh experiments/estimate_intrinsic_m_smooth_v032/exploratory_m5_ordering_b/collapsed_comparison/export_reference.R
```

Run `compare.py` with the Python interpreter in `../external_methods/.venv/bin/python` (Python 3.9, NumPy 1.26.4, SciPy 1.12.0). Positional case arguments are B_noise0, B_noise05, and B_noise1; each defaults to three repeats. For half-noise controls add `--fixed`, or `--epsilon 0.0001 --repeats 1`, or `--epsilon 0.1 --repeats 1`. Then run `report.py` with the same interpreter, followed by `audit_endpoints.R` through the R wrapper above. The result JSON files retain reference and optimizer-script SHA256 hashes, timings, stopping messages, gradient checks, and histories; NPZ files retain fitted states. No website or package implementation was changed.
