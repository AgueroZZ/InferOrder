# Fixed-curve noise reduction: replicate 6, M=1, P=4

Specified on 2026-10-03 before fitting reduced-noise observations.

Use the previously inspected four-feature replicate 6 with all 200 true
positions, spline coefficients, signal scaling, and standard-normal errors
fixed. Set observations to `signal + noise_sd * standard_noise` at SD 0.25,
0.10, 0.05, and 0.01. The SD 0.25 input must match the parent bitwise; reuse
its three saved fits. No signals, sample positions, or noise realizations are
redrawn. Average dense-grid signal variance stays one, giving variance SNRs
16, 100, 400, and 10000.

At each level use automatic Isomap kmin, fixed k=15, and fixed k=10 on all
four features. Recompute kmin from each observed dataset. Use the frozen parent
MPCurver 0.4.0.9000 implementation, fit seed 202620026, fifty position bins,
quantile initialization, RW2, ridge zero, adaptive noise/precision/position
weights, normalized ELBO tolerance 1e-6, and continuation up to 10000 sweeps.
Observed variances are not standardized, and truth does not select fits.

Retain raw Isomap coordinates before binning, converged/budget-limited model
positions, objective traces, statuses, warnings, actual counts, and hashes.
Evaluate absolute Spearman to allow global reversal and average ranks for
ties. Report every level/method, including failures. This is one fixed case;
scores are descriptive, with no replication-based error bars.

Plot raw and final scores against noise SD. Compare feature 1 versus feature
4 across noise levels, colored by truth or raw automatic Isomap positions;
only global reflection is allowed for display. Both axes are observed features,
and Isomap still uses all four features. Save inputs, compact fits, plotted
data, and provenance. Redirect full fits to level-specific ignored paths;
preserve all parent fits, inputs, and scientific outputs. No package/default
or public workflowr changes.
