# Replay random generation independently from the recorded replicate seeds.
source('experiments/m1_bspline_comparison/common.R')
loadNamespace('splines')
for (replication in seq_len(design$replications)) {
  sim <- readRDS(sprintf('%s/inputs/rep%02d.rds',study,replication))
  set.seed(sim$seed)
  truth <- runif(design$N)
  coefficients <- matrix(rnorm(design$spline_df*design$P),design$spline_df)
  noise <- matrix(rnorm(design$N*design$P),design$N)
  stopifnot(identical(truth,sim$truth),identical(coefficients,sim$coefficients),
            identical(noise,sim$noise),
            isTRUE(all.equal(attr(sim$basis,'knots'),c(.2,.4,.6,.8))),
            abs(mean(sim$dense_signal^2)-1)<1e-12)
  reconstructed <- sweep(predict(sim$basis,truth)%*%coefficients,2L,sim$centers)/sim$signal_scale
  stopifnot(identical(unname(reconstructed),unname(sim$signal)))
}
writeLines('All 30 seeds exactly regenerate sample positions, spline coefficients, and noise; spline knots and unit average signal variance verified.',
           file.path(study,'generator_verification.txt'))
