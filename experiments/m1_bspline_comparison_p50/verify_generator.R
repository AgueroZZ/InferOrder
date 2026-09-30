source('experiments/m1_bspline_comparison_p50/common.R')
invisible(loadNamespace('splines'))
for(replication in seq_len(design$replications)) {
  sim <- readRDS(sprintf('%s/inputs/rep%02d.rds',study,replication))
  parent_path <- sprintf('%s/inputs/rep%02d.rds',design$parent_study,replication)
  parent <- readRDS(parent_path)
  set.seed(sim$extra_seed)
  extra_coefficients <- matrix(rnorm(8L*38L),8L)
  extra_noise <- matrix(rnorm(design$N*38L),design$N)
  stopifnot(identical(sim$truth,parent$truth),
    identical(sim$coefficients,cbind(parent$coefficients,extra_coefficients)),
    identical(sim$noise,cbind(parent$noise,extra_noise)),
    identical(sim$parent_file_sha256,digest::digest(file=parent_path,algo='sha256')),
    isTRUE(all.equal(attr(sim$basis,'knots'),c(.2,.4,.6,.8))),
    abs(mean(sim$dense_signal^2)-1)<1e-12)
  reconstructed <- sweep(predict(sim$basis,sim$truth)%*%sim$coefficients,2L,sim$centers)/sim$signal_scale
  stopifnot(identical(unname(reconstructed),unname(sim$signal)))
}
writeLines('All 30 extensions exactly preserve parent sample positions, first 12 coefficient sets and noise; extra seeds regenerate all 38 added features; signals, basis knots, and unit average variance verified.',file.path(study,'generator_verification.txt'))
