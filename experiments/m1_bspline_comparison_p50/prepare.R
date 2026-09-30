source('experiments/m1_bspline_comparison_p50/common.R')
# Extend each existing replicate with independent features, preserving its latent
# positions and the original coefficients/noise. Renormalize average signal variance.
simulate_replicate <- function(replication) {
  parent_path <- sprintf('%s/inputs/rep%02d.rds', design$parent_study, replication)
  parent <- readRDS(parent_path)
  extra_seed <- parent$seed + design$additional_seed_offset
  set.seed(extra_seed)
  extra_p <- design$P - ncol(parent$coefficients)
  coefficients <- cbind(parent$coefficients,
                        matrix(rnorm(design$spline_df*extra_p),design$spline_df))
  noise <- cbind(parent$noise,matrix(rnorm(design$N*extra_p),design$N))
  grid <- parent$grid
  basis_grid <- splines::bs(grid, df=design$spline_df, degree=design$degree,
                            intercept=TRUE, Boundary.knots=c(0,1))
  dense_signal <- basis_grid %*% coefficients
  centers <- colMeans(dense_signal)
  signal_scale <- sqrt(mean(sweep(dense_signal,2L,centers)^2))
  signal <- sweep(predict(basis_grid,parent$truth)%*%coefficients,2L,centers)/signal_scale
  list(replication=replication,seed=parent$seed,extra_seed=extra_seed,
       parent_file_sha256=digest::digest(file=parent_path,algo='sha256'),
       truth=parent$truth,coefficients=coefficients,centers=centers,
       signal_scale=signal_scale,grid=grid,
       dense_signal=sweep(dense_signal,2L,centers)/signal_scale,
       signal=signal,noise=noise,basis=basis_grid)
}
manifest <- list()
for (replication in seq_len(design$replications)) {
  sim <- simulate_replicate(replication)
  saveRDS(sim, sprintf('%s/inputs/rep%02d.rds',study,replication))
  for (noise in names(design$noise_sd)) {
    id <- sprintf('%s_r%02d',noise,replication)
    X <- sim$signal + design$noise_sd[[noise]] * sim$noise
    X <- scale(X, center=TRUE, scale=FALSE)
    pca <- prcomp(X, center=FALSE, scale.=FALSE)
    initial <- as.numeric(scale(pca$x[,1]))
    # Use population variance one, matching the latent scaling used in the pilot.
    initial <- initial / sqrt(mean(initial^2))
    start_curve <- tcrossprod(pca$x[,1],pca$rotation[,1])[order(pca$x[,1]),,drop=FALSE]
    input <- list(id=id, replication=replication, noise=noise, X=X,
                  truth=sim$truth, initial=initial, start_curve=start_curve,
                  pca_variance_fraction=pca$sdev[1]^2/sum(pca$sdev^2))
    saveRDS(input,file.path(study,'inputs',paste0(id,'.rds')))
    write.csv(X,file.path(study,'inputs',paste0(id,'_X.csv')),row.names=FALSE)
    write.csv(data.frame(truth=sim$truth,pca=initial),
              file.path(study,'inputs',paste0(id,'_positions.csv')),row.names=FALSE)
    manifest[[length(manifest)+1L]] <- data.frame(id=id,replication=replication,noise=noise,
      noise_sd=design$noise_sd[[noise]], variance_snr=1/design$noise_sd[[noise]]^2,
      seed=sim$seed, pca_rho=recovery(sim$truth,initial),
      input_sha256=digest::digest(file=file.path(study,'inputs',paste0(id,'_X.csv')),algo='sha256'))
  }
}
write.csv(do.call(rbind,manifest),file.path(study,'manifest.csv'),row.names=FALSE)
saveRDS(design,file.path(study,'design.rds'))
writeLines(capture.output(dput(design)),file.path(study,'design.txt'))
writeLines(capture.output(sessionInfo()),file.path(study,'R_session.txt'))
writeLines(c(deparse(MPCurver:::.cavi_build_from_ordering),deparse(princurve::principal_curve)),
           file.path(study,'initialization_source.R.txt'))
