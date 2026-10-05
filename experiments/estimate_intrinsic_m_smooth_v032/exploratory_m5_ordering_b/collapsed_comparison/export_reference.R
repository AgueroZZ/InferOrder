source('experiments/estimate_intrinsic_m_smooth_v032/auto_m_v034/common.R')
base <- file.path(study_dir,'exploratory_m5_ordering_b');out <- file.path(base,'collapsed_comparison')
for(id in c('B_noise0','B_noise05','B_noise1')) {
 X <- as.matrix(read.csv(file.path(base,'external_methods/inputs',paste0(id,'_X.csv'))))
 positions <- read.csv(file.path(base,'external_methods/inputs',paste0(id,'_positions.csv')))
 initial <- MPCurver::fit_mpcurve(X,method='PCA',intrinsic_dim=1,iter=0)$fit
 meta <- MPCurver:::.rw_precision_metadata(initial$Q_K,rw_q=2L)
 stages <- list();timings <- numeric(3)
 for(i in 1:3)stages[[i]] <- MPCurver::fit_mpcurve(X,method='PCA',intrinsic_dim=1,iter=i)$fit
 for(i in 1:3) {
  started <- proc.time()[['elapsed']]
  final <- MPCurver::fit_mpcurve(X,method='PCA',intrinsic_dim=1,iter=2000)$fit
  timings[i] <- proc.time()[['elapsed']]-started
 }
 pack <- function(f)list(R=f$gamma,sigma2=f$params$sigma2,lambda=f$lambda_vec,pi=f$params$pi,
  mean=f$posterior$mean,variance=f$posterior$var,elbo=tail(f$elbo_trace,1),iterations=f$iter)
 jsonlite::write_json(list(id=id,X=X,truth=positions$truth,pca=positions$pca,Q=initial$Q_K,
  rank=meta$rank,logdet_Q=meta$logdet,initial=pack(initial),stages=lapply(stages,pack),
  native_final=pack(final),native_trace=final$elbo_trace,native_timing_seconds=timings,
  input_sha256=digest::digest(X,algo='sha256'),package_commit=design$package_source_commit),
  file.path(out,paste0(id,'_reference.json')),digits=17,auto_unbox=TRUE)
 # Paired R-package check of the same 1% interiorization used for the primary optimizers.
 R <- .99*initial$gamma+.01/ncol(initial$gamma)
 attr(R,'cavi_skip_raw_preinit') <- TRUE
 warm <- MPCurver:::cavi(X,K=50,responsibilities_init=R,sigma2_init=initial$params$sigma2,
  lambda_init=initial$lambda_vec,position_prior_init=initial$params$pi,
  max_iter=2000,tol=1e-6,convergence='normalized')
 jsonlite::write_json(pack(warm),file.path(out,paste0(id,'_warm_reference.json')),digits=17,auto_unbox=TRUE)
}
