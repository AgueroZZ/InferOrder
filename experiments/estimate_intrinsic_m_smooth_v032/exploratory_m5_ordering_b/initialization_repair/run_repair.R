source('experiments/estimate_intrinsic_m_smooth_v032/exploratory_m5_ordering_b/initialization_repair/common.R')
indices <- as.integer(commandArgs(trailingOnly=TRUE))
stopifnot(length(indices)>0,all(indices %in% seq_len(nrow(scenarios))))
initialization <- make_initial_fits()
slot <- which(vapply(initialization$init_info,function(x)setequal(x$feature_idx,B),logical(1)))
stopifnot(length(slot)==1)
# Trace is process-local and only observes return values; package files and
# the state returned to the algorithm are unchanged.
record_sweep <- function(state, temperature) {
 count <<- count+1L
 pos <- as.numeric(state$gamma[[slot]] %*% grid)
 measurement <- metrics(pos)
 prediction <- state$gamma[[slot]] %*% t(state$q_u[[slot]]$m_mat[B,,drop=FALSE])
 trajectory[[count]] <<- data.frame(sweep=count,temperature=temperature,
  rho=measurement['rho'],pair_error=measurement['pair_error'],mean_rank_error=measurement['mean_rank_error'],
  objective=state$objective,mean_max_probability=mean(apply(state$gamma[[slot]],1,max)),
  mean_entropy=mean(-rowSums(state$gamma[[slot]]*log(state$gamma[[slot]]))),
  anchor_sigma2=state$sigma2[anchor],mean_B_sigma2=mean(state$sigma2[B]),
  signal_rmse=sqrt(mean((prediction-dataset$signal[,B])^2)),
  observed_rmse=sqrt(mean((prediction-dataset$X[,B])^2)))
 if(count %in% c(1,5,10,25,50,100,200,500,1000)) snapshots[[as.character(count)]] <<- pos
}
trace('.structural_partition_sweep',exit=quote(.GlobalEnv$record_sweep(returnValue(),T_now)),
 where=asNamespace('MPCurver'),print=FALSE)
for(i in indices) {
 row <- scenarios[i,]; path <- file.path(repair_dir,'results',paste0(row$id,'.rds'))
 if(file.exists(path)) next
 raw <- make_ordering(row)
 sub <- MPCurver:::.cavi_build_from_ordering(X=dataset$X[,B],ordering_vec=raw,
  K=50L,rw_q=2L,ridge=0,max_iter=2L,tol=1e-6,discretization='quantile',strict_K=TRUE)
 fits <- initialization$fits
 fits[[slot]] <- MPCurver:::cavi(X=dataset$X,K=50L,responsibilities_init=sub$gamma,
  position_prior_init=colMeans(sub$gamma),rw_q=2L,ridge=0,max_iter=0L,
  convergence='relative',verbose=FALSE)
 warm <- as.numeric(fits[[slot]]$gamma %*% grid)
 count <- 0L; trajectory <- list(); snapshots <- list(raw=raw,warm=warm); warnings <- character()
 set.seed(fit_seed)
 fit <- withCallingHandlers({
  f <- MPCurver::fit_mpcurve(X=dataset$X,intrinsic_dim=5L,algorithm='cavi',
   method='isomap',K=50L,rw_q=2L,ridge=0,fits_init=fits,
   partition_prior='adaptive',position_prior='adaptive',lambda=1,fix_lambda=FALSE,
   discretization='quantile',num_cores=1L,iter=1500L,tol=1e-6,convergence='normalized',
   T_start=5,T_end=1,n_outer=25L,inner_iter=1L,max_converge_iter=1500L,
   tol_outer=1e-6,verbose=FALSE)
  while(!isTRUE(f$fit$converged) && f$fit$iter<9999L) f <- MPCurver::do_mpcurve(f,
   iter=min(1500L,9999L-f$fit$iter),tol=1e-6,tol_outer=1e-6,convergence='normalized',verbose=FALSE)
  f
 },warning=function(w){warnings<<-c(warnings,conditionMessage(w));invokeRestart('muffleWarning')})
 final <- as.numeric(fit$fit$gamma[[slot]] %*% grid)
 snapshots$final <- final
 result <- list(scenario=row,raw=raw,warm=warm,final=final,truth=truth,
  raw_metrics=metrics(raw),warm_metrics=metrics(warm),final_metrics=metrics(final),
  trajectory=do.call(rbind,trajectory),snapshots=snapshots,
  converged=fit$fit$converged,iterations=fit$fit$iter,objective=tail(fit$fit$objective_history,1),
  ARI=mclust::adjustedRandIndex(dataset$true_assign,fit$fit$assign),warnings=warnings,
  sigma2=fit$params$sigma2,lambda_mat=fit$lambda_mat,position_pi=fit$fit$position_pi[[slot]],
  gamma=fit$fit$gamma[[slot]],conditional_mean=fit$fit$conditional_posterior$mean[[slot]],
  slot=slot,input_hash=dataset$input_hash,fit_seed=fit_seed,
  package_version=as.character(packageVersion('MPCurver')),package_commit=design$package_source_commit,
  script_hash=digest::digest(file=file.path(repair_dir,'run_repair.R'),algo='sha256'),
  common_hash=digest::digest(file=file.path(repair_dir,'common.R'),algo='sha256'))
 atomic_save(result,path,compress=TRUE)
 if(row$id %in% c('oracle','isomap10','isomap15','pca')) atomic_save(fit,
  file.path(repair_dir,'full_fits',paste0(row$id,'.rds')),compress=TRUE)
 cat(sprintf('%s: %.4f -> %.4f -> %.4f, pair error %.4f -> %.4f; converged=%s ARI=%.3f\n',
  row$id,result$raw_metrics['rho'],result$warm_metrics['rho'],result$final_metrics['rho'],
  result$raw_metrics['pair_error'],result$final_metrics['pair_error'],result$converged,result$ARI))
 flush.console()
}
