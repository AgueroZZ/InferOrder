# Replay only the original annealing phase and evaluate ordinary T=1 ELBO.
source('experiments/estimate_intrinsic_m_smooth_v032/exploratory_m5_ordering_b/initialization_repair/common.R')
base <- file.path(study_dir,'exploratory_m5_ordering_b')
out <- file.path(base,'early_selection'); dir.create(file.path(out,'results'),showWarnings=FALSE)
paths <- c(file.path(base,'sensitivity_fits',sprintf('setting_%02d.rds',1:15)),
 file.path(base,'pca_fits',sprintf('setting_%02d.rds',1:3)))
indices <- as.integer(commandArgs(trailingOnly=TRUE))
initialization <- make_initial_fits()
slot <- which(vapply(initialization$init_info,function(x)setequal(x$feature_idx,B),logical(1)))
record_early <- function(state,temperature) {
 assignment <- MPCurver:::.structural_partition_assignment_info(state$pi_weights,1,state$control)
 elbo <- MPCurver:::.structural_partition_objective(state$gamma,state$position_pi,
  state$pi_weights,state$local_blocks,assignment)$objective
 entropy <- -sum(state$pi_weights*log(pmax(state$pi_weights,.Machine$double.xmin)))
 stopifnot(abs(elbo-(state$objective-(temperature-1)*entropy))<1e-7)
 history[[length(history)+1L]] <<- data.frame(sweep=length(history)+1L,temperature=temperature,
  annealed_objective=state$objective,standard_elbo=elbo,assignment_entropy=entropy)
}
trace('.structural_partition_sweep',exit=quote(.GlobalEnv$record_early(returnValue(),T_now)),
 where=asNamespace('MPCurver'),print=FALSE)
for(i in indices) {
 path <- file.path(out,'results',sprintf('candidate_%02d.rds',i))
 if(file.exists(path)) next
 saved <- readRDS(paths[i]); if(saved$status!='converged') next
 is_pca <- i>15L
 label <- if(is_pca) 'PCA' else paste0('k=',saved$setting$k)
 X <- dataset$X
 X[,B] <- dataset$signal[,B]+saved$setting$noise_scale*(dataset$X[,B]-dataset$signal[,B])
 stopifnot(identical(digest::digest(X,algo='sha256'),saved$modified_input_hash))
 set.seed(fit_seed)
 raw <- if(is_pca) MPCurver:::PCA_ordering(X[,B])$t else
  MPCurver:::isomap_ordering(X[,B],k=saved$setting$k)$t
 sub <- MPCurver:::.cavi_build_from_ordering(X=X[,B],ordering_vec=raw,K=50L,
  rw_q=2L,ridge=0,max_iter=2L,tol=1e-6,discretization='quantile',strict_K=TRUE)
 fits <- lapply(1:5,function(m) {
  gamma <- if(m==slot) sub$gamma else initialization$fits[[m]]$gamma
  MPCurver:::cavi(X=X,K=50L,responsibilities_init=gamma,position_prior_init=colMeans(gamma),
   rw_q=2L,ridge=0,max_iter=0L,convergence='relative',verbose=FALSE)
 })
 history <- list()
 set.seed(fit_seed)
 started <- proc.time()[['elapsed']]
 f <- MPCurver::fit_mpcurve(X=X,intrinsic_dim=5L,algorithm='cavi',method='isomap',
  K=50L,rw_q=2L,ridge=0,fits_init=fits,partition_prior='adaptive',position_prior='adaptive',
  lambda=1,fix_lambda=FALSE,discretization='quantile',num_cores=1L,iter=1500L,
  tol=1e-6,convergence='normalized',T_start=5,T_end=1,n_outer=25L,inner_iter=1L,
  max_converge_iter=0L,tol_outer=1e-6,verbose=FALSE)
 history <- do.call(rbind,history)
 stopifnot(nrow(history)==25L,
  max(abs(history$annealed_objective-saved$objective[2:26]))<1e-5)
 result <- list(method=label,noise_scale=saved$setting$noise_scale,history=history,
  final_objective=tail(saved$objective,1),final_rho=abs(cor(saved$truth,saved$final,method='spearman')),
  final_sweeps=length(saved$objective)-1L,input_hash=saved$modified_input_hash,
  seed=fit_seed,package_commit=design$package_source_commit,
  replay_seconds=proc.time()[['elapsed']]-started)
 atomic_save(result,path,compress=TRUE)
 cat(sprintf('noise=%g %s: first standard ELBO %.3f; final %.3f; replay matches\n',
  result$noise_scale,label,history$standard_elbo[1],result$final_objective));flush.console()
}
