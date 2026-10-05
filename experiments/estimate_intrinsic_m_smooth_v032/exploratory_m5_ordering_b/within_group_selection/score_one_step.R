# Score each initialization using only the feature group's single-ordering model.
source('experiments/estimate_intrinsic_m_smooth_v032/auto_m_v034/common.R')
base <- file.path(study_dir,'exploratory_m5_ordering_b')
out <- file.path(base,'within_group_selection')
d <- load_fixed_dataset(manifest[manifest$id=='main_M5_S4_r001',])
b <- which(d$true_assign=='B')
paths <- c(file.path(base,'sensitivity_fits',sprintf('setting_%02d.rds',1:15)),
 file.path(base,'pca_fits',sprintf('setting_%02d.rds',1:3)))
rows <- list(); saved_results <- list()
for(i in seq_along(paths)) {
 reference <- readRDS(paths[i]); if(reference$status!='converged') next
 X <- d$signal[,b]+reference$setting$noise_scale*(d$X[,b]-d$signal[,b])
 set.seed(d$row$seed+design$seed_offset)
 initial <- if(i>15L) MPCurver:::PCA_ordering(X)$t else
  MPCurver:::isomap_ordering(X,k=reference$setting$k)$t
 warnings <- character()
 started <- proc.time()[['elapsed']]
 fit <- withCallingHandlers(MPCurver:::.cavi_build_from_ordering(
  X=X,ordering_vec=initial,S=NULL,K=50L,rw_q=2L,ridge=0,lambda_init=1,
  max_iter=1L,tol=1e-6,discretization='quantile',strict_K=TRUE),
  warning=function(w){warnings<<-c(warnings,conditionMessage(w));invokeRestart('muffleWarning')})
 elapsed <- proc.time()[['elapsed']]-started
 stopifnot(fit$iter==1L,length(fit$elbo_trace)==2L,ncol(fit$data)==12L)
 label <- if(i>15L) 'PCA' else paste0('k=',reference$setting$k)
 position <- as.numeric(fit$gamma %*% seq(0,1,length.out=50))
 row <- data.frame(noise_scale=reference$setting$noise_scale,method=label,
  initial_group_elbo=fit$elbo_trace[1],one_step_group_elbo=fit$elbo_trace[2],
  one_step_group_rho=abs(cor(d$latent_positions[,2],position,method='spearman')),
  converged_joint_elbo=tail(reference$objective,1),
  converged_joint_rho=abs(cor(reference$truth,reference$final,method='spearman')),
  one_step_seconds=elapsed,warnings=length(warnings))
 rows[[length(rows)+1L]] <- row
 saved_results[[length(saved_results)+1L]] <- list(summary=row,initial_position=initial,
  position=position,elbo_trace=fit$elbo_trace,params=fit$params,lambda_vec=fit$lambda_vec,
  input_hash=digest::digest(X,algo='sha256'),dataset_hash=d$input_hash,
  seed=d$row$seed+design$seed_offset,package_commit=design$package_source_commit,
  feature_indices=b,feature_names=colnames(X),warnings=warnings)
}
summary <- do.call(rbind,rows)
write.csv(summary,file.path(out,'one_step_scores.csv'),row.names=FALSE)
saveRDS(saved_results,file.path(out,'one_step_results.rds'))
selection <- do.call(rbind,lapply(split(summary,summary$noise_scale),function(x) {
 best <- which.max(x$one_step_group_elbo)
 data.frame(noise_scale=x$noise_scale[1],selected=x$method[best],
  local_score=x$one_step_group_elbo[best],final_joint_rho=x$converged_joint_rho[best],
  final_joint_elbo_loss=max(x$converged_joint_elbo)-x$converged_joint_elbo[best])
}))
write.csv(selection,file.path(out,'selection.csv'),row.names=FALSE)
print(summary,row.names=FALSE);print(selection,row.names=FALSE)
