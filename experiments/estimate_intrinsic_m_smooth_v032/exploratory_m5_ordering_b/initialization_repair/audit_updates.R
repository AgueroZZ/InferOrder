# Read-only audit of the frozen installed implementation and conditional updates.
source('experiments/estimate_intrinsic_m_smooth_v032/exploratory_m5_ordering_b/initialization_repair/common.R')
ns <- asNamespace('MPCurver')
functions <- c('.structural_partition_update_gamma','.structural_partition_sweep',
 '.structural_partition_update_weights','.structural_partition_stabilize_gamma')
text <- unlist(lapply(functions,function(name)c(paste0('# ',name),deparse(get(name,envir=ns)),'')))
writeLines(text,file.path(repair_dir,'frozen_update_functions.R.txt'))
fits <- lapply(c('isomap10','isomap15','pca'),function(id)
 readRDS(file.path(repair_dir,'full_fits',paste0(id,'.rds'))))
names(fits) <- c('isomap10','isomap15','pca')
states <- lapply(fits,function(x)MPCurver:::.structural_partition_state_from_fit(x$fit))
slot <- readRDS(file.path(repair_dir,'results','isomap10.rds'))$slot
rows <- list()
for(id in names(states)) {
 bad <- states[[id]]
 # This is an oracle diagnostic of the coordinate update, not an estimator:
 # supply q(U) from a good fit while retaining this fit's noise and weights.
 for(source in unique(c(id,'isomap10'))) {
  donor <- states[[source]]
  posterior <- donor$q_u[[slot]]
  # Align donor grid orientation with the receiver before comparing jumps.
  donor_pos <- as.numeric(donor$gamma[[slot]] %*% grid)
  bad_pos <- as.numeric(bad$gamma[[slot]] %*% grid)
  if(cor(donor_pos,bad_pos,method='spearman') < 0) {
   posterior$m_mat <- posterior$m_mat[,50:1,drop=FALSE]
   posterior$sdiag_mat <- posterior$sdiag_mat[,50:1,drop=FALSE]
  }
  updated <- MPCurver:::.structural_partition_update_gamma(X=dataset$X,q_u=posterior,
   pi_vec=bad$position_pi[[slot]],sigma2=bad$sigma2,measurement_sd=bad$measurement_sd,
   weights=bad$pi_weights[,slot])
  pos <- as.numeric(updated %*% grid)
  rank_change <- abs(rank01(pos)-rank01(bad_pos))
  rows[[length(rows)+1]] <- data.frame(receiver=id,trajectory_source=source,
   old_rho=abs(cor(truth,bad_pos,method='spearman')),updated_rho=abs(cor(truth,pos,method='spearman')),
   mean_rank_change=mean(rank_change),fraction_moves_over_quarter=mean(rank_change>.25),
   min_position_probability=min(updated),mean_max_position_probability=mean(apply(updated,1,max)))
 }
}
write.csv(do.call(rbind,rows),file.path(repair_dir,'conditional_position_updates.csv'),row.names=FALSE)
feature_rows <- do.call(rbind,lapply(names(states),function(id) {
 state <- states[[id]]
 data.frame(initialization=id,feature=colnames(dataset$X)[B],anchor=B %in% anchor,
  true_noise_variance=.25,realized_noise_variance=dataset$realized_noise_variance[B],
  estimated_noise_variance=state$sigma2[B],lambda=state$lambda_mat[B,slot])
}))
write.csv(feature_rows,file.path(repair_dir,'noise_and_smoothness.csv'),row.names=FALSE)
print(do.call(rbind,rows),row.names=FALSE)
print(feature_rows[feature_rows$anchor,],row.names=FALSE)
