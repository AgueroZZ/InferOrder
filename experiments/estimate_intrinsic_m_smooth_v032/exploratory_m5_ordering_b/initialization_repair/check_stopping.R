# Continue the original converged solutions with stopping disabled for 500 sweeps.
source('experiments/estimate_intrinsic_m_smooth_v032/exploratory_m5_ordering_b/initialization_repair/common.R')
ids <- commandArgs(trailingOnly=TRUE)
stopifnot(length(ids)>0,all(ids %in% c('isomap10','isomap15','pca')))
for(id in ids) {
 path <- file.path(repair_dir,'results',paste0(id,'_forced_continuation.rds'))
 if(file.exists(path)) next
 fit <- readRDS(file.path(repair_dir,'full_fits',paste0(id,'.rds')))
 result <- readRDS(file.path(repair_dir,'results',paste0(id,'.rds')))
 state <- MPCurver:::.structural_partition_state_from_fit(fit$fit)
 slot <- result$slot
 start <- state$objective; initial_position <- as.numeric(state$gamma[[slot]] %*% grid)
 records <- list(); positions <- list(initial=initial_position)
 for(j in 0:500) {
  if(j>0) state <- MPCurver:::.structural_partition_sweep(state,T_now=1)
  pos <- as.numeric(state$gamma[[slot]] %*% grid)
  measurement <- metrics(pos)
  records[[j+1]] <- data.frame(extra_sweeps=j,objective=state$objective,
   rho=measurement['rho'],pair_error=measurement['pair_error'],
   mean_rank_error=measurement['mean_rank_error'],
   mean_rank_change=mean(abs(rank01(align(pos))-rank01(align(initial_position)))))
  if(j %in% c(1,10,50,100,200,500)) positions[[as.character(j)]] <- pos
  if(j %% 100 == 0) {
   cat(sprintf('%s +%d: rho %.6f, objective change %.6f\n',id,j,measurement['rho'],state$objective-start));flush.console()
  }
 }
 saveRDS(list(id=id,trajectory=do.call(rbind,records),positions=positions,
  original_stop_tolerance=1e-6,stopping_disabled=TRUE,additional_sweeps=500,
  input_hash=dataset$input_hash,package_commit=design$package_source_commit,
  final_state_metrics=list(sigma2=state$sigma2,lambda_mat=state$lambda_mat,
   position_pi=state$position_pi[[slot]],gamma=state$gamma[[slot]])),path)
}
