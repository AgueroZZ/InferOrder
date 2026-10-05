# Single-ordering MPCurve and principal-curve iteration diagnostics.
source('experiments/estimate_intrinsic_m_smooth_v032/auto_m_v034/common.R')
out <- file.path(study_dir,'exploratory_m5_ordering_b','external_methods')
dir.create(file.path(out,'mpcurve_single_results'),showWarnings=FALSE)
summary <- list(); paths <- list()
for(id in c('B_noise0','B_noise05','B_noise1')) {
 X <- as.matrix(read.csv(file.path(out,'inputs',paste0(id,'_X.csv'))))
 positions <- read.csv(file.path(out,'inputs',paste0(id,'_positions.csv')))
 for(start in c('pca','isomap10','isomap15')) {
  set.seed(20260929); warnings <- character()
  fit <- withCallingHandlers(MPCurver:::.cavi_build_from_ordering(
   X=X,ordering_vec=positions[[start]],S=NULL,K=50L,rw_q=2L,ridge=0,
   lambda_init=1,max_iter=2000L,tol=1e-6,discretization='quantile',strict_K=TRUE),
   warning=function(w){warnings<<-c(warnings,conditionMessage(w));invokeRestart('muffleWarning')})
  final <- as.numeric(fit$gamma %*% seq(0,1,length.out=50))
  result <- list(id=id,start=start,initial=positions[[start]],final=final,truth=positions$truth,
   fit=fit,warnings=warnings,package_version=as.character(packageVersion('MPCurver')),
   package_commit=design$package_source_commit,seed=20260929,
   input_file_sha256=digest::digest(file=file.path(out,'inputs',paste0(id,'_X.csv')),algo='sha256'))
  saveRDS(result,file.path(out,'mpcurve_single_results',paste0(id,'_',start,'.rds')))
  summary[[length(summary)+1L]] <- data.frame(id=id,start=start,
   initial_rho=abs(cor(positions$truth,positions[[start]],method='spearman')),
   final_rho=abs(cor(positions$truth,final,method='spearman')),
   converged=fit$converged,iterations=fit$iter,objective=tail(fit$elbo_trace,1),warnings=length(warnings))
 }
 for(start in c('pca',if(id=='B_noise1')'isomap15')) {
  reference <- readRDS(file.path(out,'princurve_results',paste0(id,'_',start,'_extended.rds')))
  for(budget in c(0,1,2,5,10,15,20,25,30,40,50,70)) {
   fit <- princurve::principal_curve(X,start=reference$start_curve,maxit=budget,thresh=1e-6)
   paths[[length(paths)+1L]] <- data.frame(id=id,start=start,budget=budget,
    iteration=fit$num_iterations,converged=fit$converged,distance=sum((X-fit$s)^2),package_distance=fit$dist,
    rho=abs(cor(positions$truth,fit$lambda,method='spearman')))
  }
 }
}
write.csv(do.call(rbind,summary),file.path(out,'mpcurve_single_summary.csv'),row.names=FALSE)
write.csv(do.call(rbind,paths),file.path(out,'princurve_iteration_paths.csv'),row.names=FALSE)
print(do.call(rbind,summary),row.names=FALSE)
writeLines(capture.output(sessionInfo()),file.path(out,'R_session.txt'))
