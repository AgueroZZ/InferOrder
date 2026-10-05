source('experiments/estimate_intrinsic_m_smooth_v032/auto_m_v034/common.R')
out <- file.path(study_dir,'exploratory_m5_ordering_b','external_methods')
dir.create(file.path(out,'princurve_results'),showWarnings=FALSE)
cases <- read.csv(file.path(out,'cases.csv'))
results <- list()
for(id in cases$id) {
 X <- as.matrix(read.csv(file.path(out,'inputs',paste0(id,'_X.csv'))))
 positions <- read.csv(file.path(out,'inputs',paste0(id,'_positions.csv')))
 starts <- c('pca',if(grepl('^B_',id))c('isomap10','isomap15'))
 for(start_method in starts) {
  start_curve <- NULL
  if(start_method!='pca') {
   latent <- positions[[start_method]]
   start_curve <- vapply(seq_len(ncol(X)),function(j)
    predict(smooth.spline(latent,X[,j],df=5),x=sort(latent))$y,numeric(nrow(X)))
  }
  initial <- princurve::principal_curve(X,start=start_curve,maxit=0)$lambda
  for(setting in c('package_defaults','extended')) {
   settings <- if(setting=='package_defaults') list() else list(maxit=1000L,thresh=1e-6)
   warns <- character(); started <- proc.time()[['elapsed']]
   fit <- withCallingHandlers(do.call(princurve::principal_curve,c(list(x=X,start=start_curve),settings)),
    warning=function(w){warns<<-c(warns,conditionMessage(w));invokeRestart('muffleWarning')})
   tag <- paste(id,start_method,setting,sep='_')
   result <- list(id=id,start_method=start_method,setting=setting,
    initial=initial,final=fit$lambda,truth=positions$truth,fit=fit,
    package_version=as.character(packageVersion('princurve')),warnings=warns,
    elapsed_seconds=proc.time()[['elapsed']]-started,
    input_file_sha256=digest::digest(file=file.path(out,'inputs',paste0(id,'_X.csv')),algo='sha256'),
    start_curve=start_curve)
   saveRDS(result,file.path(out,'princurve_results',paste0(tag,'.rds')))
   results[[tag]] <- data.frame(id=id,start_method=start_method,setting=setting,
    initial_rho=abs(cor(positions$truth,initial,method='spearman')),
    final_rho=abs(cor(positions$truth,fit$lambda,method='spearman')),
    initial_final_rho=abs(cor(initial,fit$lambda,method='spearman')),
    converged=fit$converged,iterations=fit$num_iterations,distance=fit$dist,
    warning_count=length(warns),elapsed_seconds=result$elapsed_seconds)
   print(results[[tag]],row.names=FALSE);flush.console()
  }
 }
}
write.csv(do.call(rbind,results),file.path(out,'princurve_summary.csv'),row.names=FALSE)
writeLines(c(deparse(princurve::principal_curve),deparse(princurve:::smoother_functions$smooth_spline)),
 file.path(out,'princurve_source_snapshot.R.txt'))
