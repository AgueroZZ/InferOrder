# Benchmark the public single-ordering API with PCA and native initialization.
source('experiments/estimate_intrinsic_m_smooth_v032/auto_m_v034/common.R')
base <- file.path(study_dir,'exploratory_m5_ordering_b')
out <- file.path(base,'single_order_benchmark')
for(id in c('B_noise0','B_noise05','B_noise1')) {
 input <- file.path(base,'external_methods','inputs',paste0(id,'_X.csv'))
 X <- as.matrix(read.csv(input))
 positions <- read.csv(file.path(base,'external_methods','inputs',paste0(id,'_positions.csv')))
 for(setting in c('package_stopping','extended_stopping')) {
  args <- list(X=X,algorithm='cavi',method='PCA',intrinsic_dim=1L)
  if(setting=='extended_stopping') args$iter <- 2000L
  set.seed(20260929); warnings <- character()
  fit <- withCallingHandlers(do.call(MPCurver::fit_mpcurve,args),
   warning=function(w){warnings<<-c(warnings,conditionMessage(w));invokeRestart('muffleWarning')})
  stopifnot(fit$intrinsic_dim==1L,ncol(fit$data)==12L,is.matrix(fit$gamma))
  final <- as.numeric(fit$gamma %*% seq(0,1,length.out=ncol(fit$gamma)))
  result <- list(id=id,setting=setting,initial=positions$pca,truth=positions$truth,
   final=final,fit=fit,warnings=warnings,seed=20260929,
   package_version=as.character(packageVersion('MPCurver')),
   package_commit=design$package_source_commit,arguments=args[names(args)!='X'],
   input_path=input,input_sha256=digest::digest(file=input,algo='sha256'),
   script_sha256=digest::digest(file=file.path(out,'fit_mpcurve.R'),algo='sha256'))
  saveRDS(result,file.path(out,paste0(id,'_mpcurve_',setting,'.rds')))
  print(data.frame(id=id,setting=setting,rho=abs(cor(positions$truth,final,method='spearman')),
   iterations=fit$fit$iter,converged=fit$converged,K=ncol(fit$gamma),warnings=length(warnings)))
 }
}
writeLines(capture.output(formals(MPCurver::fit_mpcurve)),file.path(out,'mpcurve_public_defaults.txt'))
