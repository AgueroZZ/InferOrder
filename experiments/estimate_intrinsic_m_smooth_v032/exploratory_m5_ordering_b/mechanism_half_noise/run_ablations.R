# Targeted mechanism diagnostics, all using the same PCA initialization.
source('experiments/estimate_intrinsic_m_smooth_v032/auto_m_v034/common.R')
base <- file.path(study_dir,'exploratory_m5_ordering_b')
out <- file.path(base,'mechanism_half_noise')
X <- as.matrix(read.csv(file.path(base,'external_methods/inputs/B_noise05_X.csv')))
pos <- read.csv(file.path(base,'external_methods/inputs/B_noise05_positions.csv'))
initial <- readRDS(file.path(out,'native_prefix_fits.rds'))[['0']]
Q <- initial$Q_K; Nk <- colSums(initial$gamma)
calibrate <- function(nk,sigma2,df) {
 diagonal <- nk/sigma2
 f <- function(loglambda) sum(diag(chol2inv(chol(diag(diagonal)+exp(loglambda)*Q)))*diagonal)-df
 exp(uniroot(f,c(-20,30),tol=1e-8)$root)
}
lambda5 <- vapply(initial$params$sigma2,function(s)calibrate(Nk,s,5),numeric(1))
scenarios <- list(native=list(),smooth_start_df5=list(lambda=lambda5),
 fixed_initial_df5=list(lambda=lambda5,fix_lambda=TRUE),
 fixed_lambda_100=list(lambda=100,fix_lambda=TRUE),
 fixed_lambda_1000=list(lambda=1000,fix_lambda=TRUE),
 fixed_lambda_10000=list(lambda=10000,fix_lambda=TRUE),
 fixed_lambda_100000=list(lambda=100000,fix_lambda=TRUE),
 known_noise=list(S=matrix(.25,nrow(X),ncol(X))),
 fixed_uniform_positions=list(position_prior='fixed'),
 equal_width_init=list(discretization='equal'))
rows <- list(); results <- list()
for(id in names(scenarios)) {
 set.seed(20260929); warnings <- character()
 args <- c(list(X=X,method='PCA',intrinsic_dim=1,iter=2000L),scenarios[[id]])
 fit <- withCallingHandlers(do.call(MPCurver::fit_mpcurve,args)$fit,
  warning=function(w){warnings<<-c(warnings,conditionMessage(w));invokeRestart('muffleWarning')})
 final <- as.numeric(fit$gamma %*% seq(0,1,length.out=ncol(fit$gamma)))
 variances <- if(is.null(fit$measurement_sd))fit$params$sigma2 else rep(.0625,ncol(X))
 edf <- vapply(seq_len(ncol(X)),function(j)sum(diag(fit$posterior$cov[[j]])*colSums(fit$gamma)/variances[j]),numeric(1))
 result <- list(id=id,method='MPCurve',fit=fit,final=final,truth=pos$truth,warnings=warnings,arguments=scenarios[[id]],edf=edf)
 results[[id]] <- result
 rows[[id]] <- data.frame(method='MPCurve',setting=id,rho=abs(cor(pos$truth,final,method='spearman')),
  iterations=fit$iter,converged=fit$converged,min_edf=min(edf),max_edf=max(edf),warnings=length(warnings))
 print(rows[[id]],row.names=FALSE);flush.console()
}
for(df in c(3,5,10,20)) {
 id <- paste0('pcurve_df',df)
 warnings <- character()
 fit <- withCallingHandlers(princurve::principal_curve(X,df=df,maxit=1000L,thresh=1e-6),
  warning=function(w){warnings<<-c(warnings,conditionMessage(w));invokeRestart('muffleWarning')})
 results[[id]] <- list(id=id,method='P-curve',fit=fit,final=fit$lambda,truth=pos$truth,warnings=warnings,df=df)
 rows[[id]] <- data.frame(method='P-curve',setting=id,rho=abs(cor(pos$truth,fit$lambda,method='spearman')),
  iterations=fit$num_iterations,converged=fit$converged,min_edf=df,max_edf=df,warnings=length(warnings))
 print(rows[[id]],row.names=FALSE);flush.console()
}
saveRDS(list(results=results,calibrated_lambda=lambda5,input_sha256=digest::digest(X,algo='sha256'),
 package_commit=design$package_source_commit,seed=20260929),file.path(out,'ablation_fits.rds'))
write.csv(do.call(rbind,rows),file.path(out,'ablation_summary.csv'),row.names=FALSE)
