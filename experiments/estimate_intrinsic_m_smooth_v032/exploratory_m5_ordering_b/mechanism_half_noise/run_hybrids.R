# Experimental algorithm hybrids to distinguish penalty family from curve updates.
source('experiments/estimate_intrinsic_m_smooth_v032/auto_m_v034/common.R')
base <- file.path(study_dir,'exploratory_m5_ordering_b');out <- file.path(base,'mechanism_half_noise')
X <- as.matrix(read.csv(file.path(base,'external_methods/inputs/B_noise05_X.csv')))
pos <- read.csv(file.path(base,'external_methods/inputs/B_noise05_positions.csv'))
Q <- readRDS(file.path(out,'native_prefix_fits.rds'))[['0']]$Q_K
# Linear interpolation on 50 regular nodes with the exact MPCurve RW2 precision.
# Calibration targets conditional smoother df=5 at every iteration. It uses no truth.
make_rw2_smoother <- function(df=5) {
 cache <- new.env(parent=emptyenv())
 function(lambda,xj,...) {
  ord <- order(lambda); t <- lambda[ord]
  t <- (t-min(t))/diff(range(t))
  if(is.null(cache$t) || !identical(t,cache$t)) {
   K <- nrow(Q); u <- t*(K-1)+1; left <- pmin(floor(u),K-1L); w <- u-left
   B <- matrix(0,length(t),K)
   B[cbind(seq_along(t),left)] <- 1-w
   B[cbind(seq_along(t),left+1L)] <- w
   BtB <- crossprod(B)
   solve_at <- function(a) chol2inv(chol(BtB+exp(a)*Q))
   root <- uniroot(function(a)sum(solve_at(a)*BtB)-df,c(-15,25),tol=1e-8)$root
   inverse <- solve_at(root)
   cache$B <- B; cache$operator <- inverse %*% t(B); cache$t <- t
   cache$edf <- sum(inverse*BtB)
  }
  as.numeric(cache$B %*% (cache$operator %*% xj[ord]))
 }
}
rank_spline <- function(lambda,xj,...) {
 ord <- order(lambda); t <- rank(lambda[ord],ties.method='average')
 predict(smooth.spline(t,xj[ord],df=5),x=t)$y
}
scenarios <- list(rw2_smoother_df5=list(smoother=make_rw2_smoother(5)),
 spline_rank_spacing=list(smoother=rank_spline),
 spline_no_stretch=list(stretch=0))
results <- list(); rows <- list()
for(id in names(scenarios)) {
 warnings <- character()
 fit <- withCallingHandlers(do.call(princurve::principal_curve,
  c(list(x=X,maxit=1000L,thresh=1e-6),scenarios[[id]])),
  warning=function(w){warnings<<-c(warnings,conditionMessage(w));invokeRestart('muffleWarning')})
 results[[id]] <- list(fit=fit,final=fit$lambda,truth=pos$truth,warnings=warnings)
 rows[[id]] <- data.frame(setting=id,rho=abs(cor(pos$truth,fit$lambda,method='spearman')),
  iterations=fit$num_iterations,converged=fit$converged,warnings=length(warnings))
 print(rows[[id]],row.names=FALSE);flush.console()
}
saveRDS(results,file.path(out,'hybrid_fits.rds'))
write.csv(do.call(rbind,rows),file.path(out,'hybrid_summary.csv'),row.names=FALSE)
