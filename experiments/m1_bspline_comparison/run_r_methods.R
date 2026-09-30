source('experiments/m1_bspline_comparison/common.R')
manifest <- read.csv(file.path(study,'manifest.csv'))
args <- commandArgs(TRUE)
if (length(args)) manifest <- manifest[manifest$replication %in% as.integer(args),]
for (id in manifest$id) {
  d <- readRDS(file.path(study,'inputs',paste0(id,'.rds')))
  for (method in c('MPCurve','PCurve')) {
    destination <- file.path(study,'results',paste0(id,'_',method,'.rds'))
    if (file.exists(destination)) next
    warnings <- character(); started <- proc.time()[['elapsed']]
    result <- tryCatch(withCallingHandlers({
      if (method == 'MPCurve') {
        initial_fit <- MPCurver:::.cavi_build_from_ordering(X=d$X,ordering_vec=d$initial,
          S=NULL,K=design$K,rw_q=2L,ridge=0,lambda_init=1,max_iter=0L,
          tol=design$mp_tol,discretization='quantile',strict_K=TRUE)
        fit <- MPCurver:::.cavi_build_from_ordering(X=d$X,ordering_vec=d$initial,
          S=NULL,K=design$K,rw_q=2L,ridge=0,lambda_init=1,max_iter=design$mp_max_iter,
          tol=design$mp_tol,discretization='quantile',strict_K=TRUE)
        position <- as.numeric(fit$gamma %*% seq(0,1,length.out=design$K))
        initial <- as.numeric(initial_fit$gamma %*% seq(0,1,length.out=design$K))
        iterations <- fit$iter
      } else {
        initial_fit <- princurve::principal_curve(d$X,start=d$start_curve,maxit=0)
        fit <- princurve::principal_curve(d$X,start=d$start_curve,
          maxit=design$pc_max_iter,thresh=design$pc_tol)
        initial <- initial_fit$lambda; position <- fit$lambda
        iterations <- fit$num_iterations
      }
      list(fit=fit,initial=initial,final=position,converged=isTRUE(fit$converged),
           iterations=iterations,status=if(isTRUE(fit$converged)) 'converged' else 'iteration limit')
    },warning=function(w){warnings<<-c(warnings,conditionMessage(w));invokeRestart('muffleWarning')}),
    error=function(e) list(status=paste('error:',conditionMessage(e)),converged=FALSE,
                          initial=rep(NA_real_,design$N),final=rep(NA_real_,design$N),iterations=NA_integer_))
    result$id <- id; result$method <- method; result$warnings <- warnings
    result$elapsed_seconds <- proc.time()[['elapsed']]-started
    result$rho <- recovery(d$truth,result$final)
    result$tau <- recovery(d$truth,result$final,'kendall')
    result$input_sha256 <- manifest$input_sha256[match(id,manifest$id)]
    result$script_sha256 <- digest::digest(file=file.path(study,'run_r_methods.R'),algo='sha256')
    saveRDS(result,destination)
    cat(id,method,result$status,result$rho,result$elapsed_seconds,'seconds\n');flush.console()
  }
}
