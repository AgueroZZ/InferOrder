source('experiments/estimate_intrinsic_m_smooth_v032/auto_m_v034/common.R')
base <- file.path(study_dir,'exploratory_m5_ordering_b');out <- file.path(base,'gp_kernel_sanity')
files <- list.files(out,pattern='_audit_input.json$',full.names=TRUE)
results <- lapply(files,function(path) {
 tag <- sub('_audit_input.json$','',basename(path));case <- sub('_swap.*','',tag)
 ref <- jsonlite::fromJSON(file.path(base,'collapsed_comparison',paste0(case,'_reference.json')))
 input <- jsonlite::fromJSON(path);expected <- jsonlite::fromJSON(file.path(out,paste0(tag,'.json')))
 # Only this local function copy resolves the GP precision constructor differently.
 local_cavi <- MPCurver:::cavi
 env <- new.env(parent=environment(local_cavi))
 env$make_random_walk_precision <- function(...) input$Q
 environment(local_cavi) <- env
 R <- ref$initial$R;attr(R,'cavi_skip_raw_preinit') <- TRUE
 fit <- local_cavi(ref$X,K=50,responsibilities_init=R,sigma2_init=ref$initial$sigma2,
  lambda_init=rep(1,ncol(ref$X)),position_prior_init=ref$initial$pi,
  max_iter=4000,tol=1e-6,convergence='normalized')
 saveRDS(fit,file.path(out,paste0(tag,'_package.rds')))
 result <- data.frame(tag=tag,rho=abs(cor(ref$truth,drop(fit$gamma%*%seq(0,1,length.out=50)),method='spearman')),
  iterations=fit$iter,elbo_error=tail(fit$elbo_trace,1)-expected$elbo,max_R_error=max(abs(fit$gamma-input$R)))
 print(result);result
})
result <- do.call(rbind,results);write.csv(result,file.path(out,'package_audit.csv'),row.names=FALSE)
stopifnot(max(result$max_R_error)<1e-4,max(abs(result$elbo_error))<1e-3)
