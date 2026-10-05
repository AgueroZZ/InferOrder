source('experiments/estimate_intrinsic_m_smooth_v032/auto_m_v034/common.R')
out <- file.path(study_dir,'exploratory_m5_ordering_b/collapsed_comparison')
inputs <- jsonlite::fromJSON(file.path(out,'endpoint_inputs.json'),simplifyVector=FALSE)
checks <- lapply(inputs,function(s) {
 ref <- jsonlite::fromJSON(file.path(out,paste0(s$case,'_reference.json')))
 R <- do.call(rbind,lapply(s$R,unlist));attr(R,'cavi_skip_raw_preinit') <- TRUE
 fit <- MPCurver:::cavi(ref$X,K=50,responsibilities_init=R,
  sigma2_init=unlist(s$sigma2),lambda_init=unlist(s$lambda_),
  position_prior_init=unlist(s$pi),max_iter=0)
 data.frame(case=s$case,method=s$method,python_elbo=s$expected,
  package_elbo=tail(fit$elbo_trace,1),error=tail(fit$elbo_trace,1)-s$expected)
})
result <- do.call(rbind,checks);print(result)
write.csv(result,file.path(out,'endpoint_audit.csv'),row.names=FALSE)
stopifnot(max(abs(result$error))<1e-4)
