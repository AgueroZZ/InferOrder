source('experiments/estimate_intrinsic_m_smooth_v032/auto_m_v034/common.R')
out <- file.path(study_dir,'exploratory_m5_ordering_b/mechanism_half_noise')
Q <- readRDS(file.path(out,'native_prefix_fits.rds'))[['0']]$Q_K
expressions <- parse(file.path(out,'run_hybrids.R'))
for(e in expressions) if(is.call(e) && identical(e[[1]],as.name('<-')) &&
 identical(e[[2]],as.name('make_rw2_smoother')))eval(e)
t <- seq(0,1,length.out=300)^1.5
smooth <- make_rw2_smoother(5)
constant_error <- max(abs(smooth(t,rep(1,300))-1))
linear_error <- max(abs(smooth(t,2+3*t)-(2+3*t)))
edf <- environment(smooth)$cache$edf
stopifnot(constant_error<1e-7,linear_error<1e-7,abs(edf-5)<1e-6)
checks <- data.frame(constant_error=constant_error,linear_error=linear_error,edf=edf)
write.csv(checks,file.path(out,'hybrid_smoother_checks.csv'),row.names=FALSE)
files <- c('inspect.R','run_ablations.R','run_hybrids.R','render_report.R','verify_hybrid.R','frozen_cavi_source.R.txt')
write.csv(data.frame(file=files,sha256=vapply(file.path(out,files),function(f)digest::digest(file=f,algo='sha256'),character(1))),
 file.path(out,'source_hashes.csv'),row.names=FALSE)
print(checks)
