# Shared, fixed design for the replicated single-ordering comparison.
options(stringsAsFactors = FALSE)
.libPaths(c('experiments/estimate_intrinsic_m_smooth_v032/auto_m_v034/library',
            'experiments/estimate_intrinsic_m_smooth_v032/library', .libPaths()))
stopifnot(as.character(packageVersion('MPCurver')) == '0.3.4',
          as.character(packageVersion('princurve')) == '2.1.6')
study <- 'experiments/m1_bspline_comparison'
design <- list(N=200L, P=12L, replications=30L, spline_df=8L, degree=3L,
               noise_sd=c(low=0.25, high=1), seed_base=202609290L,
               K=50L, mp_max_iter=2000L, mp_tol=1e-6,
               pc_max_iter=1000L, pc_tol=1e-6, gp_block_budget=2000L,
               gp_max_blocks=3L, inducing=50L,
               package_commit='15f2b0bbe5dfa61cd46da5160b2bc251e75a0475')
recovery <- function(truth, position, method='spearman') {
  if (any(!is.finite(position)) || length(unique(position)) < 2L) return(NA_real_)
  abs(cor(truth, position, method=method))
}
