# Compare PCA initialization with the matched Isomap controls at three noise levels.
source('experiments/estimate_intrinsic_m_smooth_v032/auto_m_v034/common.R')
suppressPackageStartupMessages(library(ggplot2))
out <- file.path(study_dir, 'exploratory_m5_ordering_b')
source(file.path(out, 'plot_ordering_comparison.R'))
pca <- lapply(file.path(out, 'pca_fits', sprintf('setting_%02d.rds', 1:3)), readRDS)
stopifnot(all(vapply(pca, function(x) x$status == 'converged', logical(1))))
make_samples <- function(x, method) data.frame(sample = seq_along(x$truth),
 truth = x$truth, isomap = if (method == 'PCA (PC1)') x$initial_coordinates else x$isomap,
 final = x$final, noise_scale = x$setting$noise_scale, method = method)
samples <- do.call(rbind, lapply(pca, make_samples, method = 'PCA (PC1)'))
comparison <- plot_ordering_comparison(samples, 'noise_scale',
 title = 'Ordering B: PCA initialization and final MPCurve ordering',
 subtitle = 'PC1 of the same 12 B features, centered without variance scaling; only B initialization changes.',
 initial_label = 'Initial PCA')
caption <- paste('Blue open circles: initial ranks; orange dots: final posterior-mean ranks. Gray segments connect the same sample.',
 'Normalized average ranks remove coordinate stretching; global reversal is aligned away.',
 'Noise scale 1 is the original data; 0.5 halves the same B noise realization; 0 uses the noiseless B signal.',
 'All full-data fits use M=5, with original group initialization and unchanged annealing/convergence controls.', sep = '\n')
p <- comparison$plot + facet_wrap(~noise_scale, nrow = 1, labeller = label_both) +
 theme(plot.caption = element_text(hjust = 0, size = 9), plot.caption.position = 'plot') + labs(caption = caption)
for (ext in c('png','pdf')) ggsave(file.path(out,paste0('pca_initial_final_comparison.',ext)),p,width=13,height=5.5,dpi=180,bg='white')
write.csv(comparison$aligned, file.path(out,'pca_initial_final_positions.csv'),row.names=FALSE)
all_results <- pca
methods <- rep('PCA (PC1)',3)
for (k in c(10,15)) {
 indices <- if (k==10) 4:6 else 7:9
 iso <- lapply(file.path(out,'sensitivity_fits',sprintf('setting_%02d.rds',indices)),readRDS)
 for (i in 1:3) stopifnot(identical(iso[[i]]$modified_input_hash,pca[[i]]$modified_input_hash))
 all_results <- c(all_results,iso)
 methods <- c(methods,rep(paste0('Isomap (k=',k,')'),3))
}
combined <- do.call(rbind,Map(make_samples,all_results,methods))
combined$method <- factor(combined$method,levels=unique(methods))
combined_plot <- plot_ordering_comparison(combined,c('method','noise_scale'),
 title='Ordering B: PCA and Isomap under the same three noise settings',
 subtitle='Only B initialization differs between methods; each panel overlays the initial and final sample ordering.',
 initial_label='Initial ordering')
p2 <- combined_plot$plot + facet_grid(method~noise_scale,labeller=label_both) +
 theme(plot.caption = element_text(hjust=0,size=9),plot.caption.position='plot') + labs(caption=caption)
for(ext in c('png','pdf')) ggsave(file.path(out,paste0('pca_vs_isomap_comparison.',ext)),p2,width=13,height=11,dpi=180,bg='white')
summary <- do.call(rbind,Map(function(x,method) {
 initial <- if(method=='PCA (PC1)') x$initial_coordinates else x$isomap
 data.frame(method=method,noise_scale=x$setting$noise_scale,
  initial_rho=abs(cor(x$truth,initial,method='spearman')),
  final_rho=abs(cor(x$truth,x$final,method='spearman')),
  initial_final_rho=abs(cor(initial,x$final,method='spearman')),
  objective=tail(x$objective,1),iterations=x$iterations,ARI=x$ARI,
  status=x$status,warning_count=length(x$warnings))
},all_results,methods))
write.csv(summary,file.path(out,'pca_vs_isomap_summary.csv'),row.names=FALSE)
print(summary,row.names=FALSE)
