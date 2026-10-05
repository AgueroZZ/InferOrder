source('experiments/estimate_intrinsic_m_smooth_v032/auto_m_v034/common.R')
suppressPackageStartupMessages(library(ggplot2))
base <- file.path(study_dir,'exploratory_m5_ordering_b');out <- file.path(base,'mechanism_half_noise')
fits <- readRDS(file.path(out,'ablation_fits.rds'))$results
hybrids <- readRDS(file.path(out,'hybrid_fits.rds'))
pos <- read.csv(file.path(base,'external_methods/inputs/B_noise05_positions.csv'))
rank01 <- function(x)(rank(x,ties.method='average')-1)/(length(x)-1)
ids <- c('native','fixed_lambda_10000','pcurve_df5','pcurve_df10','rw2_smoother_df5','spline_rank_spacing')
labels <- c('MPCurve: native','MPCurve: fixed stronger penalty','P-curve: cubic spline, df = 5',
 'P-curve: cubic spline, df = 10','P-curve: RW2 smoother, df = 5','P-curve: rank-spaced smoother')
frames <- lapply(seq_along(ids),function(i) {
 fit <- if(ids[i] %in% names(fits))fits[[ids[i]]] else hybrids[[ids[i]]]
 initial <- pos$pca; if(cor(initial,pos$truth,method='spearman')<0) initial <- -initial
 final <- fit$final; if(cor(initial,final,method='spearman')<0)final <- -final
 data.frame(panel=factor(labels[i],levels=labels),truth=pos$truth,initial=rank01(initial),final=rank01(final),
  label=sprintf('|rho|: %.3f -> %.3f',abs(cor(pos$truth,initial,method='spearman')),abs(cor(pos$truth,final,method='spearman'))))
})
d <- do.call(rbind,frames)
long <- rbind(transform(d,stage='Initial PCA',position=initial),transform(d,stage='Final',position=final))
p <- ggplot()+geom_segment(data=d,aes(x=truth,xend=truth,y=initial,yend=final),color='#888888',alpha=.18,linewidth=.2)+
 geom_point(data=long,aes(truth,position,color=stage,shape=stage),alpha=.7,size=1)+
 geom_label(data=unique(d[c('panel','label')]),aes(x=.02,y=1.13,label=label),hjust=0,vjust=1,size=3.2,linewidth=.15)+
 facet_wrap(~panel,ncol=2)+scale_color_manual(values=c('Initial PCA'='#0072B2',Final='#D55E00'),breaks=c('Initial PCA','Final'))+
 scale_shape_manual(values=c('Initial PCA'=1,Final=16),breaks=c('Initial PCA','Final'))+
 coord_cartesian(xlim=c(0,1),ylim=c(0,1.17))+scale_y_continuous(breaks=c(0,.5,1))+
 theme_minimal(base_size=12)+theme(legend.position='bottom',plot.caption=element_text(hjust=0,size=9))+
 labs(title='Half-noise B: the penalty family does not determine recovery',
  subtitle='Controlled diagnostic variants from the same PCA ordering; all six displayed runs converged.',
  x='True latent position',y='Normalized sample rank',color=NULL,shape=NULL,
  caption=paste('The RW2 hybrid retains P-curve smoothing/projection/arc-length updates, using 50 regular RW2 basis nodes and df = 5.',
   'Rank-spacing changes only the smoother input coordinates to ranks; the projection step is unchanged.',
   'MPCurve fixed penalty: lambda = 10000, final conditional smoother df = 5.05-5.91. Blue/orange show initial/final ranks.',
   'These are mechanism diagnostics on one fixed dataset, not new default methods or truth-selected benchmark settings.',sep='\n'))
for(ext in c('png','pdf'))ggsave(file.path(out,paste0('mechanism_comparison.',ext)),p,width=11,height=11,dpi=180,bg='white')
feature <- read.csv(file.path(out,'native_feature_diagnostics.csv'))
traces <- do.call(rbind,lapply(split(feature,feature$budget),function(x)
 data.frame(iteration=x$iteration[1],rho=x$rho[1],edf_min=min(x$conditional_edf),edf_max=max(x$conditional_edf))))
traces <- unique(traces[order(traces$iteration),])
write.csv(traces,file.path(out,'native_iteration_summary.csv'),row.names=FALSE)
# Verify the diagnostic hybrids preserve the principal-curve endpoint when expected.
stopifnot(abs(fits$pcurve_df5$final[1]-readRDS(file.path(base,'external_methods/princurve_results/B_noise05_pca_extended.rds'))$final[1])<1e-10)
native <- readRDS(file.path(base,'single_order_benchmark/B_noise05_mpcurve_extended_stopping.rds'))
stopifnot(max(abs(fits$native$final-native$final))<1e-12)
noise <- load_fixed_dataset(manifest[manifest$id=='main_M5_S4_r001',])$noise_variance
stopifnot(all(abs(noise*.5^2-.0625)<1e-12))
summary <- data.frame(mean_max_position_probability=mean(apply(fits$native$fit$gamma,1,max)),
 occupied_MAP_bins=length(unique(max.col(fits$native$fit$gamma))),
 elbo=tail(fits$native$fit$elbo_trace,1),nominal_half_noise_variance=noise*.5^2)
write.csv(summary,file.path(out,'additional_diagnostics.csv'),row.names=FALSE)
print(summary,row.names=FALSE)
