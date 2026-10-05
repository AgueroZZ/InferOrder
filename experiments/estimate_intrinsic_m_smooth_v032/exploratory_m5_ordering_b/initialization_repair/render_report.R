# Summarize prespecified perturbations, fixed-sweep continuation, and mechanisms.
source('experiments/estimate_intrinsic_m_smooth_v032/exploratory_m5_ordering_b/initialization_repair/common.R')
suppressPackageStartupMessages({library(ggplot2);library(patchwork)})
paths <- file.path(repair_dir,'results',paste0(scenarios$id,'.rds'))
stopifnot(all(file.exists(paths)))
results <- lapply(paths,readRDS); names(results) <- scenarios$id
stopifnot(all(vapply(results,function(x)identical(x$input_hash,dataset$input_hash),logical(1))),
 all(vapply(results,function(x)isTRUE(x$converged),logical(1))))
label_group <- function(x) {
 switch(x$scenario$type,
  oracle='Reference initializations',oracle_reversed='Reference initializations',isomap10='Reference initializations',
  isomap15='Isomap, k = 15',pca='PCA, first component',
  jitter=sprintf('Gaussian jitter, SD = %.1f',x$scenario$magnitude),
  swap_middle_blocks='Swap adjacent middle quarters',reverse_middle_half='Reverse the middle half',
  circular_shift='Circular shift by one quarter',random='Random ordering')
}
summary <- do.call(rbind,lapply(results,function(x) data.frame(
 id=x$scenario$id,group=label_group(x),raw_rho=x$raw_metrics['rho'],
 warm_rho=x$warm_metrics['rho'],final_rho=x$final_metrics['rho'],
 raw_pair_error=x$raw_metrics['pair_error'],warm_pair_error=x$warm_metrics['pair_error'],
 final_pair_error=x$final_metrics['pair_error'],raw_rank_error=x$raw_metrics['mean_rank_error'],
 final_rank_error=x$final_metrics['mean_rank_error'],objective=x$objective,
 iterations=x$iterations,converged=x$converged,ARI=x$ARI,warnings=length(x$warnings))))
rownames(summary)<-NULL
write.csv(summary,file.path(repair_dir,'repair_summary.csv'),row.names=FALSE)
# Baseline instrumentation must reproduce the previously completed runs.
previous <- read.csv(file.path(study_dir,'exploratory_m5_ordering_b','pca_vs_isomap_summary.csv'))
for(id in c('isomap10','isomap15','pca')) {
 method <- c(isomap10='Isomap (k=10)',isomap15='Isomap (k=15)',pca='PCA (PC1)')[[id]]
 row <- previous[previous$method==method & previous$noise_scale==1,]
 stopifnot(abs(results[[id]]$objective-row$objective)<1e-5,
  abs(results[[id]]$final_metrics['rho']-row$final_rho)<1e-8)
}
stopifnot(abs(results$oracle$objective-results$oracle_reversed$objective)<1e-5)
stages <- do.call(rbind,lapply(results,function(x) data.frame(id=x$scenario$id,group=label_group(x),
 stage=factor(c('Raw ordering','After 2 CAVI updates','After annealing','Final fit'),
 levels=c('Raw ordering','After 2 CAVI updates','After annealing','Final fit')),
 rho=c(x$raw_metrics['rho'],x$warm_metrics['rho'],x$trajectory$rho[25],x$final_metrics['rho']))))
p <- ggplot(stages,aes(stage,rho,group=id))+geom_line(color='#0072B2',alpha=.75)+
 geom_point(color='#0072B2',size=1.6)+facet_wrap(~group,ncol=3)+
 scale_y_continuous(limits=c(0,1),breaks=c(0,.25,.5,.75,1))+
 theme_minimal(base_size=11)+theme(axis.text.x=element_text(angle=25,hjust=1),panel.grid.minor=element_blank())+
 labs(title='Which initialization errors can MPCurve repair on this dataset?',
 subtitle='Same M=5, SNR=4 dataset; only B initialization changes. Lines are separate starts, not independent datasets.',
 x=NULL,y='Absolute Spearman correlation with true B position',
 caption='Two CAVI initialization updates are shown separately from joint structural fitting.\nOracle positions and artificial perturbations are diagnostic controls, not usable estimators.')
ggsave(file.path(repair_dir,'repair_by_error_type.png'),p,width=14,height=12,dpi=170,bg='white')
selected <- c('isomap10','isomap15','pca','jitter_030_r1','jitter_060_r1','swap_middle_blocks','random_r1')
traces <- do.call(rbind,lapply(selected,function(id) {
 x<-results[[id]]; data.frame(x$trajectory,id=id)
}))
trace_labels <- c(isomap10='Isomap (k=10)',isomap15='Isomap (k=15)',pca='PCA',
 jitter_030_r1='Jitter SD 0.3, seed 1',jitter_060_r1='Jitter SD 0.6, seed 1',
 swap_middle_blocks='Swap middle quarters',random_r1='Random, seed 1')
trace_colors <- setNames(c('#0072B2','#D55E00','#009E73','#CC79A7','#E69F00','#56B4E9','#444444'),selected)
p2 <- ggplot(traces,aes(sweep,rho,color=id))+geom_line(linewidth=.65)+
 geom_vline(xintercept=25,linetype=2,color='gray50')+
 scale_color_manual(values=trace_colors,labels=trace_labels)+
 theme_minimal(base_size=12)+theme(legend.position='bottom')+guides(color=guide_legend(nrow=2))+
 labs(title='Ordering changes during joint structural iterations',
 subtitle='Sweep 1 begins after the two subset CAVI initialization updates; dashed line marks the end of annealing.',
 x='Structural sweep',y='Absolute Spearman correlation',color=NULL)
ggsave(file.path(repair_dir,'ordering_iteration_traces.png'),p2,width=12,height=6,dpi=170,bg='white')
continuations <- lapply(c('isomap10','isomap15','pca'),function(id) {
 x <- readRDS(file.path(repair_dir,'results',paste0(id,'_forced_continuation.rds')))
 data.frame(x$trajectory,id=id)
})
continuation <- do.call(rbind,continuations)
continuation_summary <- do.call(rbind,lapply(continuations,function(x)data.frame(id=x$id[1],
 before_rho=x$rho[1],after_rho=tail(x$rho,1),objective_gain=tail(x$objective,1)-x$objective[1],
 mean_rank_change=tail(x$mean_rank_change,1),last_normalized_increment=tail(diff(x$objective),1)/(300*60))))
write.csv(continuation_summary,file.path(repair_dir,'continuation_summary.csv'),row.names=FALSE)
p3 <- ggplot(continuation,aes(extra_sweeps,rho,color=id))+geom_line(linewidth=.7)+
 scale_color_manual(values=trace_colors,labels=trace_labels)+
 theme_minimal(base_size=12)+theme(legend.position='bottom')+labs(title='Disabling early stopping does not unfold these solutions',
 subtitle='500 additional ordinary T=1 sweeps after the original convergence flag.',
 x='Additional sweeps after convergence',y='Absolute Spearman correlation',color=NULL)
ggsave(file.path(repair_dir,'forced_continuation.png'),p3,width=10,height=5,dpi=170,bg='white')
# Anchor trajectories explain how the fitted curve can accommodate a bad order.
anchor_rows <- do.call(rbind,lapply(c('isomap10','isomap15','pca'),function(id) {
 x <- results[[id]]
 pos <- as.numeric(x$gamma %*% grid)
 mu <- x$conditional_mean[anchor,]
 curve_grid <- grid
 if(cor(truth,pos,method='spearman')<0) {pos<-1-pos;curve_grid<-1-grid}
 label <- sprintf('%s: rho=%.3f, lambda=%.1f',trace_labels[[id]],x$final_metrics['rho'],x$lambda_mat[anchor,x$slot])
 data.frame(id=id,label=label,truth=truth,position=pos,observation=dataset$X[,anchor])
}))
curves <- do.call(rbind,lapply(c('isomap10','isomap15','pca'),function(id) {
 x <- results[[id]]; pos <- as.numeric(x$gamma %*% grid); curve_grid<-grid
 if(cor(truth,pos,method='spearman')<0) curve_grid<-1-grid
 data.frame(label=unique(anchor_rows$label[anchor_rows$id==id]),position=curve_grid,
 mean=x$conditional_mean[anchor,])
}))
p4 <- ggplot(anchor_rows,aes(position,observation,color=truth))+geom_point(size=.9,alpha=.65)+
 geom_line(data=curves,aes(position,mean),inherit.aes=FALSE,color='black',linewidth=.6)+
 facet_wrap(~label,nrow=1)+scale_color_viridis_c(option='D')+theme_minimal(base_size=11)+
 labs(title='The monotone anchor is fitted as a more irregular curve under bad orderings',
 subtitle='V20 observations versus inferred position; black line is its conditional posterior mean across the 50 grid points.',
 x='Inferred B position (orientation aligned)',y='Observed anchor value',color='True B position',
 caption='Lambda is the fitted RW2 precision: smaller values allow more trajectory curvature.\nThe anchor is not constrained to be monotone by the fitted model.')
ggsave(file.path(repair_dir,'anchor_trajectory_adaptation.png'),p4,width=14,height=5,dpi=170,bg='white')
print(summary,row.names=FALSE)
print(continuation_summary,row.names=FALSE)
