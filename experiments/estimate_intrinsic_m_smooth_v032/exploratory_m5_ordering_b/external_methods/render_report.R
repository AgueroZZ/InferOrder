# Assemble saved, unselected fits and plot initial versus final sample ranks.
source('experiments/estimate_intrinsic_m_smooth_v032/auto_m_v034/common.R')
suppressPackageStartupMessages(library(ggplot2))
base <- file.path(study_dir,'exploratory_m5_ordering_b')
out <- file.path(base,'external_methods')
rank01 <- function(x) (rank(x,ties.method='average')-1)/(length(x)-1)
rho <- function(x,y) abs(cor(x,y,method='spearman'))
rows <- list(); samples <- list()
add <- function(id,method,start,initial,final,truth,status,setting='',seed=NA_integer_) {
 stopifnot(length(initial)==300,length(final)==300,all(is.finite(final)))
 key <- paste(id,method,start,setting,seed,sep=':')
 rows[[key]] <<- data.frame(id=id,method=method,start=start,setting=setting,seed=seed,
  initial_rho=rho(initial,truth),final_rho=rho(final,truth),initial_final_rho=rho(initial,final),status=status)
 if(cor(initial,truth,method='spearman')<0) initial <- -initial
 if(cor(initial,final,method='spearman')<0) final <- -final
 samples[[key]] <<- data.frame(id=id,method=method,start=start,setting=setting,seed=seed,
  truth=truth,initial_rank=rank01(initial),final_rank=rank01(final),sample=seq_along(truth),status=status)
}
for(f in list.files(file.path(out,'princurve_results'),full.names=TRUE,pattern='rds$')) {
 d <- readRDS(f)
 add(d$id,'Principal curve',d$start_method,d$initial,d$final,d$truth,
  if(d$fit$converged)'Converged' else 'Iteration limit',d$setting)
}
for(f in list.files(file.path(out,'mpcurve_single_results'),full.names=TRUE,pattern='rds$')) {
 d <- readRDS(f)
 add(d$id,'MPCurve (B only)',d$start,d$initial,d$final,d$truth,
  if(d$fit$converged)'Converged' else 'Iteration limit')
}
for(f in list.files(file.path(out,'gpy_results'),full.names=TRUE,pattern='json$')) {
 d <- jsonlite::fromJSON(f)
 add(d$job$case,if(d$job$model=='GPLVM')'GPLVM' else paste0('Bayesian GPLVM (',d$job$inducing,' inducing)'),
  d$job$start,d$initial,d$final,d$truth,d$optimizer_status,seed=d$job$seed)
}
for(i in 1:3) {
 id <- c('B_noise0','B_noise05','B_noise1')[i]
 for(start in c('pca','isomap10','isomap15')) {
  f <- if(start=='pca') file.path(base,'pca_fits',sprintf('setting_%02d.rds',i)) else
   file.path(base,'sensitivity_fits',sprintf('setting_%02d.rds',i+if(start=='isomap10')3 else 6))
  d <- readRDS(f)
  initial <- if(start=='pca')d$initial_coordinates else d$isomap
  add(id,'MPCurve (five orderings)',start,initial,d$final,d$truth,
   if(d$status=='converged')'Converged' else d$status)
 }
}
summary <- do.call(rbind,rows); all_samples <- do.call(rbind,samples)
write.csv(summary,file.path(out,'all_comparisons.csv'),row.names=FALSE)
write.csv(all_samples,file.path(out,'all_rank_positions.csv'),row.names=FALSE)
# Native GPy PCA differs from centered, unscaled PCA; use explicit matched starts here.
selected <- all_samples[(all_samples$method %in% c('MPCurve (five orderings)','MPCurve (B only)')) |
 (all_samples$method=='Principal curve') |
 (all_samples$method %in% c('GPLVM','Bayesian GPLVM (50 inducing)') &
  all_samples$start!='pca' & all_samples$seed==20260929),]
selected <- selected[grepl('^B_',selected$id),]
selected$method_label <- selected$method
selected$method_label[selected$method=='Principal curve' & selected$setting=='package_defaults'] <- 'Principal curve: defaults'
selected$method_label[selected$method=='Principal curve' & selected$setting=='extended'] <- 'Principal curve: tighter stop'
selected$method_label[selected$method=='Bayesian GPLVM (50 inducing)'] <- 'Bayesian GPLVM (50 inducing)'
selected$method_label <- factor(selected$method_label,levels=c('MPCurve (five orderings)',
 'MPCurve (B only)','Principal curve: defaults','Principal curve: tighter stop','GPLVM','Bayesian GPLVM (50 inducing)'))
selected$start_label <- factor(ifelse(selected$start %in% c('pca','pca_matched'),'Centered, unscaled PCA',
 ifelse(selected$start=='isomap15','Isomap k = 15','Isomap k = 10')),
 levels=c('Centered, unscaled PCA','Isomap k = 15','Isomap k = 10'))
selected$noise_label <- factor(selected$id,levels=c('B_noise0','B_noise05','B_noise1'),
 labels=c('No noise','Half original noise','Original noise (SNR 4)'))
plot_pairs <- function(d,columns,title,subtitle,caption) {
 panels <- split(d,interaction(d[c('method_label',columns)],drop=TRUE))
 labels <- do.call(rbind,lapply(panels,function(x) {
  stopifnot(!anyDuplicated(x$sample))
  status <- if(x$status[1]=='Converged')'' else '\nOptimizer did not converge'
  cbind(x[1,c('method_label',columns),drop=FALSE],label=paste0(sprintf('|rho|: %.3f -> %.3f',
   rho(x$truth,x$initial_rank),rho(x$truth,x$final_rank)),status))
 }))
 long <- rbind(transform(d,stage='Initial',position=initial_rank),transform(d,stage='Final',position=final_rank))
 ggplot()+geom_segment(data=d,aes(x=truth,xend=truth,y=initial_rank,yend=final_rank),
  color='#777777',alpha=.16,linewidth=.2)+
  geom_point(data=long,aes(truth,position,color=stage,shape=stage),size=.8,alpha=.65)+
  geom_label(data=labels,aes(x=.02,y=1.18,label=label),hjust=0,vjust=1,size=3,linewidth=.12)+
  scale_color_manual(values=c(Initial='#0072B2',Final='#D55E00'),breaks=c('Initial','Final'))+
  scale_shape_manual(values=c(Initial=1,Final=16),breaks=c('Initial','Final'))+
  scale_y_continuous(breaks=c(0,.5,1))+coord_cartesian(xlim=c(0,1),ylim=c(0,1.2))+
  facet_grid(reformulate(columns,response='method_label'))+
  theme_minimal(base_size=11)+theme(legend.position='bottom',panel.grid.minor=element_blank(),
   plot.caption=element_text(hjust=0,size=9),strip.text.y=element_text(size=10))+
  labs(title=title,subtitle=subtitle,x='True latent position',y='Normalized sample rank',
   color=NULL,shape=NULL,caption=caption)
}
caption <- paste('Blue open circles: actual initial ranks; orange dots: final ranks. Gray lines connect the same sample.',
 'Absolute Spearman correlation removes global reversal. All single-order methods receive the same known 12 B features.',
 'MPCurve (five orderings) is the earlier full-data fit. Principal-curve tighter stop: maxit = 1000, thresh = 1e-6.',sep='\n')
p1 <- plot_pairs(selected[selected$start_label=='Centered, unscaled PCA' &
 selected$method_label!='Principal curve: defaults',],'noise_label',
 'Ordering B: can iteration repair the same PCA ordering?',
 'Three noise levels on the same signal and noise realization; GP models receive the exact R PCA coordinates.',caption)
p2 <- plot_pairs(selected[selected$id=='B_noise1',],'start_label',
 'Original challenging case: initialization changes the comparison',
 'Same upstream starts; principal-curve Isomap starts additionally require a smooth curve and an initial projection.',
 paste(caption,'Default principal curve: maxit = 10, thresh = 0.001. Its initial Isomap projection changes the raw ordering.',sep='\n'))
for(ext in c('png','pdf')) {
 ggsave(file.path(out,paste0('matched_pca_three_noise.',ext)),p1,width=13,height=14,dpi=170,bg='white')
 ggsave(file.path(out,paste0('original_noise_matched_starts.',ext)),p2,width=13,height=16.5,dpi=170,bg='white')
}
paths <- read.csv(file.path(out,'princurve_iteration_paths.csv'))
paths$panel <- ifelse(paths$start=='pca','Default PCA initialization','Isomap k = 15 initialization')
paths$noise <- factor(paths$id,levels=c('B_noise0','B_noise05','B_noise1'),labels=c('No noise','Half original noise','Original noise'))
default <- read.csv(file.path(out,'princurve_summary.csv'))
default <- default[default$setting=='package_defaults' & grepl('^B_',default$id) &
 (default$start_method=='pca' | (default$id=='B_noise1' & default$start_method=='isomap15')),]
default$panel <- ifelse(default$start_method=='pca','Default PCA initialization','Isomap k = 15 initialization')
p3 <- ggplot(paths,aes(iteration,rho,color=noise))+geom_line(linewidth=.8)+geom_point(size=2)+
 geom_point(data=default,aes(x=iterations,y=final_rho),inherit.aes=FALSE,shape=4,size=4,stroke=1.2)+
 facet_wrap(~panel,nrow=1)+scale_color_manual(values=c('#0072B2','#009E73','#D55E00'))+
 scale_y_continuous(limits=c(0,1),breaks=seq(0,1,.25))+theme_minimal(base_size=12)+
 theme(legend.position='bottom',plot.caption=element_text(hjust=0,size=10))+
 labs(title='Principal curve: recovery can appear late, or deteriorate with further updates',
  subtitle='Saved iteration prefixes from identical starts; black crosses mark the package-default endpoints.',
  x='Completed principal-curve iterations',y='Absolute Spearman correlation with truth',color=NULL,
  caption=paste('Lines join sampled iteration counts; tighter stopping uses maxit = 1000 and thresh = 1e-6.',
   'The smoother stays at its default df = 5. Truth is used only for evaluation, never to select a stopping point.',sep='\n'))
for(ext in c('png','pdf')) ggsave(file.path(out,paste0('princurve_iteration_recovery.',ext)),p3,width=12,height=5.5,dpi=180,bg='white')
# Explicit checks of matched upstream coordinates and complete output counts.
stopifnot(nrow(summary)==69L,length(list.files(file.path(out,'gpy_results'),pattern='json$'))==25L)
for(id in c('B_noise0','B_noise05','B_noise1')) {
 ref <- summary$initial_rho[summary$id==id & summary$method=='MPCurve (five orderings)' & summary$start=='pca']
 matched <- summary$initial_rho[summary$id==id & summary$start=='pca_matched']
 stopifnot(all(abs(matched-ref)<1e-12))
}
message('Rendered three figures; verified 69 fit records and exact matched-PCA rank correlations.')
