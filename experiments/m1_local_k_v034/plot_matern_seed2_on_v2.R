# Pairwise view of the second Matérn realization using V2 as the reference.
source('experiments/m2_local_k_v034/common.R')
library(ggplot2)
out <- 'experiments/m1_local_k_v034/matern_seeds'
x <- readRDS(file.path(out,'results.rds'))$results[['20260932']]
d <- readRDS(file.path(study,'data','rich_S4_r06.rds'))
features <- colnames(d$X)[d$truth==1]
reference <- match('V2',features)
ord <- order(as.integer(sub('V','',features)))
latent_order <- order(x$truth)
stopifnot(!is.na(reference),ncol(x$X)==length(features))
observed <- do.call(rbind,lapply(ord,function(j)data.frame(x=x$X[,reference],y=x$X[,j],feature=features[j])))
paths <- do.call(rbind,lapply(ord,function(j)data.frame(x=x$signal[latent_order,reference],y=x$signal[latent_order,j],feature=features[j])))
observed$feature <- factor(observed$feature,levels=features[ord])
paths$feature <- factor(paths$feature,levels=features[ord])
p <- ggplot(observed,aes(x,y))+geom_point(color='gray55',alpha=.28,size=.6)+geom_path(data=paths,color='#0072B2',linewidth=.8)+facet_wrap(~feature,ncol=4)+theme_minimal(base_size=12)+labs(title='Second Matern seed: all 12 trajectories plotted against V2',subtitle='Seed 20260932 | Matern 5/2 | length-scale = 0.2 | SNR = 4',x='V2',y='Feature value',caption='Figure 1. Blue paths: noiseless feature pairs, connected in true-position order at the 300 sampled positions.\nGray points: paired noisy observations, including noise in V2. The V2 panel is the identity comparison.')
ggsave(file.path(out,'seed2_trajectories_on_v2.png'),p,width=12,height=8,dpi=170,bg='white')
