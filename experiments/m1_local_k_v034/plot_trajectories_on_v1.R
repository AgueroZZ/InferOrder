# Plot every feature against V1 for the previously displayed failed case.
source('experiments/m2_local_k_v034/common.R')
library(ggplot2)
out <- 'experiments/m1_local_k_v034'
d <- readRDS(file.path(study,'data','rich_S4_r06.rds'))
cols <- which(d$truth==1)
cols <- cols[order(as.integer(sub('V','',colnames(d$X)[cols])))]
v1 <- match('V1',colnames(d$X));t <- seq(0,1,length.out=1001)
curves <- sapply(cols,function(j) {
 c <- d$coefficients[[j]];a <- outer(t,pi*d$frequencies)
 (as.numeric(sin(a)%*%(c$sine*d$attenuation)+cos(a)%*%(c$cosine*d$attenuation))-c$center)/c$scale
})
stopifnot(cols[1]==v1)
observed <- do.call(rbind,lapply(cols,function(j)data.frame(x=d$X[,v1],y=d$X[,j],feature=colnames(d$X)[j])))
paths <- do.call(rbind,lapply(seq_along(cols),function(j)data.frame(x=curves[,1],y=curves[,j],feature=colnames(d$X)[cols[j]])))
observed$feature <- factor(observed$feature,levels=colnames(d$X)[cols])
paths$feature <- factor(paths$feature,levels=colnames(d$X)[cols])
p <- ggplot(observed,aes(x,y))+geom_point(color='gray55',alpha=.28,size=.6)+geom_path(data=paths,color='#0072B2',linewidth=.8)+facet_wrap(~feature,ncol=4)+theme_minimal(base_size=12)+labs(title='All 12 trajectories plotted against V1',subtitle='rich_S4_r06_A | M = 1 | 300 samples | SNR = 4',x='V1',y='Feature value',caption='Figure 1. Blue paths: exact noiseless feature pairs, connected in true latent-position order.\nGray points: paired noisy observations, including noise in V1. The V1 panel is the identity comparison.')
ggsave(file.path(out,'failed_case_trajectories_on_v1.png'),p,width=12,height=8,dpi=170,bg='white')
