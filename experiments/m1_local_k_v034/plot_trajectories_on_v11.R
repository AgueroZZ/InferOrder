# Plot every feature against V11 for the previously displayed failed case.
source('experiments/m2_local_k_v034/common.R')
library(ggplot2)
out <- 'experiments/m1_local_k_v034'
d <- readRDS(file.path(study,'data','rich_S4_r06.rds'))
cols <- which(d$truth==1)
cols <- cols[order(as.integer(sub('V','',colnames(d$X)[cols])))]
v11 <- match('V11',colnames(d$X));t <- seq(0,1,length.out=1001)
curves <- sapply(cols,function(j) {
 c <- d$coefficients[[j]];a <- outer(t,pi*d$frequencies)
 (as.numeric(sin(a)%*%(c$sine*d$attenuation)+cos(a)%*%(c$cosine*d$attenuation))-c$center)/c$scale
})
reference_index <- match(v11,cols)
stopifnot(!is.na(reference_index))
observed <- do.call(rbind,lapply(cols,function(j)data.frame(x=d$X[,v11],y=d$X[,j],feature=colnames(d$X)[j])))
paths <- do.call(rbind,lapply(seq_along(cols),function(j)data.frame(x=curves[,reference_index],y=curves[,j],feature=colnames(d$X)[cols[j]])))
observed$feature <- factor(observed$feature,levels=colnames(d$X)[cols])
paths$feature <- factor(paths$feature,levels=colnames(d$X)[cols])
p <- ggplot(observed,aes(x,y))+geom_point(color='gray55',alpha=.28,size=.6)+geom_path(data=paths,color='#0072B2',linewidth=.8)+facet_wrap(~feature,ncol=4)+theme_minimal(base_size=12)+labs(title='All 12 trajectories plotted against V11',subtitle='rich_S4_r06_A | M = 1 | 300 samples | SNR = 4',x='V11',y='Feature value',caption='Figure 1. Blue paths: exact noiseless feature pairs, connected in true latent-position order.\nGray points: paired noisy observations, including noise in V11. The V11 panel is the identity comparison.')
ggsave(file.path(out,'failed_case_trajectories_on_v11.png'),p,width=12,height=8,dpi=170,bg='white')
