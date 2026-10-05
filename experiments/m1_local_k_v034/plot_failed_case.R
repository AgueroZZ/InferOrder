# Visualize the previously displayed largest regression, without changing fits.
source('experiments/m2_local_k_v034/common.R')
library(ggplot2)
out <- 'experiments/m1_local_k_v034'
d <- readRDS(file.path(study,'data','rich_S4_r06.rds'))
r <- readRDS(file.path(out,'results','rich_S4_r06_A.rds'))
cols <- which(d$truth==1)
cols <- cols[order(as.integer(sub('V','',colnames(d$X)[cols])))]
t <- seq(0,1,length.out=1001)
curves <- sapply(cols,function(j) {
 c <- d$coefficients[[j]];a <- outer(t,pi*d$frequencies)
 (as.numeric(sin(a)%*%(c$sine*d$attenuation)+cos(a)%*%(c$cosine*d$attenuation))-c$center)/c$scale
})
colnames(curves) <- colnames(d$X)[cols]
observed <- do.call(rbind,lapply(cols,function(j)data.frame(t=d$latent[,1],value=d$X[,j],feature=colnames(d$X)[j])))
lines <- do.call(rbind,lapply(seq_along(cols),function(j)data.frame(t=t,value=curves[,j],feature=colnames(curves)[j])))
observed$feature <- factor(observed$feature,levels=colnames(curves))
lines$feature <- factor(lines$feature,levels=colnames(curves))
theme_set(theme_minimal(base_size=12))
p <- ggplot(observed,aes(t,value))+geom_point(color='gray55',alpha=.28,size=.55)+geom_line(data=lines,color='#0072B2',linewidth=.8)+facet_wrap(~feature,ncol=4)+labs(title='A failed M=1 case: all 12 generating trajectories',subtitle='rich_S4_r06_A | 300 samples | SNR = 4 | default k=15: recovery 0.806; selected k=10: 0.513',x='True latent position',y='Feature value',caption='Figure 1. Blue curves: exact noiseless trajectories; gray points: observations (noise SD = 0.5).\nEach signal has unit sample variance. The case is the largest regression in the preceding M=1 comparison.')
ggsave(file.path(out,'failed_case_trajectories.png'),p,width=12,height=8,dpi=170,bg='white')
# Choose the first two feature IDs, independently of their projected geometry.
pair <- cols[1:2];stopifnot(identical(colnames(d$X)[pair],c('V1','V2')))
path <- data.frame(x=curves[,1],y=curves[,2],t=t)
obs <- data.frame(x=d$X[,pair[1]],y=d$X[,pair[2]],t=d$latent[,1])
p2 <- ggplot()+geom_path(data=transform(path,panel='Noiseless trajectories'),aes(x,y,color=t),linewidth=1.1)+geom_point(data=transform(obs,panel='Observed features'),aes(x,y,color=t),size=1.5,alpha=.8)+facet_wrap(~panel,nrow=1)+scale_color_viridis_c(name='True latent\nposition',limits=c(0,1))+coord_equal()+labs(title='Two trajectories plotted against each other: V1 versus V2',subtitle='The same samples and generating functions as Figure 1; colors track progression along the true ordering.',x='V1',y='V2',caption='Figure 2. Left: exact noiseless path; right: all 300 noisy observations. V1 and V2 were chosen by feature ID.\nThis two-feature projection illustrates the geometry; crossings here do not establish ambiguity in the full 12-feature space.')
ggsave(file.path(out,'failed_case_feature_pair.png'),p2,width=11,height=5.8,dpi=170,bg='white')
stopifnot(identical(r$input_hash,digest(d$X[,d$truth==1,drop=FALSE],algo='sha256')))
cat('Verified case input hash; plotted 12 trajectories and prespecified feature pair V1/V2.\n')
