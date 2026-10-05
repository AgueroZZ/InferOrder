# Detailed view of the second prespecified Matérn seed, without refitting.
source('experiments/m2_local_k_v034/common.R')
library(ggplot2)
out<-'experiments/m1_local_k_v034/matern_seeds'
x<-readRDS(file.path(out,'results.rds'))$results[['20260932']]
f<-x$fit;t<-x$truth
rho<-c(abs(cor(t,x$before,method='spearman')),abs(cor(t,x$after,method='spearman')))
elbo<-c(head(f$elbo_trace,1),tail(f$elbo_trace,1))
labels<-c(sprintf('PCA initialization\nrho = %.3f | ELBO = %.2f',rho[1],elbo[1]),sprintf('MPCurve: %d iterations\nrho = %.3f | ELBO = %.2f',f$iter,rho[2],elbo[2]))
theme_set(theme_minimal(base_size=12))
p<-ggplot(data.frame(truth=rep(t,2),estimate=c(x$before,x$after),stage=factor(rep(labels,each=300),levels=labels)),aes(truth,estimate))+geom_abline(slope=1,intercept=0,color='gray70',linetype=2)+geom_point(size=1.2,alpha=.65,color='gray25')+facet_wrap(~stage,nrow=1)+coord_equal(xlim=c(0,1),ylim=c(0,1))+labs(title='Second Matern seed: the folded ordering persists',subtitle='Seed 20260932 | Matern 5/2 | length-scale 0.2 | SNR 4 | K50 | fixed uniform prior',x='True latent position',y='Estimated latent position',caption='Figure 1. Left: scaled PC scores. Right: posterior mean positions. rho is absolute Spearman correlation.\nInitial ELBO follows quantile/variational initialization; both panels share orientation.')
ggsave(file.path(out,'seed2_ordering.png'),p,width=11,height=6,dpi=170,bg='white')
d<-readRDS(file.path(study,'data','rich_S4_r06.rds'));names<-colnames(d$X)[d$truth==1];ord<-order(as.integer(sub('V','',names)))
obs<-do.call(rbind,lapply(ord,function(j)data.frame(truth=t,value=x$X[,j],feature=names[j])))
lines<-do.call(rbind,lapply(ord,function(j)data.frame(truth=t,value=x$signal[,j],feature=names[j])))
for(n in c('obs','lines')){a<-get(n);a$feature<-factor(a$feature,levels=names[ord]);assign(n,a)}
p2<-ggplot(obs,aes(truth,value))+geom_point(color='gray55',alpha=.28,size=.55)+geom_line(data=lines,color='#0072B2',linewidth=.8)+facet_wrap(~feature,ncol=4)+labs(title='Second Matern seed: all 12 true trajectories',subtitle='Seed 20260932 | blue: noiseless signal | gray: observed data',x='True latent position',y='Feature value',caption='Figure 2. Noiseless values are joined in true-position order at the 300 sampled positions.\nEach signal has unit sample variance; noise SD is 0.5. No monotone anchor was imposed.')
ggsave(file.path(out,'seed2_trajectories.png'),p2,width=12,height=8,dpi=170,bg='white')
obs$estimate<-rep(x$after,12)
p3<-ggplot(obs,aes(estimate,value,color=truth))+geom_point(size=.8,alpha=.65)+facet_wrap(~feature,ncol=4)+scale_color_viridis_c(limits=c(0,1),name='True latent\nposition')+scale_x_continuous(limits=c(0,1))+labs(title='Second Matern seed: observations along the inferred positions',x='MPCurve posterior mean position',y='Observed feature value',caption='Figure 3. Same observations as Figure 2, plotted against the converged positions; colors show true latent position.')
ggsave(file.path(out,'seed2_columns.png'),p3,width=12,height=8,dpi=170,bg='white')
cat('rho:',rho,'ELBO:',elbo,'iterations:',f$iter,'converged:',f$converged,'before-after rho:',cor(x$before,x$after,method='spearman'),'\n')
