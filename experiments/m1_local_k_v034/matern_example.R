# Exploratory single Matern realization paired with the Fourier case's positions/noise.
source('experiments/m2_local_k_v034/common.R')
library(ggplot2)
out <- 'experiments/m1_local_k_v034/matern_example'
dir.create(out,recursive=TRUE,showWarnings=FALSE)
d <- readRDS(file.path(study,'data','rich_S4_r06.rds'))
cols<-which(d$truth==1);truth<-d$latent[,1];ell<-.2;seed<-20260930L
# Sample on observed and plotting positions jointly; exact Matern-5/2 kernel.
grid<-seq(0,1,length.out=601);tt<-c(truth,grid)
z<-sqrt(5)*abs(outer(tt,tt,'-'))/ell
covariance<-(1+z+z^2/3)*exp(-z)
set.seed(seed);draw<-t(chol(covariance+diag(1e-10,length(tt))))%*%matrix(rnorm(length(tt)*12),length(tt),12)
center<-colMeans(draw[1:300,]);scale<-apply(draw[1:300,],2,sd)
draw<-sweep(sweep(draw,2,center,'-'),2,scale,'/')
signal<-draw[1:300,];smooth<-draw[-(1:300),];X<-signal+d$unit_noise[,cols]/2
colnames(X)<-colnames(signal)<-colnames(smooth)<-colnames(d$X)[cols]
raw<-as.numeric(MPCurver:::isomap_ordering(X,k=15L)$t)
init<-MPCurver:::.cavi_exact_k_init_from_ordering(X=X,ordering_vec=raw,K=50L,discretization='quantile',strict_K=TRUE)
gamma0<-MPCurver:::.cavi_hard_gamma_from_cluster_rank(init$cluster_rank,300,50L)
f<-MPCurver:::cavi(X=X,K=50L,responsibilities_init=gamma0,position_prior_init=rep(1/50,50),position_prior='fixed',sigma2_init=init$sigma2,lambda_init=rep(1,12),rw_q=2L,ridge=0,discretization='quantile',max_iter=2000L,tol=1e-6,convergence='relative',verbose=FALSE)
while(!f$converged && f$iter<10000L)f<-MPCurver:::do_cavi(f,iter=2000L,tol=1e-6)
stopifnot(f$converged,max(abs(f$params$pi-.02))<1e-12)
before<-(raw-min(raw))/diff(range(raw));after<-as.numeric(f$gamma%*%seq(0,1,length.out=50))
if(cor(before,truth,method='spearman')<0){before<-1-before;after<-1-after}
rho<-c(abs(cor(before,truth,method='spearman')),abs(cor(after,truth,method='spearman')))
elbo<-c(head(f$elbo_trace,1),tail(f$elbo_trace,1))
labels<-c(sprintf('Isomap k=15 (iteration 0)\nrho = %.3f | ELBO = %.2f',rho[1],elbo[1]),sprintf('MPCurve (iteration %d)\nrho = %.3f | ELBO = %.2f',f$iter,rho[2],elbo[2]))
theme_set(theme_minimal(base_size=12))
p<-ggplot(data.frame(truth=rep(truth,2),estimate=c(before,after),stage=factor(rep(labels,each=300),levels=labels)),aes(truth,estimate))+geom_abline(slope=1,intercept=0,color='gray70',linetype=2)+geom_point(color='gray25',alpha=.65,size=1.2)+facet_wrap(~stage,nrow=1)+coord_equal(xlim=c(0,1),ylim=c(0,1))+labs(title='Matern trajectories: Isomap initialization before and after MPCurve',subtitle='Matern 5/2 | length-scale = 0.2 | SNR = 4 | M = 1 | K = 50 | fixed uniform prior',x='True latent position',y='Estimated latent position',caption='Figure 1. Left: scaled Isomap coordinates. Right: posterior mean positions. rho is absolute Spearman correlation.\nInitial ELBO follows quantile/variational initialization. One realization, seed 20260930; relative tolerance 1e-6.')
ggsave(file.path(out,'ordering.png'),p,width=11,height=6.3,dpi=170,bg='white')
ord<-order(as.integer(sub('V','',colnames(X))))
long<-do.call(rbind,lapply(ord,function(j)data.frame(truth=truth,before=before,after=after,value=X[,j],feature=colnames(X)[j])))
long$feature<-factor(long$feature,levels=colnames(X)[ord])
for(stage in c('before','after')) {
 q<-ggplot(long,aes(x=.data[[stage]],y=value,color=truth))+geom_point(size=.8,alpha=.65)+facet_wrap(~feature,ncol=4)+scale_color_viridis_c(limits=c(0,1),name='True latent\nposition')+scale_x_continuous(limits=c(0,1),breaks=seq(0,1,.25))+labs(title=paste('Matern observations along',if(stage=='before')'Isomap initialization' else 'MPCurve inferred positions'),subtitle='Matern 5/2 | length-scale = 0.2 | SNR = 4',x=if(stage=='before')'Estimated t (scaled Isomap coordinate)' else 'Estimated t (posterior mean)',y='Observed feature value',caption='Figure 2. All 300 observed values per feature, colored by true latent position; no fitted smoother is applied.')
 ggsave(file.path(out,paste0('columns_',stage,'.png')),q,width=12,height=8,dpi=170,bg='white')
}
lines<-do.call(rbind,lapply(ord,function(j)data.frame(truth=grid,value=smooth[,j],feature=colnames(X)[j])))
lines$feature<-factor(lines$feature,levels=colnames(X)[ord])
q<-ggplot(long,aes(truth,value))+geom_point(color='gray55',alpha=.28,size=.55)+geom_line(data=lines,color='#0072B2',linewidth=.8)+facet_wrap(~feature,ncol=4)+labs(title='All 12 Matern generating trajectories',subtitle='Matern 5/2 | length-scale = 0.2 | seed = 20260930 | SNR = 4',x='True latent position',y='Feature value',caption='Figure 3. Blue: noiseless GP draws; gray: observations. Each signal is scaled to unit sample variance.\nLatent positions and Gaussian noise are reused from rich_S4_r06_A; no monotone anchor or outcome filtering.')
ggsave(file.path(out,'trajectories.png'),q,width=12,height=8,dpi=170,bg='white')
saveRDS(list(X=X,signal=signal,grid=grid,smooth=smooth,truth=truth,seed=seed,nu=2.5,length_scale=ell,jitter=1e-10,fit=f,before=before,after=after,rho=rho,elbo=elbo,source_hash=digest(file='experiments/m1_local_k_v034/matern_example.R',algo='sha256'),original_input_hash=d$input_hash,package_version=as.character(packageVersion('MPCurver'))),file.path(out,'result.rds'))
print(data.frame(stage=c('initial','converged'),rho=rho,elbo=elbo,iterations=c(0,f$iter)))
