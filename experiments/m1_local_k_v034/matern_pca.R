# PCA comparison on the exact saved Matérn realization.
source('experiments/m2_local_k_v034/common.R')
library(ggplot2)
out<-'experiments/m1_local_k_v034/matern_example'
d<-readRDS(file.path(out,'result.rds'));X<-d$X;truth<-d$truth
raw<-as.numeric(MPCurver:::PCA_ordering(X)$t)
init<-MPCurver:::.cavi_exact_k_init_from_ordering(X=X,ordering_vec=raw,K=50L,discretization='quantile',strict_K=TRUE)
gamma0<-MPCurver:::.cavi_hard_gamma_from_cluster_rank(init$cluster_rank,300,50L)
f<-MPCurver:::cavi(X=X,K=50L,responsibilities_init=gamma0,position_prior_init=rep(1/50,50),position_prior='fixed',sigma2_init=init$sigma2,lambda_init=rep(1,12),rw_q=2L,ridge=0,discretization='quantile',max_iter=2000L,tol=1e-6,convergence='relative',verbose=FALSE)
while(!f$converged && f$iter<10000L)f<-MPCurver:::do_cavi(f,iter=2000L,tol=1e-6)
stopifnot(f$converged,max(abs(f$params$pi-.02))<1e-12)
before<-(raw-min(raw))/diff(range(raw));after<-as.numeric(f$gamma%*%seq(0,1,length.out=50))
if(cor(before,truth,method='spearman')<0){before<-1-before;after<-1-after}
rho<-c(abs(cor(before,truth,method='spearman')),abs(cor(after,truth,method='spearman')))
elbo<-c(head(f$elbo_trace,1),tail(f$elbo_trace,1))
labels<-c(sprintf('PCA initialization (iteration 0)\nrho = %.3f | ELBO = %.2f',rho[1],elbo[1]),sprintf('MPCurve (iteration %d)\nrho = %.3f | ELBO = %.2f',f$iter,rho[2],elbo[2]))
theme_set(theme_minimal(base_size=12))
p<-ggplot(data.frame(truth=rep(truth,2),estimate=c(before,after),stage=factor(rep(labels,each=300),levels=labels)),aes(truth,estimate))+geom_abline(slope=1,intercept=0,color='gray70',linetype=2)+geom_point(color='gray25',alpha=.65,size=1.2)+facet_wrap(~stage,nrow=1)+coord_equal(xlim=c(0,1),ylim=c(0,1))+labs(title='Matern trajectories: PCA initialization before and after MPCurve',subtitle='Matern 5/2 | length-scale = 0.2 | SNR = 4 | M = 1 | K = 50 | fixed uniform prior',x='True latent position',y='Estimated latent position',caption='Figure 1. Left: PC scores scaled to [0,1]. Right: posterior mean positions. rho is absolute Spearman correlation.\nInitial ELBO follows quantile/variational initialization. Same dataset as the Isomap comparison; relative tolerance 1e-6.')
ggsave(file.path(out,'pca_ordering.png'),p,width=11,height=6.3,dpi=170,bg='white')
saveRDS(list(fit=f,before=before,after=after,rho=rho,elbo=elbo,input_hash=digest(X,algo='sha256'),source_hash=digest(file='experiments/m1_local_k_v034/matern_pca.R',algo='sha256')),file.path(out,'pca_result.rds'))
print(data.frame(stage=c('initial','converged'),rho=rho,elbo=elbo,iterations=c(0,f$iter)))
