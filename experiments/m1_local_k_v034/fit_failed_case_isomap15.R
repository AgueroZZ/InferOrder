# Isomap k=15 with the same fixed-uniform single-ordering controls as PCA.
source('experiments/m2_local_k_v034/common.R')
library(ggplot2)
out <- 'experiments/m1_local_k_v034'
d <- readRDS(file.path(study,'data','rich_S4_r06.rds'))
X <- d$X[,d$truth==1,drop=FALSE];truth <- d$latent[,1]
raw <- as.numeric(MPCurver:::isomap_ordering(X,k=15L)$t)
init <- MPCurver:::.cavi_exact_k_init_from_ordering(X=X,ordering_vec=raw,K=50L,discretization='quantile',strict_K=TRUE)
gamma0 <- MPCurver:::.cavi_hard_gamma_from_cluster_rank(init$cluster_rank,nrow(X),50L)
f <- MPCurver:::cavi(X=X,K=50L,responsibilities_init=gamma0,position_prior_init=rep(1/50,50),position_prior='fixed',sigma2_init=init$sigma2,lambda_init=rep(1,ncol(X)),rw_q=2L,ridge=0,discretization='quantile',max_iter=2000L,tol=1e-6,convergence='relative',verbose=FALSE)
while(!f$converged && f$iter<10000L)f<-MPCurver:::do_cavi(f,iter=2000L,tol=1e-6)
stopifnot(f$converged,max(abs(f$params$pi-1/50))<1e-12)
before <- (raw-min(raw))/diff(range(raw));after <- as.numeric(f$gamma%*%seq(0,1,length.out=50))
if(cor(before,truth,method='spearman')<0){before<-1-before;after<-1-after}
rho <- c(abs(cor(before,truth,method='spearman')),abs(cor(after,truth,method='spearman')))
elbo <- c(head(f$elbo_trace,1),tail(f$elbo_trace,1))
labels <- c(sprintf('Isomap k=15 (iteration 0)\n|Spearman rho| = %.3f | ELBO = %.2f',rho[1],elbo[1]),sprintf('MPCurve: fixed uniform prior (iteration %d)\n|Spearman rho| = %.3f | ELBO = %.2f',f$iter,rho[2],elbo[2]))
pts<-data.frame(truth=rep(truth,2),estimate=c(before,after),stage=factor(rep(labels,each=300),levels=labels))
p<-ggplot(pts,aes(truth,estimate))+geom_abline(slope=1,intercept=0,color='gray70',linetype=2)+geom_point(color='gray25',alpha=.65,size=1.2)+facet_wrap(~stage,nrow=1)+coord_equal(xlim=c(0,1),ylim=c(0,1))+theme_minimal(base_size=12)+labs(title='Isomap k=15 initialization before and after MPCurve',subtitle='rich_S4_r06_A | M = 1 | K = 50 | fixed uniform position prior',x='True latent position',y='Estimated latent position',caption='Figure 1. Left: Isomap coordinates scaled to [0,1]. Right: posterior mean position; both share orientation.\nInitial ELBO follows quantile discretization and variational initialization, before CAVI sweeps. Relative tolerance: 1e-6.')
ggsave(file.path(out,'failed_case_isomap15_uniform_before_after.png'),p,width=11,height=6.3,dpi=170,bg='white')
saveRDS(list(fit=f,raw=raw,before=before,after=after,truth=truth,rho=rho,elbo=elbo,input_hash=digest(X,algo='sha256'),source_hash=digest(file=file.path(out,'fit_failed_case_isomap15.R'),algo='sha256')),file.path(out,'failed_case_isomap15_uniform_fit.rds'))
print(data.frame(stage=c('initial','converged'),rho=rho,elbo=elbo,iterations=c(0,f$iter)))
