# Controlled comparison: only freeze the position prior at uniform weights.
source('experiments/m2_local_k_v034/common.R')
library(ggplot2)
out <- 'experiments/m1_local_k_v034'
d <- readRDS(file.path(study,'data','rich_S4_r06.rds'))
a <- readRDS(file.path(out,'failed_case_pca_fit.rds'))
X <- d$X[,d$truth==1,drop=FALSE];truth<-d$latent[,1]
stopifnot(identical(a$input_hash,digest(X,algo='sha256')))
init <- MPCurver:::.cavi_exact_k_init_from_ordering(X=X,ordering_vec=a$raw_pca,K=50L,discretization='quantile',strict_K=TRUE)
gamma0 <- MPCurver:::.cavi_hard_gamma_from_cluster_rank(init$cluster_rank,nrow(X),50L)
stopifnot(max(abs(init$pi-1/50))<1e-12)
args <- list(X=X,K=50L,responsibilities_init=gamma0,position_prior_init=rep(1/50,50),sigma2_init=init$sigma2,lambda_init=rep(1,ncol(X)),rw_q=2L,ridge=0,discretization='quantile',max_iter=2000L,tol=1e-6,convergence='relative',verbose=FALSE)
# Verify the explicit reconstruction reproduces the earlier adaptive fit.
control <- do.call(MPCurver:::cavi,c(args,list(position_prior='adaptive')))
stopifnot(max(abs(control$gamma-a$fit$gamma))<1e-10)
f <- do.call(MPCurver:::cavi,c(args,list(position_prior='fixed')))
while(!f$converged && f$iter<10000L)f<-MPCurver:::do_cavi(f,iter=2000L,tol=1e-6)
stopifnot(f$converged,max(abs(f$params$pi-1/50))<1e-12)
pos <- as.numeric(f$gamma%*%seq(0,1,length.out=50))
if(cor(a$raw_pca,truth,method='spearman')<0)pos<-1-pos
positions <- list(PCA=a$before,Adaptive=a$after,Uniform=pos)
summary <- data.frame(method=names(positions),rho=vapply(positions,function(v)abs(cor(truth,v,method='spearman')),numeric(1)),iterations=c(0,a$fit$iter,f$iter))
labels<-sprintf('%s\n|Spearman rho| = %.3f',summary$method,summary$rho)
pts<-data.frame(truth=rep(truth,3),estimate=unlist(positions,use.names=FALSE),stage=factor(rep(labels,each=300),levels=labels))
p<-ggplot(pts,aes(truth,estimate))+geom_abline(slope=1,intercept=0,color='gray70',linetype=2)+geom_point(color='gray25',alpha=.65,size=1)+facet_wrap(~stage,nrow=1)+coord_equal(xlim=c(0,1),ylim=c(0,1))+theme_minimal(base_size=12)+labs(title='PCA initialization: adaptive versus fixed uniform position prior',subtitle=sprintf('rich_S4_r06_A | K = 50 | converged in %d adaptive / %d uniform updates',a$fit$iter,f$iter),x='True latent position',y='Estimated latent position',caption='Figure 1. PCA scores are scaled to [0,1]; fitted positions are posterior grid means. All panels share the same orientation.\nBoth fits use identical initial responsibilities and parameters, quantile discretization, and relative ELBO tolerance 1e-6.')
ggsave(file.path(out,'failed_case_pca_uniform_comparison.png'),p,width=13,height=5.6,dpi=170,bg='white')
saveRDS(list(fit=f,position=pos,summary=summary,input_hash=a$input_hash,args=args,source_hash=digest(file=file.path(out,'fit_failed_case_uniform.R'),algo='sha256')),file.path(out,'failed_case_pca_uniform_fit.rds'))
write.csv(summary,file.path(out,'failed_case_pca_prior_comparison.csv'),row.names=FALSE)
print(summary);cat('Uniform maximum responsibility median:',median(apply(f$gamma,1,max)),' occupied argmax bins:',length(unique(max.col(f$gamma))),'\n')
