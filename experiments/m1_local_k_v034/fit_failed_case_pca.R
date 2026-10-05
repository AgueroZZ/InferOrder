# Fit the same M=1 case from raw PCA using the study's single-ordering controls.
source('experiments/m2_local_k_v034/common.R')
library(ggplot2)
out <- 'experiments/m1_local_k_v034'
d <- readRDS(file.path(study,'data','rich_S4_r06.rds'))
X <- d$X[,d$truth==1,drop=FALSE];truth <- d$latent[,1]
ref <- readRDS(file.path(out,'results','rich_S4_r06_A.rds'))
stopifnot(identical(ref$input_hash,digest(X,algo='sha256')))
raw <- as.numeric(MPCurver:::PCA_ordering(X)$t)
fit <- local_fit(X,raw,2000L)
while(!fit$converged && fit$iter<10000L)fit<-MPCurver:::do_cavi(fit,iter=2000L,tol=1e-6)
stopifnot(fit$converged,all(is.finite(fit$gamma)))
before <- (raw-min(raw))/diff(range(raw))
after <- as.numeric(fit$gamma%*%seq(0,1,length.out=50))
# Use the same orientation for both panels; truth is used only for display.
if(cor(before,truth,method='spearman')<0){before<-1-before;after<-1-after}
rho <- c(before=abs(cor(truth,before,method='spearman')),after=abs(cor(truth,after,method='spearman')))
labels <- c(sprintf('Before: PCA initialization\n|Spearman rho| = %.3f',rho[1]),sprintf('After: MPCurve (%d iterations)\n|Spearman rho| = %.3f',fit$iter,rho[2]))
pts <- data.frame(truth=rep(truth,2),estimate=c(before,after),stage=factor(rep(labels,each=length(truth)),levels=labels))
p<-ggplot(pts,aes(truth,estimate))+geom_abline(slope=1,intercept=0,color='gray70',linetype=2)+geom_point(color='gray25',alpha=.65,size=1.2)+facet_wrap(~stage,nrow=1)+coord_equal(xlim=c(0,1),ylim=c(0,1))+theme_minimal(base_size=12)+labs(title='PCA initialization before and after MPCurve',subtitle='rich_S4_r06_A | M = 1 | 300 samples | 12 features | SNR = 4',x='True latent position',y='Estimated latent position',caption='Figure 1. Left: PC scores linearly scaled to [0,1]. Right: posterior mean position on the 50-point grid.\nThe fit uses the study\'s quantile discretization and relative ELBO tolerance 1e-6. Both panels share the same orientation.')
ggsave(file.path(out,'failed_case_pca_before_after.png'),p,width=11,height=6,dpi=170,bg='white')
saveRDS(list(fit=fit,raw_pca=raw,before=before,after=after,truth=truth,rho=rho,input_hash=digest(X,algo='sha256'),source_hash=digest(file=file.path(out,'fit_failed_case_pca.R'),algo='sha256'),package_version=as.character(packageVersion('MPCurver'))),file.path(out,'failed_case_pca_fit.rds'))
print(data.frame(before=rho[1],after=rho[2],iterations=fit$iter,converged=fit$converged,initial_final_rho=cor(before,after,method='spearman')))
