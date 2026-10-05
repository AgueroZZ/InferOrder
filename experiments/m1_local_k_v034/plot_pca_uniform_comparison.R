# Two-panel display of saved PCA and fixed-uniform MPCurve estimates.
source('experiments/m2_local_k_v034/common.R')
library(ggplot2)
out <- 'experiments/m1_local_k_v034'
a <- readRDS(file.path(out,'failed_case_pca_fit.rds'))
u <- readRDS(file.path(out,'failed_case_pca_uniform_fit.rds'))
stopifnot(identical(a$input_hash,u$input_hash))
initial_elbo <- head(u$fit$elbo_trace,1)
final_elbo <- tail(u$fit$elbo_trace,1)
stopifnot(length(u$fit$elbo_trace)==u$fit$iter+1L)
labels <- c(sprintf('PCA initialization (iteration 0)\n|Spearman rho| = %.3f | ELBO = %.2f',a$rho['before'],initial_elbo),sprintf('MPCurve: fixed uniform prior (iteration %d)\n|Spearman rho| = %.3f | ELBO = %.2f',u$fit$iter,u$summary$rho[u$summary$method=='Uniform'],final_elbo))
pts <- data.frame(truth=rep(a$truth,2),estimate=c(a$before,u$position),stage=factor(rep(labels,each=length(a$truth)),levels=labels))
p <- ggplot(pts,aes(truth,estimate))+geom_abline(slope=1,intercept=0,color='gray70',linetype=2)+geom_point(color='gray25',alpha=.65,size=1.2)+facet_wrap(~stage,nrow=1)+coord_equal(xlim=c(0,1),ylim=c(0,1))+theme_minimal(base_size=12)+labs(title='PCA initialization before and after MPCurve',subtitle=sprintf('rich_S4_r06_A | M = 1 | K = 50 | fixed uniform position prior | %d updates',u$fit$iter),x='True latent position',y='Estimated latent position',caption='Figure 1. Left: PC scores scaled to [0,1]. Right: posterior mean position under the fixed uniform prior.\nInitial ELBO is evaluated after PCA quantile discretization and variational initialization, before CAVI sweeps.\nBoth panels share orientation; convergence uses relative ELBO tolerance 1e-6.')
ggsave(file.path(out,'failed_case_pca_uniform_two_panel.png'),p,width=11,height=6.3,dpi=170,bg='white')
