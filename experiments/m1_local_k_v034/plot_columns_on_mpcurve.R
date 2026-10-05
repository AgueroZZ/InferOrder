# Display the same observed columns against converged fixed-uniform positions.
source('experiments/m2_local_k_v034/common.R')
library(ggplot2)
out <- 'experiments/m1_local_k_v034'
d <- readRDS(file.path(study,'data','rich_S4_r06.rds'))
u <- readRDS(file.path(out,'failed_case_pca_uniform_fit.rds'))
X <- d$X[,d$truth==1,drop=FALSE]
stopifnot(identical(u$input_hash,digest(X,algo='sha256')))
cols <- order(as.integer(sub('V','',colnames(X))))
long <- do.call(rbind,lapply(cols,function(j)data.frame(position=u$position,value=X[,j],truth=d$latent[,1],feature=colnames(X)[j])))
long$feature <- factor(long$feature,levels=colnames(X)[cols])
theme_set(theme_minimal(base_size=12))
p <- ggplot(long,aes(position,value,color=truth))+geom_point(size=.8,alpha=.65)+facet_wrap(~feature,ncol=4)+scale_color_viridis_c(limits=c(0,1),name='True latent\nposition')+scale_x_continuous(limits=c(0,1),breaks=seq(0,1,.25))+labs(title='All 12 observed columns along the MPCurve inferred positions',subtitle=sprintf('rich_S4_r06_A | SNR = 4 | fixed uniform prior | recovery = %.3f',u$summary$rho[u$summary$method=='Uniform']),x='Estimated t (MPCurve posterior mean)',y='Observed feature value',caption='Figure 2. The same 300 observations per column, now placed along converged MPCurve posterior mean positions.\nColors indicate true latent position. Layout and axis ranges match the preceding PCA display; no smoother is applied.')
ggsave(file.path(out,'failed_case_columns_on_mpcurve_uniform.png'),p,width=12,height=8,dpi=170,bg='white')
