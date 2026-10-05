# Raw PCA initialization for the previously displayed failed M=1 dataset.
source('experiments/m2_local_k_v034/common.R')
library(ggplot2)
out <- 'experiments/m1_local_k_v034'
d <- readRDS(file.path(study,'data','rich_S4_r06.rds'))
r <- readRDS(file.path(out,'results','rich_S4_r06_A.rds'))
X <- d$X[,d$truth==1,drop=FALSE]
stopifnot(identical(r$input_hash,digest(X,algo='sha256')))
truth <- d$latent[,1]
raw <- as.numeric(MPCurver:::PCA_ordering(X)$t)
# PCA direction is arbitrary; use truth only to orient the display.
if(cor(raw,truth,method='spearman')<0) raw <- -raw
ordering <- (raw-min(raw))/diff(range(raw))
rho <- cor(ordering,truth,method='spearman')
points <- data.frame(truth=truth,ordering=ordering)
theme_set(theme_minimal(base_size=12))
p <- ggplot(points,aes(truth,ordering,color=truth))+geom_abline(slope=1,intercept=0,color='gray75',linetype=2)+geom_point(size=1.5,alpha=.8)+scale_color_viridis_c(limits=c(0,1),name='True latent\nposition')+coord_equal()+labs(title='Raw PCA initialization versus true latent position',subtitle=sprintf('rich_S4_r06_A | 300 samples | 12 features | Spearman correlation = %.3f',rho),x='True latent position',y='Estimated t (PC scores scaled to [0,1])',caption='Figure 1. PC scores linearly scaled to [0,1], before MPCurve updates.\nPCA direction is aligned with truth for display. Dashed line: identity.')
ggsave(file.path(out,'failed_case_pca_ordering.png'),p,width=7.5,height=6,dpi=170,bg='white')
cols <- order(as.integer(sub('V','',colnames(X))))
long <- do.call(rbind,lapply(cols,function(j)data.frame(ordering=ordering,value=X[,j],truth=truth,feature=colnames(X)[j])))
long$feature <- factor(long$feature,levels=colnames(X)[cols])
p2 <- ggplot(long,aes(ordering,value,color=truth))+geom_point(size=.8,alpha=.65)+facet_wrap(~feature,ncol=4)+scale_color_viridis_c(limits=c(0,1),name='True latent\nposition')+labs(title='All 12 observed columns along the raw PCA ordering',subtitle=sprintf('rich_S4_r06_A | SNR = 4 | PCA ordering recovery = %.3f',rho),x='Estimated t (PC scores scaled to [0,1])',y='Observed feature value',caption='Figure 2. Each panel contains all 300 observed values, placed along the same PCA ordering as Figure 1.\nColors indicate true latent position; no fitted smoother or MPCurve update is applied.')
ggsave(file.path(out,'failed_case_columns_on_pca.png'),p2,width=12,height=8,dpi=170,bg='white')
write.csv(data.frame(sample_index=seq_len(nrow(X)),true_position=truth,pca_position_aligned=raw,pca_ordering=ordering),file.path(out,'failed_case_pca_positions.csv'),row.names=FALSE)
cat('PCA Spearman correlation:',rho,'\n')
