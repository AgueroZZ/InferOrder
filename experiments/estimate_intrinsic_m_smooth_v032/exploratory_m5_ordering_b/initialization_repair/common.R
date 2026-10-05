# Controlled basin-of-attraction diagnostics; truth is used only for experiments.
source('experiments/estimate_intrinsic_m_smooth_v032/auto_m_v034/common.R')
repair_dir <- file.path(study_dir, 'exploratory_m5_ordering_b', 'initialization_repair')
for (directory in c('results', 'full_fits')) dir.create(file.path(repair_dir,directory),showWarnings=FALSE)
dataset_id <- 'main_M5_S4_r001'
dataset <- load_fixed_dataset(manifest[manifest$id == dataset_id, ])
B <- which(dataset$true_assign == 'B')
truth <- dataset$latent_positions[,2]
anchor <- intersect(B,dataset$anchor_indices)
fit_seed <- dataset$row$seed + design$seed_offset
grid <- seq(0,1,length.out=50)
rank01 <- function(x) (rank(x,ties.method='average')-1)/(length(x)-1)
align <- function(x) if (cor(truth,x,method='spearman') < 0) -x else x
truth_ranks <- rank01(truth)
pair_mask <- upper.tri(matrix(0,300,300))
truth_pair_sign <- sign(outer(truth,truth,'-'))[pair_mask]
metrics <- function(position) {
 position <- align(position)
 pred_rank <- rank01(position)
 c(rho=cor(truth,position,method='spearman'),
  pair_error=mean((1-truth_pair_sign*sign(outer(position,position,'-'))[pair_mask])/2),
  mean_rank_error=mean(abs(pred_rank-truth_ranks)))
}
make_initial_fits <- function() {
 set.seed(fit_seed)
 MPCurver:::init_m_trajectories_cavi(X=dataset$X,M=5L,methods=rep('isomap',5),
  K=50L,rw_q=2L,ridge=0,discretization='quantile',partition_init='similarity',
  similarity_metric='spline_r2',spline_r2_df=5L,cluster_linkage='single',num_iter=2L)
}
scenarios <- data.frame(id=c('oracle','oracle_reversed','isomap10','isomap15','pca',
 paste0('jitter_',rep(c('010','030','060'),each=3),'_r',rep(1:3,3)),
 'swap_middle_blocks','reverse_middle_half','circular_shift','random_r1','random_r2'),
 type=c('oracle','oracle_reversed','isomap10','isomap15','pca',rep('jitter',9),
 'swap_middle_blocks','reverse_middle_half','circular_shift',rep('random',2)),
 magnitude=c(rep(NA_real_,5),rep(c(.1,.3,.6),each=3),rep(NA_real_,5)),
 seed=91000L+seq_len(19))
make_ordering <- function(row) {
 set.seed(row$seed)
 switch(row$type,
  oracle=truth,
  oracle_reversed=1-truth,
  isomap10=MPCurver:::isomap_ordering(dataset$X[,B],k=10L)$t,
  isomap15=MPCurver:::isomap_ordering(dataset$X[,B],k=15L)$t,
  pca=MPCurver:::PCA_ordering(dataset$X[,B])$t,
  jitter=truth + rnorm(length(truth),sd=row$magnitude),
  swap_middle_blocks={
   v <- truth_ranks; left <- v>=.25 & v<.5; right <- v>=.5 & v<.75
   v[left] <- v[left]+.25; v[right] <- v[right]-.25; v
  },
  reverse_middle_half={v <- truth_ranks; mid <- v>=.25 & v<=.75; v[mid] <- 1-v[mid]; v},
  circular_shift=(truth_ranks+.25) %% 1,
  random=runif(length(truth)))
}
