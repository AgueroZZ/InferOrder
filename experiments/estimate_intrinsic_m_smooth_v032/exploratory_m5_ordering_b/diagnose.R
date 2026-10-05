# Isolated diagnostic of the current-package M=5, SNR=4 example.
source('experiments/estimate_intrinsic_m_smooth_v032/auto_m_v034/common.R')
suppressPackageStartupMessages({library(ggplot2); library(patchwork)})
out <- file.path(study_dir, 'exploratory_m5_ordering_b')
id <- 'main_M5_S4_r001'
d <- load_fixed_dataset(manifest[manifest$id == id, ])
r <- readRDS(file.path(extension_dir,'results',paste0(id,'_auto_adaptive.rds')))
s <- readRDS(file.path(extension_dir,'full_fits',paste0(id,'.rds')))
rho <- function(x,y) abs(cor(x,y,method='spearman'))
set.seed(d$row$seed + design$seed_offset)
ini <- MPCurver:::init_m_trajectories_cavi(X=d$X, M=5L, methods=rep('isomap',5), K=50L,
 rw_q=2L, ridge=0, discretization='quantile', partition_init='similarity',
 similarity_metric='spline_r2', spline_r2_df=5L, cluster_linkage='single', num_iter=2L)
grid <- seq(0,1,length.out=50)
rows <- list(); points <- list()
for(m in 1:5) {
 label <- d$ordering_labels[m]; jj <- which(d$true_assign == label)
 slot <- which(vapply(ini$init_info,function(z) setequal(z$feature_idx,jj),logical(1)))
 stopifnot(length(slot)==1L, setequal(s$fit$init_info[[slot]]$feature_idx,jj))
 t <- d$latent_positions[,m]; raw <- MPCurver:::isomap_ordering(d$X[,jj])$t
 initial <- as.numeric(ini$fits[[slot]]$gamma %*% grid)
 final <- r$compact$positions[,slot]
 rows[[m]] <- data.frame(ordering=label,slot=slot,features=length(jj),
  raw_isomap=rho(t,raw),initial=rho(t,initial),final=rho(t,final),
  initial_final=rho(initial,final),noise_variance=mean(d$realized_noise_variance[jj]))
 for(stage in c('raw','initial','final')) {
  y <- get(stage); if(cor(t,y)<0) y <- max(y)+min(y)-y
  points[[length(points)+1L]] <- data.frame(truth=t,position=y,ordering=label,stage=stage)
 }
}
summary <- do.call(rbind,rows); print(summary)
write.csv(summary,file.path(out,'ordering_summary.csv'),row.names=FALSE)
b <- which(d$true_assign=='B'); t <- d$latent_positions[,2]; anchor <- intersect(b,d$anchor_indices)
print(data.frame(anchor=colnames(d$X)[anchor],rho=rho(t,d$X[,anchor])))
noise <- d$X[,b]-d$signal[,b]
geometry <- function(X,k) {
 nn <- RANN::nn2(X,X,k=k+1L)
 delta <- abs(t[row(nn$nn.idx[,-1,drop=FALSE])]-t[nn$nn.idx[,-1,drop=FALSE]])
 edges <- data.frame(from=rep(seq_len(nrow(X)),each=k),to=as.vector(t(nn$nn.idx[,-1,drop=FALSE])),weight=as.vector(t(nn$nn.dists[,-1,drop=FALSE])))
 g <- igraph::graph_from_data_frame(edges,directed=FALSE,vertices=seq_len(nrow(X)))
 gd <- igraph::distances(g); td <- as.matrix(dist(t)); up <- upper.tri(td)
 pos <- MPCurver:::isomap_ordering(X,k=k)$t
 list(pos=pos,stats=data.frame(rho=rho(t,pos),components=igraph::components(g)$no,
  geodesic_rho=rho(gd[up],td[up]),long_edge_fraction=mean(delta>0.25),max_edge_span=max(delta)))
}
settings <- expand.grid(noise_scale=c(0,0.5,1),k=c(5L,10L,15L,20L,30L))
sensitivity <- list(); bpoints <- list()
for(i in seq_len(nrow(settings))) {
 z <- settings[i,]; g <- geometry(d$signal[,b]+z$noise_scale*noise,z$k)
 sensitivity[[i]] <- cbind(z,g$stats)
 y <- g$pos; if(cor(t,y,use='complete.obs')<0) y <- 1-y
 bpoints[[i]] <- data.frame(truth=t,position=y,noise_scale=z$noise_scale,k=z$k)
}
sensitivity <- do.call(rbind,sensitivity); print(sensitivity)
write.csv(sensitivity,file.path(out,'geometry_sensitivity.csv'),row.names=FALSE)
point_frame <- do.call(rbind,points)
point_frame$stage <- factor(point_frame$stage,levels=c('raw','initial','final'))
p <- ggplot(point_frame,aes(truth,position))+geom_point(size=.65,alpha=.5,color='#0072B2')+
 facet_grid(stage~ordering,scales='free_y')+theme_minimal(base_size=11)+labs(x='True latent position',y='Inferred position (orientation aligned)',title='M=5, SNR=4, replicate 1: initialization and final fit',subtitle='Feature groups are exact; raw = Isomap, initial = two initialization CAVI sweeps')
ggsave(file.path(out,'initial_vs_final.png'),p,width=13,height=7,dpi=160)
p2 <- ggplot(do.call(rbind,bpoints),aes(truth,position))+geom_point(size=.65,alpha=.5,color='#009E73',na.rm=TRUE)+
 facet_grid(noise_scale~k,labeller=label_both)+theme_minimal(base_size=11)+labs(x='True B position',y='Isomap position (orientation aligned)',title='Ordering B: sensitivity to neighborhood size and measurement noise',subtitle='Same 12 B features and same noise realization; noise scale 1 is the observed data',
 caption='Noiseless k=5: disconnected graph; only the largest component (244/300 samples) is plotted. Other panels include all samples.')
ggsave(file.path(out,'geometry_sensitivity.png'),p2,width=13,height=7,dpi=160)
features <- do.call(rbind,lapply(b,function(j) data.frame(truth=t,signal=d$signal[,j],observed=d$X[,j],feature=paste0(colnames(d$X)[j],if(j%in%anchor)' (anchor)' else ''))))
p3 <- ggplot(features,aes(truth,observed))+geom_point(alpha=.25,size=.5)+geom_line(aes(y=signal),color='#D55E00')+facet_wrap(~feature,ncol=4)+theme_minimal()+labs(x='True B position',y='Feature value',title='All 12 features in true ordering B',subtitle='Black: observations; orange: noiseless trajectories')
ggsave(file.path(out,'b_features.png'),p3,width=12,height=8,dpi=160)
saveRDS(list(dataset=id,input_hash=d$input_hash,provenance=r$provenance,seed=d$row$seed+design$seed_offset,session=sessionInfo()),file.path(out,'provenance.rds'))
# Identify the actual long-range k=15 graph edges; threshold is diagnostic only.
nn <- RANN::nn2(d$X[,b],d$X[,b],k=16L)
edges <- data.frame(from=rep(1:300,each=15),to=as.vector(t(nn$nn.idx[,-1])),weight=as.vector(t(nn$nn.dists[,-1])))
edges$true_from <- t[edges$from]; edges$true_to <- t[edges$to]
edges$span <- abs(edges$true_from-edges$true_to)
write.csv(edges[edges$span>.25,],file.path(out,'long_range_edges.csv'),row.names=FALSE)
if ('--diagnostics-only' %in% commandArgs(trailingOnly=TRUE)) quit(status=0L)
# Same full-data structural fit, changing only B's neighborhood size to 10.
slot_b <- summary$slot[summary$ordering=='B']
controlled <- list()
for(k in c(15L,10L)) {
 fits <- ini$fits
 if(k==10L) {
  pos <- MPCurver:::isomap_ordering(d$X[,b],k=k)$t
  sub <- MPCurver:::.cavi_build_from_ordering(X=d$X[,b],ordering_vec=pos,S=NULL,K=50L,rw_q=2L,ridge=0,
   lambda_sd_prior_rate=NULL,lambda_min=1e-10,lambda_max=1e10,sigma_min=1e-10,sigma_max=1e10,
   max_iter=2L,tol=1e-6,discretization='quantile',strict_K=TRUE)
  fits[[slot_b]] <- MPCurver:::cavi(X=d$X,K=50L,responsibilities_init=sub$gamma,
   position_prior_init=colMeans(sub$gamma),rw_q=2L,ridge=0,max_iter=0L,convergence='relative',verbose=FALSE)
 }
 set.seed(d$row$seed+design$seed_offset)
 f <- MPCurver::fit_mpcurve(X=d$X,intrinsic_dim=5L,algorithm='cavi',method='isomap',K=50L,rw_q=2L,ridge=0,
  fits_init=fits,partition_prior='adaptive',position_prior='adaptive',
  lambda=1,fix_lambda=FALSE,discretization='quantile',num_cores=1L,iter=1500L,tol=1e-6,convergence='normalized',
  T_start=5,T_end=1,n_outer=25L,inner_iter=1L,max_converge_iter=1500L,tol_outer=1e-6,verbose=FALSE)
 positions <- vapply(f$fit$fits, function(z) as.numeric(z$gamma %*% grid), numeric(300))
 controlled[[length(controlled)+1L]] <- data.frame(B_neighbors=k,B_rho=rho(t,positions[,slot_b]),
  objective=tail(f$fit$objective_history,1),iterations=f$fit$iter,converged=f$fit$converged,
  ARI=mclust::adjustedRandIndex(d$true_assign,f$fit$assign))
 print(controlled[[length(controlled)]])
}
write.csv(do.call(rbind,controlled),file.path(out,'controlled_restarts.csv'),row.names=FALSE)
