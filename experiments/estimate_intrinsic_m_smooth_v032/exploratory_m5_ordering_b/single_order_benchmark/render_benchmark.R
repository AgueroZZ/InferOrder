# A fixed single-ordering benchmark using each package's native initialization.
source('experiments/estimate_intrinsic_m_smooth_v032/auto_m_v034/common.R')
suppressPackageStartupMessages(library(ggplot2))
base <- file.path(study_dir,'exploratory_m5_ordering_b')
external <- file.path(base,'external_methods')
out <- file.path(base,'single_order_benchmark')
cases <- c('B_noise0','B_noise05','B_noise1')
methods <- c('MPCurve','P-curve','GPLVM','Bayesian GPLVM')
records <- list(); sources <- list()
load_checked <- function(path,id,gp=FALSE) {
 d <- if(gp) jsonlite::fromJSON(path) else readRDS(path)
 input <- file.path(external,'inputs',paste0(id,'_X.csv'))
 expected <- if(gp) d$input_file_sha256 else if(!is.null(d$input_sha256)) d$input_sha256 else d$input_file_sha256
 stopifnot(identical(expected,digest::digest(file=input,algo='sha256')))
 sources[[length(sources)+1L]] <<- data.frame(case=id,result_path=path,
  result_sha256=digest::digest(file=path,algo='sha256'),input_path=input,input_sha256=expected)
 d
}
add <- function(id,method,protocol,d,status,version,initialization,iterations=NA_integer_) {
 key <- paste(id,method,protocol,sep=':')
 stopifnot(length(d$initial)==300L,length(d$final)==300L,length(d$truth)==300L,
  all(is.finite(c(d$initial,d$final,d$truth))))
 records[[key]] <<- list(id=id,method=method,protocol=protocol,initial=d$initial,
  final=d$final,truth=d$truth,status=status,version=version,initialization=initialization,iterations=iterations)
}
for(i in seq_along(cases)) {
 id <- cases[i]
 for(protocol in c('package_stopping','extended_stopping')) {
  d <- load_checked(file.path(out,paste0(id,'_mpcurve_',protocol,'.rds')),id)
  stopifnot(d$fit$intrinsic_dim==1L,d$fit$fit$control$method=='PCA',
   d$fit$fit$control$gamma_preinit=='ordering_init')
  add(id,'MPCurve',protocol,d,if(d$fit$converged)'Converged' else 'Iteration limit',
   d$package_version,'Centered, unscaled PCA',d$fit$fit$iter)
  setting <- if(protocol=='package_stopping')'package_defaults' else 'extended'
  d <- load_checked(file.path(external,'princurve_results',paste0(id,'_pca_',setting,'.rds')),id)
  stopifnot(is.null(d$start_curve),d$start_method=='pca')
  add(id,'P-curve',protocol,d,if(d$fit$converged)'Converged' else 'Iteration limit',
   d$package_version,'Centered, unscaled PCA',d$fit$num_iterations)
 }
 for(method in c('GPLVM','Bayesian GPLVM')) {
  index <- if(method=='GPLVM')i else i+7L
  file <- list.files(file.path(external,'gpy_results'),pattern=sprintf('^%02d_.*json$',index),full.names=TRUE)
  stopifnot(length(file)==1L)
  d <- load_checked(file,id,gp=TRUE)
  stopifnot(d$job$start=='pca',d$job$seed==20260929,
   d$job$inducing==if(method=='GPLVM')0L else 10L)
  add(id,method,'extended_stopping',d,d$optimizer_status,d$gpy_version,
   'Native GPy PCA (features standardized)')
 }
}
rank01 <- function(x) (rank(x,ties.method='average')-1)/(length(x)-1)
metrics <- function(x,truth) {
 rho <- cor(x,truth,method='spearman')
 if(rho<0)x <- -x
 upper <- upper.tri(matrix(0,length(x),length(x)))
 products <- (sign(outer(x,x,'-'))*sign(outer(truth,truth,'-')))[upper]
 c(rho=abs(rho),pair_error=mean((1-products)/2),rank_mae=mean(abs(rank01(x)-rank01(truth))))
}
summary <- do.call(rbind,lapply(records,function(d) {
 initial <- metrics(d$initial,d$truth); final <- metrics(d$final,d$truth)
 data.frame(case=d$id,noise_scale=c(0,.5,1)[match(d$id,cases)],method=d$method,
  protocol=d$protocol,initialization=d$initialization,package_version=d$version,
  initial_rho=initial['rho'],final_rho=final['rho'],
  initial_pair_error=initial['pair_error'],final_pair_error=final['pair_error'],
  initial_rank_mae=initial['rank_mae'],final_rank_mae=final['rank_mae'],
  initial_final_rho=abs(cor(d$initial,d$final,method='spearman')),
  status=d$status,iterations=d$iterations,row.names=NULL)
}))
main <- summary[summary$protocol=='extended_stopping',]
write.csv(main,file.path(out,'benchmark_summary.csv'),row.names=FALSE)
write.csv(summary,file.path(out,'all_stopping_results.csv'),row.names=FALSE)
write.csv(do.call(rbind,sources),file.path(out,'source_manifest.csv'),row.names=FALSE)
saveRDS(records,file.path(out,'benchmark_positions.rds'))
make_samples <- function(package_stop=FALSE) {
 selected <- Filter(function(d) d$protocol==if(package_stop && d$method %in% c('MPCurve','P-curve'))
  'package_stopping' else 'extended_stopping',records)
 do.call(rbind,lapply(selected,function(d) {
  initial <- d$initial; if(cor(initial,d$truth,method='spearman')<0)initial <- -initial
  final <- d$final; if(cor(initial,final,method='spearman')<0) final <- -final
  data.frame(case=d$id,method=factor(d$method,levels=methods),
   noise=factor(d$id,levels=cases,labels=c('No noise','Half original noise','Original noise (SNR 4)')),
   sample=seq_along(initial),truth=d$truth,initial_rank=rank01(initial),final_rank=rank01(final),status=d$status)
 }))
}
plot_pairs <- function(d,title,subtitle,caption) {
 labels <- do.call(rbind,lapply(split(d,interaction(d$method,d$noise,drop=TRUE)),function(x) {
  stopifnot(!anyDuplicated(x$sample))
  data.frame(method=x$method[1],noise=x$noise[1],
   label=paste0(sprintf('|rho|: %.3f -> %.3f',abs(cor(x$truth,x$initial_rank,method='spearman')),
   abs(cor(x$truth,x$final_rank,method='spearman'))),if(x$status[1]!='Converged')'\nBudget reached' else ''))
 }))
 long <- rbind(transform(d,stage='Initial',position=initial_rank),transform(d,stage='Final',position=final_rank))
 ggplot()+geom_segment(data=d,aes(x=truth,xend=truth,y=initial_rank,yend=final_rank),
  color='#888888',alpha=.2,linewidth=.2)+
  geom_point(data=long,aes(truth,position,color=stage,shape=stage),size=1,alpha=.7)+
  geom_label(data=labels,aes(x=.02,y=1.19,label=label),hjust=0,vjust=1,size=3.1,linewidth=.15)+
  facet_grid(method~noise)+scale_color_manual(values=c(Initial='#0072B2',Final='#D55E00'),breaks=c('Initial','Final'))+
  scale_shape_manual(values=c(Initial=1,Final=16),breaks=c('Initial','Final'))+
  scale_y_continuous(breaks=c(0,.5,1))+coord_cartesian(xlim=c(0,1),ylim=c(0,1.21))+
  theme_minimal(base_size=12)+theme(legend.position='bottom',panel.grid.minor=element_blank(),
   strip.text.y=element_text(size=12),plot.caption=element_text(hjust=0,size=9))+
  labs(title=title,subtitle=subtitle,x='True latent position',y='Normalized sample rank',color=NULL,shape=NULL,caption=caption)
}
caption <- paste('All methods fit one ordering to the same 300 samples and 12 features; each uses its native PCA start.',
 'Blue: initial ranks; orange: final ranks; gray lines link the same sample. Global reversal is aligned away.',
 'MPCurve: up to 2000 sweeps, tol = 1e-6. P-curve: maxit = 1000, thresh = 1e-6; default df = 5.',
 'Both GP models use native optimization with up to three 2000-evaluation blocks; Bayesian GPLVM uses the default 10 inducing points.',
 'The two noiseless GP fits reached their evaluation budgets; those endpoints are not reported as converged.',sep='\n')
samples <- make_samples()
p <- plot_pairs(samples,'B as a single-ordering benchmark: native initialization',
 'MPCurver 0.3.4 | princurve 2.1.6 | GPy 1.13.2; default initialization, extended optimization budgets.',caption)
for(ext in c('png','pdf')) ggsave(file.path(out,paste0('native_initialization_benchmark.',ext)),p,width=13,height=12,dpi=180,bg='white')
write.csv(samples,file.path(out,'benchmark_rank_positions.csv'),row.names=FALSE)
p_default <- plot_pairs(make_samples(TRUE),'Stopping sensitivity with the same native initialization',
 'MPCurve and P-curve use their package-default stopping settings; GP endpoints are unchanged.',
 paste('MPCurve: default 100 sweeps, tol = 1e-6. P-curve: default maxit = 10, thresh = 0.001.',
 'All four methods still fit only the same 12 features. Colors and rank conventions match the main figure.',
 'Bayesian GPLVM uses 10 inducing points. GP optimization budgets match the main figure.',sep='\n'))
for(ext in c('png','pdf')) ggsave(file.path(out,paste0('package_stopping_sensitivity.',ext)),p_default,width=13,height=12,dpi=180,bg='white')
stopifnot(nrow(main)==12L,nrow(summary)==18L,nrow(samples)==3600L,
 all(table(main$method)==3L),all(main$method %in% methods))
for(id in cases) {
 x <- as.matrix(read.csv(file.path(external,'inputs',paste0(id,'_X.csv'))))
 stopifnot(identical(dim(x),c(300L,12L)))
 ref <- records[[paste(id,'MPCurve','extended_stopping',sep=':')]]$initial
 stopifnot(max(abs(MPCurver:::PCA_ordering(x)$t-ref))<1e-12)
}
writeLines(c('12 primary fits; 18 total endpoints; 300 finite positions per fit.',
 'All source fit input SHA-256 hashes match the same per-case matrix.',
 'MPCurve uses the public intrinsic_dim=1 PCA initializer; no injected coordinates.',
 'P-curve starts are NULL. GP starts are native PCA; Bayesian inducing count is 10.',
 'All main runs at half/original noise converged. Both zero-noise GP runs exhausted the budget.'),
 file.path(out,'verification.txt'))
print(main[,c('case','method','initial_rho','final_rho','status')],row.names=FALSE)
