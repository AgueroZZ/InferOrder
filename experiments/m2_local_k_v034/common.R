# Prespecified fixed-M study of within-group one-step neighborhood selection.
options(stringsAsFactors=FALSE)
.libPaths(c('experiments/estimate_intrinsic_m_smooth_v032/auto_m_v034/library',
 'experiments/estimate_intrinsic_m_smooth_v032/library',.libPaths()))
suppressPackageStartupMessages({library(MPCurver);library(digest);library(mclust);library(clue)})
stopifnot(as.character(packageVersion('MPCurver'))=='0.3.4')
study <- 'experiments/m2_local_k_v034'
design <- list(n=300L,D=24L,M=2L,features_per_group=12L,K=50L,
 families=c('broad','rich'),snr=c(1,4,16),replicates=6L,k_candidates=c(5L,10L,15L,20L,30L),
 default_k=15L,anchor_features=0L,screening_sweeps=1L,warmup_sweeps=2L,
 failure_threshold=.9,recovery_threshold=.95,material_drop=.05,
 package_commit='15f2b0bbe5dfa61cd46da5160b2bc251e75a0475',
 archive_sha256='58d0c99c6170994c82eedba190fb1c28a63eb47530c575705323c3cf3992b7c6')
manifest <- expand.grid(family=design$families,replicate=1:6,snr=design$snr)
manifest$id <- sprintf('%s_S%d_r%02d',manifest$family,manifest$snr,manifest$replicate)
manifest$seed <- 730000L+1000L*match(manifest$family,design$families)+manifest$replicate
manifest$fit_seed <- manifest$seed+1000000L
atomic_save <- function(x,path) {tmp<-paste0(path,'.tmp-',Sys.getpid());saveRDS(x,tmp);stopifnot(file.rename(tmp,path))}
generate_data <- function(row) {
 set.seed(row$seed)
 latent <- matrix(runif(600),300,2,dimnames=list(NULL,c('A','B')))
 frequencies <- if(row$family=='broad') 2:4 else 2:6
 attenuation <- frequencies^if(row$family=='broad') -2 else -1.5
 signal <- matrix(NA_real_,300,24); coefficient <- vector('list',24)
 truth <- rep(1:2,each=12)
 for(j in 1:24) {
  a<-rnorm(length(frequencies));b<-rnorm(length(frequencies));angles<-outer(latent[,truth[j]],pi*frequencies)
  raw <- as.numeric(sin(angles)%*%(a*attenuation)+cos(angles)%*%(b*attenuation))
  signal[,j]<-(raw-mean(raw))/sd(raw)
  coefficient[[j]]<-list(sine=a,cosine=b,center=mean(raw),scale=sd(raw))
 }
 noise <- matrix(rnorm(300*24),300,24)
 permutation <- sample.int(24)
 signal<-signal[,permutation];noise<-noise[,permutation];truth<-truth[permutation]
 colnames(signal)<-colnames(noise)<-paste0('V',permutation)
 X<-signal+noise/sqrt(row$snr)
 list(row=row,X=X,signal=signal,unit_noise=noise,latent=latent,truth=truth,
  frequencies=frequencies,attenuation=attenuation,coefficients=coefficient[permutation],
  input_hash=digest(X,algo='sha256'),signal_hash=digest(signal,algo='sha256'),
  truth_hash=digest(list(latent,truth),algo='sha256'))
}
local_fit <- function(X,position,iterations) MPCurver:::.cavi_build_from_ordering(
 X=X,ordering_vec=position,K=50L,rw_q=2L,ridge=0,lambda_init=1,
 max_iter=iterations,tol=1e-6,discretization='quantile',strict_K=TRUE)
expand_fit <- function(X,fit) MPCurver:::cavi(X=X,K=50L,
 responsibilities_init=fit$gamma,position_prior_init=colMeans(fit$gamma),
 rw_q=2L,ridge=0,max_iter=0L,convergence='relative',verbose=FALSE)
fit_joint <- function(X,fits,seed) {
 set.seed(seed)
 f<-fit_mpcurve(X=X,intrinsic_dim=2L,algorithm='cavi',method='isomap',fits_init=fits,
  K=50L,rw_q=2L,ridge=0,lambda=1,fix_lambda=FALSE,partition_prior='adaptive',
  position_prior='adaptive',discretization='quantile',iter=1500L,num_cores=1L,
  tol=1e-6,convergence='normalized',T_start=5,T_end=1,n_outer=25L,inner_iter=1L,
  max_converge_iter=1500L,tol_outer=1e-6,verbose=FALSE)
 while(!isTRUE(f$fit$converged) && f$fit$iter<9999L)
  f<-do_mpcurve(f,iter=min(1500L,9999L-f$fit$iter),tol=1e-6,tol_outer=1e-6,convergence='normalized',verbose=FALSE)
 f
}
