source('experiments/m2_local_k_v034/common.R')
for(i in seq_len(nrow(manifest))) {
 d<-readRDS(file.path(study,'data',paste0(manifest$id[i],'.rds')))
 x<-readRDS(file.path(study,'results',paste0(manifest$id[i],'.rds')))
 stopifnot(identical(d,generate_data(manifest[i,])),identical(x$input_hash,d$input_hash),
  max(abs(apply(d$signal,2,var)-1))<1e-12,
  identical(x$common_hash,digest(file=file.path(study,'common.R'),algo='sha256')),
  identical(x$driver_hash,digest(file=file.path(study,'run_dataset.R'),algo='sha256')))
 for(z in x$records) {
  stopifnot(z$converged,all(is.finite(z$positions)),all(is.finite(z$objective)),
   max(abs(rowSums(z$weights)-1))<1e-10,z$initial_ARI==1,z$ARI==1,!any(z$fallback))
  if(z$method=='one_step_k')for(g in 1:2) {
   tab<-subset(x$candidates,regime==z$regime & group==g & eligible)
   stopifnot(z$selected_k[g]==tab$k[which.max(tab$score)])
  }
 }
}
for(seed in unique(manifest$seed)) {
 ds<-lapply(manifest$id[manifest$seed==seed],function(id)readRDS(file.path(study,'data',paste0(id,'.rds'))))
 stopifnot(length(unique(vapply(ds,function(d)d$signal_hash,'')))==1,
  length(unique(vapply(ds,function(d)d$truth_hash,'')))==1,
  identical(ds[[1]]$unit_noise,ds[[2]]$unit_noise),identical(ds[[1]]$unit_noise,ds[[3]]$unit_noise))
}
# Recompute local one-sweep scores for a successful repair and a regression.
for(id in c('rich_S16_r03','rich_S4_r06')) {
 d<-readRDS(file.path(study,'data',paste0(id,'.rds')))
 x<-readRDS(file.path(study,'results',paste0(id,'.rds')))
 for(g in 1:2)for(k in design$k_candidates) {
  tab<-x$candidates[x$candidates$regime=='truth' & x$candidates$group==g & x$candidates$k==k,]
  X<-d$X[,d$truth==g,drop=FALSE];iso<-suppressWarnings(MPCurver:::isomap_ordering(X,k=k))
  if(tab$eligible) {
   one<-local_fit(X,iso$t,1L)
   stopifnot(one$iter==1,length(one$elbo_trace)==2,abs(tail(one$elbo_trace,1)-tab$score)<1e-7)
  }
 }
}
cat('PASS: 36 regenerated inputs, 12 paired realizations, 144 converged fits, all selections and 20 recomputed one-sweep candidates.\n')
