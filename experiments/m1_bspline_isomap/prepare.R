source('experiments/m1_bspline_isomap/common.R')
manifest <- list()
for(P in design$P) {
  parent <- if(P==12L)'experiments/m1_bspline_comparison' else 'experiments/m1_bspline_comparison_p50'
  parent_manifest <- read.csv(file.path(parent,'manifest.csv'))
  for(i in seq_len(nrow(parent_manifest))) {
    m <- parent_manifest[i,]
    d <- readRDS(file.path(parent,'inputs',paste0(m$id,'.rds')))
    id <- sprintf('p%d_%s',P,m$id)
    started <- proc.time()[['elapsed']]
    warnings <- character()
    embedding <- withCallingHandlers(MPCurver:::isomap_ordering(d$X,k=design$isomap_k,
      ndim=1L,landmark=design$N,seed=m$seed,keep='all'),
      warning=function(w){warnings<<-c(warnings,conditionMessage(w));invokeRestart('muffleWarning')})
    embedding_seconds <- proc.time()[['elapsed']]-started
    # Disconnected graphs are never silently restricted to their largest component.
    stopifnot(embedding$n_components==1L,length(embedding$keep_idx)==design$N,
              all(is.finite(embedding$t)),all(is.finite(embedding$geodesic_to_landmark)))
    initial <- embedding$t-mean(embedding$t)
    initial <- initial/sqrt(mean(initial^2))
    started <- proc.time()[['elapsed']]
    start_curve <- vapply(seq_len(P),function(j)
      predict(smooth.spline(initial,d$X[,j],df=5),x=sort(initial))$y,numeric(design$N))
    curve_seconds <- proc.time()[['elapsed']]-started
    input <- list(id=id,parent_id=m$id,parent_study=parent,P=P,replication=m$replication,
      noise=m$noise,X=d$X,truth=d$truth,initial=initial,pca=d$initial,
      start_curve=start_curve,embedding=embedding,embedding_warnings=warnings,
      embedding_seconds=embedding_seconds,curve_seconds=curve_seconds)
    saveRDS(input,file.path(study,'inputs',paste0(id,'.rds')))
    stopifnot(file.copy(file.path(parent,'inputs',paste0(m$id,'_X.csv')),
      file.path(study,'inputs',paste0(id,'_X.csv')),overwrite=TRUE))
    write.csv(data.frame(truth=d$truth,isomap=initial),
      file.path(study,'inputs',paste0(id,'_positions.csv')),row.names=FALSE)
    manifest[[length(manifest)+1L]] <- data.frame(id=id,parent_id=m$id,parent_study=parent,
      P=P,replication=m$replication,noise=m$noise,noise_sd=m$noise_sd,variance_snr=m$variance_snr,
      seed=m$seed,input_sha256=m$input_sha256,n_components=embedding$n_components,
      embedding_seconds=embedding_seconds,curve_seconds=curve_seconds,
      pca_rho=recovery(d$truth,d$initial),isomap_rho=recovery(d$truth,initial))
  }
}
write.csv(do.call(rbind,manifest),file.path(study,'manifest.csv'),row.names=FALSE)
saveRDS(design,file.path(study,'design.rds'))
writeLines(capture.output(dput(design)),file.path(study,'design.txt'))
writeLines(capture.output(sessionInfo()),file.path(study,'R_session.txt'))
writeLines(deparse(MPCurver:::isomap_ordering),file.path(study,'isomap_source.R.txt'))
print(aggregate(isomap_rho~P+noise,do.call(rbind,manifest),median))
