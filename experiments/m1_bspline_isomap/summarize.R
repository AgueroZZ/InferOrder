source('experiments/m1_bspline_isomap/common.R')
manifest <- read.csv(file.path(study,'manifest.csv'))
methods <- c('Isomap','MPCurve','PCurve','GPLVM','BayesianGPLVM')
labels <- c('Initial\nIsomap','MPCurve','Principal\ncurve','GPLVM','Bayesian\nGPLVM')
colors <- c('#999999','#0072B2','#D55E00','#009E73','#CC79A7')
rows <- list(); positions <- list(); paired <- list()
for(i in seq_len(nrow(manifest))) {
  m <- manifest[i,];d <- readRDS(file.path(study,'inputs',paste0(m$id,'.rds')))
  parent <- readRDS(file.path(m$parent_study,'inputs',paste0(m$parent_id,'.rds')))
  stopifnot(identical(d$X,parent$X),identical(d$truth,parent$truth),
    identical(m$input_sha256,digest::digest(file=file.path(study,'inputs',paste0(m$id,'_X.csv')),algo='sha256')),
    identical(m$input_sha256,digest::digest(file=file.path(m$parent_study,'inputs',paste0(m$parent_id,'_X.csv')),algo='sha256')),
    d$embedding$n_components==1L,length(d$embedding$keep_idx)==design$N,
    all(is.finite(d$embedding$geodesic_to_landmark)),
    abs(mean(d$initial))<1e-12,abs(mean(d$initial^2)-1)<1e-12)
  # Independent classical-MDS audit of the complete saved geodesic matrix.
  D <- d$embedding$geodesic_to_landmark
  stopifnot(identical(d$embedding$landmark_idx,seq_len(design$N)),max(abs(D-t(D)))<1e-8)
  mds <- as.numeric(cmdscale(as.dist(D),k=1))
  stopifnot(recovery(mds,d$initial)>1-1e-10)
  pc_initial <- princurve::principal_curve(d$X,start=d$start_curve,maxit=0)$lambda
  for(method in methods) {
    if(method=='Isomap') {
      f <- list(initial=d$initial,final=d$initial,converged=TRUE,status='baseline',
        elapsed_seconds=0,warnings=character(),input_sha256=m$input_sha256)
    } else if(method %in% c('MPCurve','PCurve')) {
      f <- readRDS(file.path(study,'results',paste0(m$id,'_',method,'.rds')))
      stopifnot(identical(f$script_sha256,digest::digest(file=file.path(study,'run_r_methods.R'),algo='sha256')))
    } else {
      f <- jsonlite::fromJSON(file.path(study,'results',paste0(m$id,'_',method,'.json')))
      stopifnot(identical(f$script_sha256,digest::digest(file=file.path(study,'run_gpy.py'),algo='sha256')))
    }
    stopifnot(identical(f$input_sha256,m$input_sha256),length(f$final)==design$N,all(is.finite(f$final)))
    start_agreement <- recovery(d$initial,f$initial)
    if(method %in% c('GPLVM','BayesianGPLVM')) stopifnot(max(abs(f$initial-d$initial))<1e-10)
    if(method=='MPCurve')stopifnot(start_agreement>0.999,length(unique(f$initial))==design$K,
      all(diff(f$initial[order(d$initial)])>=0))
    if(method=='PCurve')stopifnot(max(abs(f$initial-pc_initial))<1e-10)
    rho <- recovery(d$truth,f$final);tau <- recovery(d$truth,f$final,'kendall')
    if(!is.null(f$rho))stopifnot(abs(rho-f$rho)<1e-10)
    row <- data.frame(id=m$id,parent_id=m$parent_id,P=m$P,replication=m$replication,noise=m$noise,
      method=method,rho=rho,tau=tau,initial_rho=recovery(d$truth,d$initial),
      actual_initial_rho=recovery(d$truth,f$initial),start_agreement=start_agreement,
      gain=rho-recovery(d$truth,d$initial),gain_from_actual=rho-recovery(d$truth,f$initial),
      converged=isTRUE(f$converged),status=f$status,fit_seconds=f$elapsed_seconds,
      embedding_seconds=m$embedding_seconds,
      curve_seconds=if(method=='PCurve')m$curve_seconds else 0,
      total_seconds=f$elapsed_seconds+m$embedding_seconds+if(method=='PCurve')m$curve_seconds else 0,
      warning_count=length(f$warnings))
    rows[[length(rows)+1L]] <- row
    positions[[length(positions)+1L]] <- data.frame(id=m$id,method=method,sample=seq_len(design$N),
      truth=d$truth,upstream=d$initial,actual_initial=f$initial,final=f$final)
  }
  pca <- read.csv(file.path(m$parent_study,'metrics.csv'))
  pca <- pca[pca$id==m$parent_id,]
  pca$method[pca$method=='PCA'] <- 'Isomap'
  current <- do.call(rbind,tail(rows,5))
  comparison <- merge(current,pca[,c('method','rho','converged','elapsed_seconds')],by='method',suffixes=c('_isomap','_pca'))
  comparison$delta <- comparison$rho_isomap-comparison$rho_pca
  paired[[length(paired)+1L]] <- comparison
}
metrics <- do.call(rbind,rows);comparison <- do.call(rbind,paired)
stopifnot(nrow(metrics)==600L,nrow(comparison)==600L,!anyDuplicated(metrics[c('id','method')]))
write.csv(metrics,file.path(study,'metrics.csv'),row.names=FALSE)
write.csv(do.call(rbind,positions),file.path(study,'positions.csv'),row.names=FALSE)
write.csv(comparison,file.path(study,'initialization_pairs.csv'),row.names=FALSE)
summary <- do.call(rbind,lapply(design$P,function(P)do.call(rbind,lapply(c('low','high'),function(noise)
  do.call(rbind,lapply(methods,function(method){
    x <- metrics[metrics$P==P & metrics$noise==noise & metrics$method==method,]
    z <- comparison[comparison$P==P & comparison$noise==noise & comparison$method==method,]
    data.frame(P=P,noise=noise,method=method,n=nrow(x),median_rho=median(x$rho),mean_rho=mean(x$rho),
      median_pca_rho=median(z$rho_pca),mean_change=mean(z$delta),median_change=median(z$delta),
      mcse_change=sd(z$delta)/sqrt(nrow(z)),median_tau=median(x$tau),
      median_actual_initial_rho=median(x$actual_initial_rho),median_start_agreement=median(x$start_agreement),
      median_gain=median(x$gain),gains_over05=sum(x$gain>.05),losses_over05=sum(x$gain< -.05),
      converged=sum(x$converged),recovery95=sum(x$rho>=.95),
      median_fit_seconds=median(x$fit_seconds),median_total_seconds=median(x$total_seconds),
      median_pca_fit_seconds=median(z$elapsed_seconds))
  }))))))
write.csv(summary,file.path(study,'summary.csv'),row.names=FALSE)
# Compare fitted methods within the same Isomap dataset.
method_pairs <- merge(metrics[metrics$method=='MPCurve',c('id','rho')],metrics,by='id',suffixes=c('_mpcurve',''))
method_pairs <- method_pairs[method_pairs$method %in% c('PCurve','GPLVM','BayesianGPLVM'),]
method_pairs$delta <- method_pairs$rho-method_pairs$rho_mpcurve
method_summary <- do.call(rbind,lapply(split(method_pairs,list(method_pairs$P,method_pairs$noise,method_pairs$method),drop=TRUE),function(x)
  data.frame(P=x$P[1],noise=x$noise[1],method=x$method[1],mean_delta=mean(x$delta),median_delta=median(x$delta),
    mcse=sd(x$delta)/sqrt(nrow(x)),wins=sum(x$delta>0))))
write.csv(method_summary,file.path(study,'method_comparison.csv'),row.names=FALSE)
writeLines(c('Verified all 120 input matrices and hashes against the PCA studies; all Isomap graphs connected;',
  'independent classical MDS ranks; shared GP coordinates; monotone MPCurve bins; PCurve initial projection replay;',
  'all 600 baseline/method records finite; metrics and driver hashes verified; all endpoints retained.',
  sprintf('Converged method fits: %d/480.',sum(metrics$converged[metrics$method!='Isomap']))),file.path(study,'verification.txt'))

save_plot <- function(name,draw,width=12,height=10) {
  png(file.path(study,'figures',paste0(name,'.png')),width=width,height=height,units='in',res=160)
  draw();dev.off()
  pdf(file.path(study,'figures',paste0(name,'.pdf')),width=width,height=height,useDingbats=FALSE)
  draw();dev.off()
}
save_plot('ordering_recovery',function(){
  par(mfrow=c(2,2),mar=c(4.1,4.8,2.9,1),oma=c(2.5,0,0,0),las=1)
  for(P in design$P) for(noise in c('low','high')) {
    d <- metrics[metrics$P==P & metrics$noise==noise,]
    boxplot(lapply(methods,function(method)d$rho[d$method==method]),names=labels,
      col=adjustcolor(colors,.28),border=colors,outline=FALSE,ylim=c(0,1.055),
      ylab='Absolute Spearman correlation',cex.axis=.83,
      main=sprintf('P = %d | %s noise (SNR %d)',P,noise,if(noise=='low')16 else 1))
    abline(h=.95,col='#bbbbbb',lty=2)
    for(j in seq_along(methods)) {
      x <- d[d$method==methods[j],]
      points(j+(x$replication-15.5)/150,x$rho,pch=ifelse(!x$converged,4,ifelse(x$method=='PCurve' & x$warning_count>0,17,16)),col=colors[j],cex=.65)
      text(j,1.035,sprintf('%.3f',median(x$rho)),cex=.8)
    }
  }
  mtext('Isomap k = 15; N = 200; 30 replicates per panel | numbers: medians',side=1,outer=TRUE,cex=.8)
  mtext('Dots: converged; crosses: stopping rule unmet; triangle: numerical smoothing warning',side=1,outer=TRUE,line=1.1,cex=.75)
})
save_plot('initialization_effect',function(){
  par(mfrow=c(2,2),mar=c(4.1,5,2.9,1),oma=c(2.5,0,0,0),las=1)
  for(P in design$P) for(noise in c('low','high')) {
    d <- comparison[comparison$P==P & comparison$noise==noise,]
    boxplot(lapply(methods,function(method)d$delta[d$method==method]),
      names=c('Initial\nembedding',labels[-1]),col=adjustcolor(colors,.28),border=colors,
      outline=FALSE,ylim=c(-1,1),ylab='Recovery change: Isomap start minus PCA start',cex.axis=.83,
      main=sprintf('P = %d | %s noise',P,noise))
    abline(h=0,lty=2,col='#555555')
    for(j in seq_along(methods)) {
      x <- d[d$method==methods[j],]
      points(j+(x$replication-15.5)/150,x$delta,pch=ifelse(!(x$converged_isomap & x$converged_pca),4,ifelse(x$method=='PCurve' & x$warning_count>0,17,16)),col=colors[j],cex=.65)
    }
  }
  mtext('Same observations and settings; only initialization changes',side=1,outer=TRUE,cex=.8)
  mtext('Crosses: either stopping rule unmet; triangle: numerical smoothing warning',side=1,outer=TRUE,line=1.1,cex=.75)
})
print(summary,row.names=FALSE)
print(method_summary,row.names=FALSE)
