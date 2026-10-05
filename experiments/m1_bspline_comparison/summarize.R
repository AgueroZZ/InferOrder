source('experiments/m1_bspline_comparison/common.R')
manifest <- read.csv(file.path(study,'manifest.csv'))
methods <- c('PCA','MPCurve','PCurve','GPLVM','BayesianGPLVM')
labels <- c('Initial PCA','MPCurve','Principal\ncurve','GPLVM','Bayesian\nGPLVM')
colors <- c('#999999','#0072B2','#D55E00','#009E73','#CC79A7')
rows <- list(); positions <- list(); audits <- list()
for (i in seq_len(nrow(manifest))) {
  m <- manifest[i,]; d <- readRDS(file.path(study,'inputs',paste0(m$id,'.rds')))
  sim <- readRDS(sprintf('%s/inputs/rep%02d.rds',study,m$replication))
  expected <- scale(sim$signal+design$noise_sd[[m$noise]]*sim$noise,center=TRUE,scale=FALSE)
  stopifnot(identical(unname(expected),unname(d$X)),
            identical(d$truth,sim$truth),
            abs(mean(d$initial))<1e-12,abs(mean(d$initial^2)-1)<1e-12,
            identical(m$input_sha256,digest::digest(file=file.path(study,'inputs',paste0(m$id,'_X.csv')),algo='sha256')))
  for (method in methods) {
    if (method=='PCA') {
      f <- list(final=d$initial,initial=d$initial,converged=TRUE,status='baseline',
                elapsed_seconds=0,warnings=character(),input_sha256=m$input_sha256)
    } else if (method %in% c('MPCurve','PCurve')) {
      f <- readRDS(file.path(study,'results',paste0(m$id,'_',method,'.rds')))
    } else {
      f <- jsonlite::fromJSON(file.path(study,'results',paste0(m$id,'_',method,'.json')))
    }
    stopifnot(identical(f$input_sha256,m$input_sha256))
    valid <- length(f$final)==design$N && all(is.finite(f$final))
    # All returned endpoints remain in the comparison, including budget limits.
    rho <- if(valid) recovery(d$truth,f$final) else NA_real_
    tau <- if(valid) recovery(d$truth,f$final,'kendall') else NA_real_
    start_agreement <- if(length(f$initial)==design$N) recovery(d$initial,f$initial) else NA_real_
    if(method %in% c('PCurve','GPLVM','BayesianGPLVM')) stopifnot(start_agreement>1-1e-10)
    if(method=='MPCurve') stopifnot(start_agreement>0.999,
      length(unique(f$initial))==design$K,all(diff(f$initial[order(d$initial)])>=0))
    if(!is.null(f$rho) && is.finite(f$rho)) stopifnot(abs(rho-f$rho)<1e-10)
    row <- data.frame(id=m$id,replication=m$replication,noise=m$noise,method=method,
      rho=rho,tau=tau,initial_rho=recovery(d$truth,d$initial),
      actual_initial_rho=if(is.finite(start_agreement))recovery(d$truth,f$initial) else NA_real_,
      start_agreement=start_agreement,gain=rho-recovery(d$truth,d$initial),
      converged=isTRUE(f$converged),status=f$status,valid=valid,
      elapsed_seconds=f$elapsed_seconds,warning_count=length(f$warnings))
    rows[[length(rows)+1L]] <- row
    if(valid) positions[[length(positions)+1L]] <- data.frame(id=m$id,method=method,
       sample=seq_len(design$N),truth=d$truth,initial=f$initial,final=f$final)
  }
}
results <- do.call(rbind,rows)
stopifnot(nrow(results)==300L,!anyDuplicated(results[c('id','method')]))
write.csv(results,file.path(study,'metrics.csv'),row.names=FALSE)
write.csv(do.call(rbind,positions),file.path(study,'positions.csv'),row.names=FALSE)
summary <- do.call(rbind,lapply(c('low','high'),function(noise)do.call(rbind,lapply(methods,function(method){
  x <- results[results$noise==noise & results$method==method,]
  data.frame(noise=noise,method=method,n=nrow(x),valid=sum(x$valid),converged=sum(x$converged),
             median_rho=median(x$rho,na.rm=TRUE),mean_rho=mean(x$rho,na.rm=TRUE),
             q25=unname(quantile(x$rho,.25,na.rm=TRUE)),q75=unname(quantile(x$rho,.75,na.rm=TRUE)),
             median_tau=median(x$tau,na.rm=TRUE),median_gain=median(x$gain,na.rm=TRUE),
             improved=sum(x$gain>0.05,na.rm=TRUE),worsened=sum(x$gain< -0.05,na.rm=TRUE),
             recovery95=sum(x$rho>=0.95,na.rm=TRUE),median_seconds=median(x$elapsed_seconds))
}))))
write.csv(summary,file.path(study,'summary.csv'),row.names=FALSE)
paired <- merge(results[results$method=='MPCurve',c('id','rho')],results,by='id',suffixes=c('_mpcurve',''))
paired$delta <- paired$rho-paired$rho_mpcurve
paired <- paired[paired$method %in% c('PCurve','GPLVM','BayesianGPLVM'),]
write.csv(paired,file.path(study,'paired_differences.csv'),row.names=FALSE)
paired_summary <- do.call(rbind,lapply(split(paired,list(paired$noise,paired$method),drop=TRUE),function(x)
  data.frame(noise=x$noise[1],method=x$method[1],n=sum(is.finite(x$delta)),
             mean_delta=mean(x$delta,na.rm=TRUE),median_delta=median(x$delta,na.rm=TRUE),
             mcse=sd(x$delta,na.rm=TRUE)/sqrt(sum(is.finite(x$delta))),
             wins=sum(x$delta>0,na.rm=TRUE))))
write.csv(paired_summary,file.path(study,'paired_summary.csv'),row.names=FALSE)

save_plot <- function(name, draw, width=12,height=5.8) {
  png(file.path(study,'figures',paste0(name,'.png')),width=width,height=height,units='in',res=160)
  draw();dev.off()
  pdf(file.path(study,'figures',paste0(name,'.pdf')),width=width,height=height,useDingbats=FALSE)
  draw();dev.off()
}
plot_recovery <- function(metric='rho') {
  par(mfrow=c(1,2),mar=c(4.5,4.8,3.3,1),oma=c(1.5,0,0,0),las=1)
  for(noise in c('low','high')) {
    d <- results[results$noise==noise,]
    values <- lapply(methods,function(method)d[d$method==method,metric])
    boxplot(values,names=labels,col=adjustcolor(colors,.28),border=colors,
      ylim=c(0,1.04),outline=FALSE,ylab=if(metric=='rho')'Absolute Spearman correlation' else 'Absolute Kendall correlation',
      main=if(noise=='low')'Low noise: SD 0.25 (SNR 16)' else 'High noise: SD 1 (SNR 1)',cex.axis=.85)
    abline(h=c(.5,.95),col='#dddddd',lty=c(3,2))
    for(j in seq_along(methods)) {
      x <- d[d$method==methods[j],];offset <- (x$replication-15.5)/150
      points(j+offset,x[[metric]],pch=ifelse(x$converged,16,4),col=colors[j],cex=.72)
      text(j,1.025,sprintf('%.3f',median(x[[metric]],na.rm=TRUE)),cex=.8)
    }
  }
  mtext('30 paired replicates per noise level | dots: converged; crosses: stopping rule unmet | numbers: medians',side=1,outer=TRUE,cex=.8)
}
save_plot('ordering_recovery',function()plot_recovery('rho'))
save_plot('kendall_recovery',function()plot_recovery('tau'))
save_plot('paired_differences',function(){
  par(mfrow=c(1,2),mar=c(4.5,5,3.2,1),oma=c(1.5,0,0,0),las=1)
  comparison <- c('PCurve','GPLVM','BayesianGPLVM')
  for(noise in c('low','high')) {
    d <- paired[paired$noise==noise,]
    boxplot(lapply(comparison,function(method)d$delta[d$method==method]),
      names=c('Principal curve','GPLVM','Bayesian GPLVM'),col=adjustcolor(colors[3:5],.28),
      border=colors[3:5],outline=FALSE,ylim=c(-1,1),
      main=if(noise=='low')'Low noise (SNR 16)' else 'High noise (SNR 1)',
      ylab='Recovery difference: comparator minus MPCurve',cex.axis=.85)
    abline(h=0,lty=2,col='#555555')
    for(j in seq_along(comparison)) {
      x <- d[d$method==comparison[j],]
      points(j+(x$replication-15.5)/150,x$delta,pch=ifelse(x$converged,16,4),col=colors[j+2],cex=.8)
    }
  }
  mtext('Each point compares methods on the same data; positive values favor the comparator.',side=1,outer=TRUE,cex=.85)
})
save_plot('example_trajectories',function(){
  sim <- readRDS(file.path(study,'inputs/rep01.rds'))
  par(mfrow=c(1,3),mar=c(4.4,4.4,3,1),las=1)
  matplot(sim$grid,sim$dense_signal,type='l',lty=1,col=hcl.colors(design$P,'Dark 3'),
    xlab='True latent position',ylab='Noiseless feature value',main='Replicate 1: all 12 trajectories')
  for(noise in c('low','high')) {
    d <- readRDS(file.path(study,'inputs',paste0(noise,'_r01.rds')))
    # Use the saved shared initialization, preserving score spacing and sign.
    plot(d$truth,d$initial,pch=16,col='#0072B2',cex=.7,
         xlab='True latent position',ylab='PC1 score (unit variance)',
         main=if(noise=='low')'Low noise (SNR 16)' else 'High noise (SNR 1)')
    legend('topleft',sprintf('Absolute Spearman = %.3f',recovery(d$truth,d$initial)),
           bty='n',cex=.85)
  }
},width=13,height=4.6)
writeLines(c('Verified: 60 paired inputs; 300 method/baseline records; input SHA256 hashes;',
  'exact signal/noise reconstruction; shared PCA ranks; MPCurve monotone 50-bin quantization;',
  'R/Python recovery metric agreement; all endpoints retained regardless of convergence.',
  sprintf('Finite endpoints: %d/300. Converged method fits: %d/240.',sum(results$valid),sum(results$converged[results$method!='PCA']))),
  file.path(study,'verification.txt'))
print(summary,row.names=FALSE)
print(paired_summary,row.names=FALSE)
