# Compare the paired P=12 and P=50 experiments without selecting endpoints.
source('experiments/m1_bspline_comparison_p50/common.R')
small <- read.csv(file.path(design$parent_study,'metrics.csv'))
large <- read.csv(file.path(study,'metrics.csv'))
paired <- merge(small,large,by=c('id','replication','noise','method'),suffixes=c('_p12','_p50'))
stopifnot(nrow(paired)==300L,!anyDuplicated(paired[c('id','method')]))
paired$recovery_change <- paired$rho_p50-paired$rho_p12
paired$runtime_ratio <- ifelse(paired$method=='PCA',NA_real_,paired$elapsed_seconds_p50/paired$elapsed_seconds_p12)
write.csv(paired,file.path(study,'dimension_pairs.csv'),row.names=FALSE)
methods <- c('PCA','MPCurve','PCurve','GPLVM','BayesianGPLVM')
labels <- c('Initial PCA','MPCurve','Principal\ncurve','GPLVM','Bayesian\nGPLVM')
colors <- c('#999999','#0072B2','#D55E00','#009E73','#CC79A7')
comparison <- do.call(rbind,lapply(c('low','high'),function(noise)do.call(rbind,lapply(methods,function(method){
  x <- paired[paired$noise==noise & paired$method==method,]
  data.frame(noise=noise,method=method,n=nrow(x),
    median_rho_p12=median(x$rho_p12),median_rho_p50=median(x$rho_p50),
    mean_change=mean(x$recovery_change),median_change=median(x$recovery_change),
    mcse_change=sd(x$recovery_change)/sqrt(nrow(x)),
    gains_over05=sum(x$recovery_change>.05),losses_over05=sum(x$recovery_change< -.05),
    median_seconds_p12=median(x$elapsed_seconds_p12),median_seconds_p50=median(x$elapsed_seconds_p50),
    median_paired_runtime_ratio=if(method=='PCA')NA_real_ else median(x$runtime_ratio),
    converged_p50=sum(x$converged_p50))
}))))
write.csv(comparison,file.path(study,'dimension_summary.csv'),row.names=FALSE)
save_plot <- function(name,draw,width=12,height=5.8) {
  png(file.path(study,'figures',paste0(name,'.png')),width=width,height=height,units='in',res=160)
  draw();dev.off()
  pdf(file.path(study,'figures',paste0(name,'.pdf')),width=width,height=height,useDingbats=FALSE)
  draw();dev.off()
}
save_plot('feature_count_effect',function(){
  par(mfrow=c(1,2),mar=c(4.5,4.8,3.3,1),oma=c(1.5,0,0,0),las=1)
  for(noise in c('low','high')) {
    d <- paired[paired$noise==noise,]
    boxplot(lapply(methods,function(method)d$recovery_change[d$method==method]),
      names=labels,col=adjustcolor(colors,.28),border=colors,ylim=c(-1,1),outline=FALSE,
      ylab='Absolute Spearman change: P = 50 minus P = 12',cex.axis=.85,
      main=if(noise=='low')'Low noise (SNR 16)' else 'High noise (SNR 1)')
    abline(h=0,lty=2,col='#555555')
    for(j in seq_along(methods)) {
      x <- d[d$method==methods[j],]
      points(j+(x$replication-15.5)/150,x$recovery_change,
             pch=ifelse(x$converged_p12 & x$converged_p50,16,4),col=colors[j],cex=.75)
    }
  }
  mtext('30 paired feature extensions | positive values favor P = 50 | crosses: either fit did not meet its stopping rule',side=1,outer=TRUE,cex=.8)
})
save_plot('runtime_comparison',function(){
  par(mfrow=c(1,2),mar=c(5,4.8,3.2,1),oma=c(1,0,0,0),las=1)
  for(noise in c('low','high')) {
    d <- comparison[comparison$noise==noise & comparison$method!='PCA',]
    heights <- rbind(d$median_seconds_p12,d$median_seconds_p50)
    bars <- barplot(heights,beside=TRUE,names.arg=labels[-1],col=c('#B8C8D9','#0072B2'),
      ylim=c(0,max(heights)*1.22),ylab='Median fitting time (seconds)',cex.names=.85,
      main=if(noise=='low')'Low noise (SNR 16)' else 'High noise (SNR 1)')
    text(bars,heights,labels=sprintf('%.2f',heights),pos=3,cex=.8)
    legend('topleft',c('P = 12','P = 50'),fill=c('#B8C8D9','#0072B2'),bty='n',cex=.85)
  }
  mtext('All 30 endpoints per condition; descriptive timings under concurrent single-threaded fits.',side=1,outer=TRUE,cex=.8)
})
print(comparison,row.names=FALSE)
# Optimization effort helps interpret changes in observed fitting time.
evaluations <- list()
for (P in c(12L,50L)) {
  folder <- if(P==12L)design$parent_study else study
  for (noise in c('low','high')) for(method in c('GPLVM','BayesianGPLVM')) {
    files <- file.path(folder,'results',sprintf('%s_r%02d_%s.json',noise,seq_len(design$replications),method))
    counts <- vapply(files,function(path)sum(jsonlite::fromJSON(path)$blocks$evaluations),numeric(1))
    evaluations[[length(evaluations)+1L]] <- data.frame(P=P,noise=noise,method=method,
      median_evaluations=median(counts))
  }
}
write.csv(do.call(rbind,evaluations),file.path(study,'gp_evaluation_summary.csv'),row.names=FALSE)
