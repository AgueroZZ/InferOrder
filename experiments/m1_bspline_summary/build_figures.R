# Three-figure overview of the existing replicated single-ordering experiments.
# This script reads saved inputs and endpoints; it does not fit any models.
source('experiments/m1_bspline_comparison/common.R')
summary_dir <- 'experiments/m1_bspline_summary'
methods <- c('Initial','MPCurve','PCurve','GPLVM','BayesianGPLVM')
colors <- c('#969696','#0072B2','#D55E00','#009E73','#CC79A7')

pca12 <- read.csv('experiments/m1_bspline_comparison/metrics.csv')
pca50 <- read.csv('experiments/m1_bspline_comparison_p50/metrics.csv')
pca12$P <- 12L;pca50$P <- 50L
pca <- rbind(pca12,pca50)
pca$initialization <- 'PCA'
pca$method[pca$method=='PCA'] <- 'Initial'
isomap <- read.csv('experiments/m1_bspline_isomap/metrics.csv')
isomap$initialization <- 'Isomap'
isomap$method[isomap$method=='Isomap'] <- 'Initial'
columns <- c('P','noise','replication','method','initialization','rho','converged','warning_count')
plot_data <- rbind(pca[,columns],isomap[,columns])
stopifnot(nrow(plot_data)==1200L,all(is.finite(plot_data$rho)),
          all(table(plot_data$P,plot_data$noise,plot_data$method,plot_data$initialization)==30L),
          !anyDuplicated(plot_data[c('P','noise','replication','method','initialization')]))
plot_data$numerical_warning <- plot_data$method=='PCurve' & plot_data$warning_count>0
stopifnot(sum(plot_data$numerical_warning)==1L)
write.csv(plot_data,file.path(summary_dir,'ordering_plot_data.csv'),row.names=FALSE)

# Fixed example selected by index and shared across noise panels. Display raw
# generated values (before observation-column centering) against generating t.
sim <- readRDS('experiments/m1_bspline_comparison/inputs/rep01.rds')
features <- 1:3
noise_levels <- c(low=.25,high=1)
observations <- do.call(rbind,lapply(names(noise_levels),function(noise)
  do.call(rbind,lapply(features,function(j)data.frame(noise=noise,feature=j,
    sample=seq_along(sim$truth),t=sim$truth,
    value=sim$signal[,j]+noise_levels[[noise]]*sim$noise[,j])))))
trajectories <- do.call(rbind,lapply(features,function(j)
  data.frame(feature=j,t=sim$grid,value=sim$dense_signal[,j])))
write.csv(observations,file.path(summary_dir,'example_observations.csv'),row.names=FALSE)
write.csv(trajectories,file.path(summary_dir,'example_trajectories.csv'),row.names=FALSE)
write.csv(data.frame(P=12,replication=1,seed=sim$seed,features='1,2,3',
  selection='Fixed first replicate and first three features; shared signals and errors across noise levels'),
  file.path(summary_dir,'example_selection.csv'),row.names=FALSE)

save_figure <- function(name,draw,width,height) {
  png(file.path(summary_dir,'figures',paste0(name,'.png')),width=width,height=height,units='in',res=180)
  draw();dev.off()
  pdf(file.path(summary_dir,'figures',paste0(name,'.pdf')),width=width,height=height,useDingbats=FALSE)
  draw();dev.off()
}
save_figure('figure1_trajectories',function(){
  par(mfrow=c(1,2),mar=c(4.1,4.4,3.1,.8),oma=c(.1,0,0,0),las=1)
  feature_colors <- colors[2:4]
  limits <- range(observations$value,trajectories$value)
  for(noise in names(noise_levels)) {
    plot(NA,xlim=c(0,1),ylim=limits,xlab='True sample position',ylab='Feature value',
      main=if(noise=='low')'Simpler: low noise (SD 0.25)' else 'Challenging: high noise (SD 1)',
      cex.lab=1.05,cex.axis=.9)
    for(j in features) {
      d <- observations[observations$noise==noise & observations$feature==j,]
      points(d$t,d$value,pch=16,cex=.48,col=adjustcolor(feature_colors[j],.36))
    }
    for(j in features)lines(sim$grid,sim$dense_signal[,j],col=feature_colors[j],lwd=2.1)
    legend('topright',paste('Feature',features),col=feature_colors,lty=1,lwd=2,pch=16,
      bty='n',cex=.82,inset=.015)
  }
},width=11,height=4.4)

plot_recovery <- function(initialization) {
  par(mfrow=c(2,2),mar=c(4.1,4.5,3.1,.8),oma=c(2.1,0,0,0),las=1)
  labels <- c(paste('Initial',initialization,sep='\n'),'MPCurve','Principal\ncurve','GPLVM','Bayesian\nGPLVM')
  for(P in c(12L,50L))for(noise in c('low','high')) {
    panel <- plot_data[plot_data$P==P & plot_data$noise==noise & plot_data$initialization==initialization,]
    boxplot(lapply(methods,function(method)panel$rho[panel$method==method]),
      names=labels,col=adjustcolor(colors,.25),border=colors,outline=FALSE,
      ylim=c(0,1.06),cex.axis=.85,ylab='Absolute Spearman correlation',
      main=sprintf('P = %d | %s noise (SNR %d)',P,if(noise=='low')'Low' else 'High',if(noise=='low')16 else 1))
    abline(h=.95,col='#C8C8C8',lty=2)
    for(j in seq_along(methods)) {
      d <- panel[panel$method==methods[j],]
      point_shape <- ifelse(!d$converged,4,ifelse(d$numerical_warning,17,16))
      points(j+(d$replication-15.5)/150,d$rho,pch=point_shape,cex=.65,col=colors[j])
      center <- median(d$rho)
      text(j,1.04,if(center>=.995)sprintf('%.4f',center) else sprintf('%.3f',center),cex=.8)
    }
  }
  mtext('30 replicates per panel | numbers: medians | dashed line: 0.95 recovery',side=1,outer=TRUE,cex=.8)
  mtext('Dots: converged or initial baseline; crosses: stopping rule unmet; triangle: numerical warning',
        side=1,outer=TRUE,line=1,cex=.73)
}
save_figure('figure2_pca',function()plot_recovery('PCA'),width=11,height=8.8)
save_figure('figure3_isomap',function()plot_recovery('Isomap'),width=11,height=8.8)
input_files <- c('experiments/m1_bspline_comparison/inputs/rep01.rds',
 'experiments/m1_bspline_comparison/metrics.csv','experiments/m1_bspline_comparison_p50/metrics.csv',
 'experiments/m1_bspline_isomap/metrics.csv')
write.csv(data.frame(path=input_files,sha256=vapply(input_files,function(path)
  digest::digest(file=path,algo='sha256'),character(1))),file.path(summary_dir,'input_hashes.csv'),row.names=FALSE)
writeLines(c('All 1,200 displayed baseline/model endpoints match the saved source metrics.',
 'Thirty replicates per method, feature count, noise level and initialization; no endpoints filtered.',
 'Figure 1 uses fixed replicate 1, P=12, features 1-3, identical signal and standard errors at both noise levels.'),
 file.path(summary_dir,'verification.txt'))
