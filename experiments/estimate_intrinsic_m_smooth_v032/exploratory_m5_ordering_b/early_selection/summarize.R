source('experiments/estimate_intrinsic_m_smooth_v032/auto_m_v034/common.R')
suppressPackageStartupMessages(library(ggplot2))
out <- file.path(study_dir,'exploratory_m5_ordering_b','early_selection')
paths <- list.files(file.path(out,'results'),pattern='[.]rds$',full.names=TRUE)
stopifnot(length(paths)==17L)
results <- lapply(paths,readRDS)
long <- do.call(rbind,lapply(results,function(x) data.frame(method=x$method,noise_scale=x$noise_scale,
 x$history,final_elbo=x$final_objective,final_rho=x$final_rho,final_sweeps=x$final_sweeps)))
write.csv(long,file.path(out,'early_scores.csv'),row.names=FALSE)
selected <- do.call(rbind,lapply(split(long,interaction(long$noise_scale,long$sweep,drop=TRUE)),function(x) {
 best <- which.max(x$standard_elbo); final_best <- max(x$final_elbo)
 data.frame(noise_scale=x$noise_scale[1],sweep=x$sweep[1],selected=x$method[best],
  eventual_elbo_loss=final_best-x$final_elbo[best],eventual_rho=x$final_rho[best],
  final_best=x$method[which.max(x$final_elbo)],
  final_best_early_rank=rank(-x$standard_elbo,ties.method='min')[which.max(x$final_elbo)],
  top_two_contains_final_best=rank(-x$standard_elbo,ties.method='min')[which.max(x$final_elbo)]<=2)
}))
selected <- selected[order(selected$noise_scale,selected$sweep),]
write.csv(selected,file.path(out,'selection_by_budget.csv'),row.names=FALSE)
comparison <- do.call(rbind,lapply(results,function(x) data.frame(noise_scale=x$noise_scale,method=x$method,
 first_elbo=x$history$standard_elbo[1],fifth_elbo=x$history$standard_elbo[5],
 tenth_elbo=x$history$standard_elbo[10],after_annealing_elbo=x$history$standard_elbo[25],
 final_elbo=x$final_objective,final_rho=x$final_rho)))
write.csv(comparison,file.path(out,'comparison.csv'),row.names=FALSE)
observed <- long[long$noise_scale==1,]
observed$loss <- ave(observed$standard_elbo,observed$sweep,FUN=function(x)max(x)-x)
observed$method <- factor(observed$method,levels=c('k=5','k=10','k=15','k=20','k=30','PCA'))
colors <- c('k=5'='#56B4E9','k=10'='#0072B2','k=15'='#D55E00','k=20'='#CC79A7','k=30'='#E69F00','PCA'='#444444')
p <- ggplot(observed,aes(sweep,loss,color=method))+geom_line(linewidth=.8)+geom_point(size=1.2)+
 scale_color_manual(values=colors)+scale_x_continuous(breaks=c(1,2,5,10,15,20,25))+
 theme_minimal(base_size=12)+theme(legend.position='bottom',plot.caption=element_text(hjust=0))+
 labs(title='The first-step ELBO selects a different Isomap neighborhood',
 subtitle='Original M=5, SNR=4 observations; smaller gap is better. Every state is scored with the ordinary T=1 ELBO.',
 x='Joint structural sweep (after the two subset initialization updates)',
 y='ELBO gap to the best candidate at the same sweep',color='Initialization',
 caption='The fitting path keeps its original T=5 to 1 annealing schedule; only evaluation uses T=1.\nThese are paired starts on one selected dataset, not an independent validation of a screening rule.')
for(ext in c('png','pdf')) ggsave(file.path(out,paste0('early_elbo_selection.',ext)),p,width=12,height=6,dpi=170,bg='white')
print(comparison[comparison$noise_scale==1,],row.names=FALSE)
print(selected[selected$sweep %in% c(1,2,3,5,10,25),],row.names=FALSE)
