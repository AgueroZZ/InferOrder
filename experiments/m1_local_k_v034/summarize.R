source('experiments/m2_local_k_v034/common.R')
library(ggplot2)
out<-'experiments/m1_local_k_v034'
paths<-list.files(file.path(out,'results'),full.names=TRUE,pattern='rds$')
stopifnot(length(paths)==72)
results<-lapply(paths,readRDS);rows<-list();all_candidates<-list()
for(x in results) {
 d<-readRDS(file.path(study,'data',paste0(x$row$id,'.rds')))
 stopifnot(identical(x$input_hash,digest(d$X[,d$truth==x$group,drop=FALSE],algo='sha256')),
 identical(x$source_hash,digest(file=file.path(out,'run.R'),algo='sha256')))
 t<-x$candidates;b<-t[t$k==15,];s<-t[t$k==x$selected_k,];v<-t[t$eligible,]
 stopifnot(b$eligible,all(v$converged),all(is.finite(v$score)),x$selected_k==t$k[which.max(t$score)])
 rows[[x$id]]<-data.frame(id=x$id,family=x$row$family,snr=x$row$snr,replicate=x$row$replicate,group=x$group,
  initial=b$initial_rho,default=b$final_rho,selected=s$final_rho,k=s$k,oracle=max(v$final_rho),
  broken=b$final_rho<.9,rescued=b$final_rho<.9 & s$final_rho>=.95,
  oracle_rescue=b$final_rho<.9 & max(v$final_rho)>=.95,
  initial_broken=b$initial_rho<.9,iteration_repair=b$initial_rho<.9 & b$final_rho>=.95,
  default_seconds=b$embedding+b$fitting,selected_seconds=sum(t$embedding+t$scoring)+s$fitting,
  screening_seconds=sum(t$embedding+t$scoring),scoring_seconds=sum(t$scoring))
 all_candidates[[x$id]]<-cbind(id=x$id,t)
}
pairs<-do.call(rbind,rows);candidates<-do.call(rbind,all_candidates)
pairs$delta<-pairs$selected-pairs$default
stats<-function(x)data.frame(n=nrow(x),initial_broken=sum(x$initial_broken),iteration_repair=sum(x$iteration_repair),broken=sum(x$broken),rescued=sum(x$rescued),oracle_rescue=sum(x$oracle_rescue),remaining=sum(x$selected<.9),gains=sum(x$delta>.05),drops=sum(x$delta< -.05),good_to_bad=sum(x$default>=.95 & x$selected<.9),mean_default=mean(x$default),mean_selected=mean(x$selected))
summary<-do.call(rbind,lapply(split(pairs,pairs$snr),function(x)cbind(snr=x$snr[1],stats(x))))
summary<-rbind(summary,cbind(snr='all',stats(pairs)))
for(n in c('pairs','candidates','summary'))write.csv(get(n),file.path(out,paste0(n,'.csv')),row.names=FALSE)
print(summary,row.names=FALSE)
print(data.frame(median_ratio=median(pairs$selected_seconds/pairs$default_seconds),median_extra=median(pairs$selected_seconds-pairs$default_seconds),median_screening=median(pairs$screening_seconds),median_scoring=median(pairs$scoring_seconds)))
cat('Eligible converged candidates:',sum(candidates$eligible),'/',nrow(candidates),'\n')
print(pairs[order(pairs$delta),c('id','default','selected','k','oracle','delta')][c(1:3,70:72),],row.names=FALSE)
print(subset(pairs,oracle_rescue & !rescued))
theme_set(theme_minimal(base_size=12))
p<-ggplot(pairs,aes(default,selected,color=family))+geom_abline(slope=1,intercept=0,color='gray60')+geom_vline(xintercept=.9,linetype=2,color='gray65')+geom_hline(yintercept=.95,linetype=2,color='gray65')+geom_point(size=2.2,alpha=.8)+facet_wrap(~snr,labeller=label_both)+coord_equal(xlim=c(0,1),ylim=c(0,1))+scale_color_manual(values=c(broad='#0072B2',rich='#D55E00'))+labs(title='Standalone M=1: which default failures can one-step screening rescue?',x='Converged recovery with default Isomap k=15',y='Final recovery: one-step selected k',color='Trajectory family',caption='Recovery is absolute Spearman correlation. Dashed thresholds: default failure <0.90; rescue >=0.95.\nAll 72 datasets retained; 24 underlying single-ordering realizations are paired across SNR levels.')+theme(legend.position='bottom')
ggsave(file.path(out,'recovery.png'),p,width=12,height=5,dpi=170,bg='white')
examples<-pairs[c(which.max(pairs$delta),which.min(pairs$delta)),];points<-list()
for(id in examples$id) {
 x<-results[[match(id,vapply(results,function(x)x$id,''))]]
 for(method in c('Default k=15','One-step selected k')) {
 k<-if(method=='Default k=15')15 else x$selected_k
 for(stage in c('initial','final')) {
 pos<-x$positions[[as.character(k)]][[stage]]
 if(cor(x$truth,pos,method='spearman')<0)pos<- -pos
 pos<-(pos-min(pos))/diff(range(pos))
 points[[length(points)+1]]<-data.frame(id=id,method=if(method=='Default k=15')method else paste0(method,' (k=',k,')'),stage=stage,truth=x$truth,position=pos)
 }
 }
}
p2<-ggplot(do.call(rbind,points),aes(truth,position,color=stage))+geom_point(alpha=.65,size=.8)+facet_wrap(~id+method,ncol=2)+scale_color_manual(values=c(initial='#999999',final='#0072B2'))+labs(title='Largest improvement and regression, with raw Isomap initialization',x='True latent position',y='Estimated position (orientation aligned, rescaled)',color='Stage',caption='Cases selected after fitting by largest and smallest recovery change; all cases appear in the complete comparison.')+theme(legend.position='bottom')
ggsave(file.path(out,'examples.png'),p2,width=12,height=8,dpi=170,bg='white')
writeLines(capture.output(sessionInfo()),file.path(out,'sessionInfo.txt'))
