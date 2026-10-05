source('experiments/m2_local_k_v034/common.R')
suppressPackageStartupMessages({library(ggplot2);library(patchwork)})
paths<-file.path(study,'results',paste0(manifest$id,'.rds'))
stopifnot(all(file.exists(paths)))
results<-lapply(paths,readRDS)
runs<-list();orders<-list();candidates<-list()
for(x in results) {
 stopifnot(identical(x$design_hash,digest(design)))
 for(z in x$records) {
  runs[[length(runs)+1]]<-data.frame(id=z$id,family=x$row$family,snr=x$row$snr,replicate=x$row$replicate,
   regime=z$regime,method=z$method,initial_ARI=z$initial_ARI,ARI=z$ARI,
   mean_rho=mean(z$recovery),min_rho=min(z$recovery),failures=sum(z$recovery<.9),
   strong_recovery=sum(z$recovery>=.95),effective_M=sum(colMeans(z$weights)>1e-12),
   objective=tail(z$objective,1),converged=z$converged,iterations=z$iterations,
   joint_warnings=length(z$warnings),fallbacks=sum(z$fallback),
   clustering_seconds=z$timing['clustering'],embedding_seconds=z$timing['embedding'],
   scoring_seconds=z$timing['scoring'],warmup_seconds=z$timing['warmup'],
   joint_seconds=z$timing['joint'],total_seconds=z$timing['total'])
  for(m in 1:2) orders[[length(orders)+1]]<-data.frame(id=z$id,family=x$row$family,snr=x$row$snr,
   replicate=x$row$replicate,regime=z$regime,method=z$method,ordering=LETTERS[m],rho=z$recovery[m],
   selected_k=z$selected_k[z$matching[m]],initial_ARI=z$initial_ARI,ARI=z$ARI)
 }
 candidates[[length(candidates)+1]]<-x$candidates
}
runs<-do.call(rbind,runs);orders<-do.call(rbind,orders);candidates<-do.call(rbind,candidates)
base<-orders[orders$method=='default15',];selected<-orders[orders$method=='one_step_k',]
pairs<-merge(base,selected,by=c('id','family','snr','replicate','regime','ordering'),suffixes=c('_default','_selected'))
pairs$delta<-pairs$rho_selected-pairs$rho_default
pairs$default_failure<-pairs$rho_default<.9
pairs$rescued<-pairs$default_failure & pairs$rho_selected>=.95
pairs$gain_005<-pairs$delta>.05;pairs$drop_005<-pairs$delta< -.05
pairs$good_to_bad<-pairs$rho_default>=.95 & pairs$rho_selected<.9
base_runs<-runs[runs$method=='default15',];selected_runs<-runs[runs$method=='one_step_k',]
timing<-merge(base_runs,selected_runs,by=c('id','family','snr','replicate','regime'),suffixes=c('_default','_selected'))
timing$total_ratio<-timing$total_seconds_selected/timing$total_seconds_default
summarize_pairs<-function(x)data.frame(orderings=nrow(x),mean_default=mean(x$rho_default),mean_selected=mean(x$rho_selected),
 default_failures=sum(x$default_failure),rescued=sum(x$rescued),selected_failures=sum(x$rho_selected<.9),
 gains_over_005=sum(x$gain_005),drops_over_005=sum(x$drop_005),good_to_bad=sum(x$good_to_bad),
 median_delta=median(x$delta),mean_delta=mean(x$delta))
overall<-do.call(rbind,lapply(split(pairs,pairs$regime),function(x)cbind(regime=x$regime[1],summarize_pairs(x))))
condition<-do.call(rbind,lapply(split(pairs,interaction(pairs$regime,pairs$family,pairs$snr)),
 function(x)cbind(regime=x$regime[1],family=x$family[1],snr=x$snr[1],summarize_pairs(x))))
timing_summary<-do.call(rbind,lapply(split(timing,timing$regime),function(x)data.frame(regime=x$regime[1],
 median_default_total=median(x$total_seconds_default),median_selected_total=median(x$total_seconds_selected),
 median_paired_ratio=median(x$total_ratio),median_extra_seconds=median(x$total_seconds_selected-x$total_seconds_default),
 median_screening_seconds=median(x$embedding_seconds_selected+x$scoring_seconds_selected),
 median_scoring_seconds=median(x$scoring_seconds_selected))))
# Oracle candidate coverage is an evaluation only; it never enters selection.
coverage<-do.call(rbind,lapply(results,function(x)do.call(rbind,lapply(1:2,function(g) {
 tab<-x$candidates[x$candidates$regime=='truth' & x$candidates$group==g & x$candidates$eligible,]
 scores<-tab[[paste0('rho_',LETTERS[g])]]
 best<-if(length(scores))max(scores) else NA_real_
 chosen<-if(nrow(tab))which.max(tab$score) else NA_integer_
 data.frame(id=x$row$id,family=x$row$family,snr=x$row$snr,ordering=LETTERS[g],
  best_candidate_rho=best,selected_candidate_rho=if(nrow(tab))scores[chosen] else NA_real_,
  selected_k=if(nrow(tab))tab$k[chosen] else NA_integer_,
  default_candidate_rho=if(any(tab$k==15))scores[tab$k==15] else NA_real_)
}))))
for(name in c('runs','orders','candidates','pairs','overall','condition','timing','timing_summary','coverage'))
 write.csv(get(name),file.path(study,paste0(name,'.csv')),row.names=FALSE)
theme_set(theme_minimal(base_size=12))
colors<-c('broad'='#0072B2','rich'='#D55E00')
p<-ggplot(pairs,aes(rho_default,rho_selected,color=family))+geom_abline(slope=1,intercept=0,color='gray60')+
 geom_vline(xintercept=.9,linetype=2,color='gray70')+geom_hline(yintercept=.95,linetype=2,color='gray70')+
 geom_point(size=2,alpha=.8)+facet_grid(regime~snr)+scale_color_manual(values=colors)+
 coord_equal(xlim=c(0,1),ylim=c(0,1))+theme(legend.position='bottom')+
 labs(title='Can one local MPCurve update select a better Isomap neighborhood?',
 subtitle='Known M=2; 36 datasets, two orderings each. Feature groups are inferred (top) or supplied only for initialization (bottom).',
 x='Final ordering recovery: default k=15',y='Final ordering recovery: one-step k selection',color='Trajectory family',
 caption='Points above the diagonal improve. Dashed lines mark default failure (<0.90) and rescue (>=0.95).\nOrderings within a dataset and the three SNR versions of each base realization are dependent; these are descriptive counts.')
ggsave(file.path(study,'recovery_comparison.png'),p,width=13,height=9,dpi=170,bg='white')
p2<-ggplot(timing,aes(total_seconds_default,total_seconds_selected,color=family))+geom_abline(slope=1,intercept=0,color='gray50')+
 geom_point(size=2,alpha=.8)+facet_wrap(~regime)+scale_color_manual(values=colors)+theme(legend.position='bottom')+
 labs(title='Required pipeline time, including candidate search',x='Default pipeline time (seconds)',y='Selected-k pipeline time (seconds)',color='Family',
 caption='Times sum separately measured required components, including embeddings, local scores, warmup and joint fitting.\nShared computations are charged to each pipeline as needed; concurrent single-thread timings are descriptive.')
ggsave(file.path(study,'runtime_comparison.png'),p2,width=11,height=5.5,dpi=170,bg='white')
# Display largest rescue and largest regression, selected by the stated rule.
inferred_pairs<-pairs[pairs$regime=='inferred',]
examples<-rbind(inferred_pairs[which.max(inferred_pairs$delta),],inferred_pairs[which.min(inferred_pairs$delta),])
example_points<-list();example_signals<-list()
for(i in 1:nrow(examples)) {
 row<-examples[i,];x<-results[[match(row$id,manifest$id)]];d<-readRDS(file.path(study,'data',paste0(row$id,'.rds')))
 m<-match(row$ordering,LETTERS);label<-sprintf('%s, ordering %s: change %+.3f',row$id,row$ordering,row$delta)
 for(method in c('default15','one_step_k')) {
  z<-x$records[[paste('inferred',method,sep='_')]];pos<-z$positions[,z$matching[m]]
  if(cor(pos,d$latent[,m],method='spearman')<0)pos<-1-pos
  example_points[[length(example_points)+1]]<-data.frame(truth=d$latent[,m],position=pos,example=label,method=method)
 }
 for(j in which(d$truth==m))example_signals[[length(example_signals)+1]]<-data.frame(truth=d$latent[,m],signal=d$signal[,j],example=label,feature=colnames(d$X)[j])
}
p3<-ggplot(do.call(rbind,example_points),aes(truth,position,color=method))+geom_point(size=.85,alpha=.65)+
 facet_wrap(~example,ncol=2)+scale_color_manual(values=c(default15='#0072B2',one_step_k='#D55E00'))+
 theme(legend.position='bottom')+labs(title='Largest observed improvement and regression in the inferred-group regime',
 x='True latent position',y='Inferred position (orientation aligned)',color='Method',
 caption='Examples are selected after fitting by maximum and minimum change in recovery; the complete results include all 36 datasets.')
ggsave(file.path(study,'selected_examples.png'),p3,width=13,height=5.5,dpi=170,bg='white')
p4<-ggplot(do.call(rbind,example_signals),aes(truth,signal,group=feature,color=feature))+geom_line(alpha=.65)+facet_wrap(~example,ncol=2)+
 theme(legend.position='none')+labs(title='All 12 generating trajectories for the displayed orderings',
 x='True latent position',y='Noiseless feature signal (unit sample variance)',
 caption='Each line is an independently drawn Fourier mixture. There are no imposed monotone anchors.')
ggsave(file.path(study,'example_trajectories.png'),p4,width=13,height=5,dpi=170,bg='white')
print(overall,row.names=FALSE);print(condition,row.names=FALSE);print(timing_summary,row.names=FALSE)
print(examples[,c('id','ordering','rho_default','rho_selected','selected_k_selected','delta')],row.names=FALSE)
cat('converged=',sum(runs$converged),'/',nrow(runs),' joint warnings=',sum(runs$joint_warnings),'\n')
