source('experiments/m2_local_k_v034/common.R')
indices<-as.integer(commandArgs(trailingOnly=TRUE))
stopifnot(length(indices)>0,all(indices %in% 1:nrow(manifest)))
for(index in indices) {
 row<-manifest[index,];output<-file.path(study,'results',paste0(row$id,'.rds'))
 if(file.exists(output)) next
 d<-generate_data(row);atomic_save(d,file.path(study,'data',paste0(row$id,'.rds')))
 started<-proc.time()[['elapsed']]
 similarity<-MPCurver:::.compute_same_ordering_similarity(d$X,metric='spline_r2',spline_r2_df=5L)
 hc<-hclust(as.dist(similarity$distance),method='single')
 inferred<-MPCurver:::.cavi_canonicalize_feature_clusters(cutree(hc,k=2))$feature_cluster
 cluster_seconds<-proc.time()[['elapsed']]-started
 groups<-list(inferred=inferred,truth=d$truth)
 records<-list();candidate_rows<-list();cache<-new.env(parent=emptyenv())
 for(regime in names(groups)) {
  assignment<-groups[[regime]];members<-split(1:24,assignment)
  stopifnot(length(members)==2)
  local<-lapply(seq_along(members),function(g) {
   cols<-members[[g]];key<-paste(cols,collapse=',')
   if(exists(key,envir=cache,inherits=FALSE)) return(get(key,envir=cache))
   X<-d$X[,cols,drop=FALSE];candidates<-list();times<-numeric();diagnostics<-list()
   for(k in design$k_candidates) {
    warnings<-character();error<-'';start<-proc.time()[['elapsed']]
    iso<-tryCatch(withCallingHandlers(MPCurver:::isomap_ordering(X,k=k),
     warning=function(w){warnings<<-c(warnings,conditionMessage(w));invokeRestart('muffleWarning')}),
     error=function(e){error<<-conditionMessage(e);NULL})
    embed_seconds<-proc.time()[['elapsed']]-start;score_seconds<-0
    valid<-!is.null(iso) && all(is.finite(iso$t)) && length(cols)>1
    one<-NULL
    if(valid) {
     start<-proc.time()[['elapsed']]
     one<-tryCatch(local_fit(X,iso$t,1L),error=function(e){error<<-conditionMessage(e);NULL})
     score_seconds<-proc.time()[['elapsed']]-start;valid<-!is.null(one)
    }
    if(valid) stopifnot(one$iter==1L,length(one$elbo_trace)==2)
    candidates[[as.character(k)]]<-list(raw=if(!is.null(iso))iso$t else NULL,fit=one)
    diagnostics[[as.character(k)]]<-data.frame(k=k,eligible=valid,
     score=if(valid)tail(one$elbo_trace,1) else NA_real_,embed_seconds=embed_seconds,
     score_seconds=score_seconds,components=if(!is.null(iso))iso$n_components else NA_integer_,
     warnings=paste(warnings,collapse=' | '),error=error)
   }
   tab<-do.call(rbind,diagnostics)
   selected<-if(any(tab$eligible))tab$k[which.max(replace(tab$score,!tab$eligible,-Inf))] else NA_integer_
   value<-list(candidates=candidates,table=tab,selected=selected,cols=cols)
   assign(key,value,envir=cache);value
  })
  for(g in 1:2) {
   tab<-local[[g]]$table
   for(j in 1:nrow(tab)) {
    one<-local[[g]]$candidates[[as.character(tab$k[j])]]$fit
    tab$rho_A[j]<-if(!is.null(one))abs(cor(d$latent[,1],as.numeric(one$gamma%*%seq(0,1,length.out=50)),method='spearman')) else NA_real_
    tab$rho_B[j]<-if(!is.null(one))abs(cor(d$latent[,2],as.numeric(one$gamma%*%seq(0,1,length.out=50)),method='spearman')) else NA_real_
   }
   candidate_rows[[length(candidate_rows)+1]]<-cbind(id=row$id,regime=regime,group=g,
    true_A_features=sum(d$truth[members[[g]]]==1),true_B_features=sum(d$truth[members[[g]]]==2),tab)
  }
  # Alternate method order across datasets to reduce timing-order artifacts.
  methods<-if(index%%2) c('default15','one_step_k') else c('one_step_k','default15')
  for(method in methods) {
   k_used<-integer(2);fallback<-logical(2);warm_seconds<-0;fits<-list();initial_positions<-matrix(NA_real_,300,2)
   for(g in 1:2) {
    info<-local[[g]];k<-if(method=='default15')15L else info$selected
    valid<-!is.na(k) && info$table$eligible[match(k,info$table$k)]
    start<-proc.time()[['elapsed']]
    if(valid) sub<-local_fit(d$X[,info$cols,drop=FALSE],info$candidates[[as.character(k)]]$raw,2L) else {
     sub<-MPCurver:::.cavi_similarity_subset_fit(X_sub=d$X[,info$cols,drop=FALSE],S=NULL,
      K=50L,method='isomap',pca_component=NA_integer_,rw_q=2L,ridge=0,lambda_sd_prior_rate=NULL,
      lambda_min=1e-10,lambda_max=1e10,sigma_min=1e-10,sigma_max=1e10,max_iter=2L,
      tol=1e-6,discretization='quantile',verbose=FALSE)$fit
     fallback[g]<-TRUE
    }
    fits[[g]]<-expand_fit(d$X,sub);k_used[g]<-if(valid)k else NA_integer_
    initial_positions[,g]<-as.numeric(sub$gamma%*%seq(0,1,length.out=50))
    warm_seconds<-warm_seconds+proc.time()[['elapsed']]-start
   }
   embedding_seconds<-sum(vapply(local,function(x)if(method=='default15')x$table$embed_seconds[x$table$k==15] else sum(x$table$embed_seconds),numeric(1)))
   scoring_seconds<-if(method=='default15')0 else sum(vapply(local,function(x)sum(x$table$score_seconds),numeric(1)))
   warnings<-character();start<-proc.time()[['elapsed']]
   f<-withCallingHandlers(fit_joint(d$X,fits,row$fit_seed),warning=function(w){warnings<<-c(warnings,conditionMessage(w));invokeRestart('muffleWarning')})
   joint_seconds<-proc.time()[['elapsed']]-start
   position<-vapply(f$fit$gamma,function(x)as.numeric(x%*%seq(0,1,length.out=50)),numeric(300))
   correlations<-abs(cor(d$latent,position,method='spearman'))
   matching<-as.integer(clue::solve_LSAP(correlations,maximum=TRUE))
   recovery<-correlations[cbind(1:2,matching)]
   cluster_cost<-if(regime=='inferred')cluster_seconds else 0
   item<-list(id=row$id,regime=regime,method=method,selected_k=k_used,fallback=fallback,
    initial_ARI=adjustedRandIndex(d$truth,assignment),ARI=adjustedRandIndex(d$truth,f$fit$assign),
    recovery=recovery,matching=matching,initial_positions=initial_positions,positions=position,
    weights=f$fit$pi_weights,sigma2=f$params$sigma2,objective=f$fit$objective_history,
    iterations=f$fit$iter,converged=f$fit$converged,warnings=warnings,
    timing=c(clustering=cluster_cost,embedding=embedding_seconds,scoring=scoring_seconds,
     warmup=warm_seconds,joint=joint_seconds,total=cluster_cost+embedding_seconds+scoring_seconds+warm_seconds+joint_seconds))
   records[[paste(regime,method,sep='_')]]<-item
   cat(sprintf('%s %s %s: rho=(%.3f,%.3f) initialARI=%.3f finalARI=%.3f k=%s converged=%s\n',
    row$id,regime,method,recovery[1],recovery[2],item$initial_ARI,item$ARI,paste(k_used,collapse=','),item$converged));flush.console()
  }
 }
 atomic_save(list(row=row,input_hash=d$input_hash,records=records,candidates=do.call(rbind,candidate_rows),
  inferred_groups=inferred,design=design,design_hash=digest(design),
  common_hash=digest(file=file.path(study,'common.R'),algo='sha256'),
  driver_hash=digest(file=file.path(study,'run_dataset.R'),algo='sha256')),output)
}
