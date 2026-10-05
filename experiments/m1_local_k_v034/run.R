# Each source ordering is now a standalone M=1 dataset; no clustering or joint fit.
source('experiments/m2_local_k_v034/common.R')
out <- 'experiments/m1_local_k_v034'
indices <- as.integer(commandArgs(trailingOnly=TRUE))
for(i in indices) {
 d <- readRDS(file.path(study,'data',paste0(manifest$id[i],'.rds')))
 for(g in 1:2) {
  id <- paste0(manifest$id[i],'_',LETTERS[g]);path<-file.path(out,'results',paste0(id,'.rds'))
  if(file.exists(path))next
  X<-d$X[,d$truth==g,drop=FALSE];truth<-d$latent[,g]
  rows<-list();positions<-list();traces<-list()
  for(k in design$k_candidates) {
   warnings<-character();start<-proc.time()[['elapsed']]
   iso<-withCallingHandlers(MPCurver:::isomap_ordering(X,k=k),warning=function(w){warnings<<-c(warnings,conditionMessage(w));invokeRestart('muffleWarning')})
   embedding<-proc.time()[['elapsed']]-start
   valid<-all(is.finite(iso$t)) && iso$n_components==1
   if(!valid){rows[[as.character(k)]]<-data.frame(k=k,eligible=FALSE,score=NA,initial_rho=NA,one_rho=NA,final_rho=NA,converged=FALSE,iterations=NA,embedding=embedding,scoring=0,fitting=0);next}
   start<-proc.time()[['elapsed']];one<-local_fit(X,iso$t,1L);scoring<-proc.time()[['elapsed']]-start
   stopifnot(one$iter==1L,length(one$elbo_trace)==2L)
   start<-proc.time()[['elapsed']];fit<-local_fit(X,iso$t,2000L)
   while(!fit$converged && fit$iter<10000L)fit<-MPCurver:::do_cavi(fit,iter=2000L,tol=1e-6)
   fitting<-proc.time()[['elapsed']]-start
   onepos<-as.numeric(one$gamma%*%seq(0,1,length.out=50));final<-as.numeric(fit$gamma%*%seq(0,1,length.out=50))
   rows[[as.character(k)]]<-data.frame(k=k,eligible=TRUE,score=tail(one$elbo_trace,1),initial_rho=abs(cor(truth,iso$t,method='spearman')),one_rho=abs(cor(truth,onepos,method='spearman')),final_rho=abs(cor(truth,final,method='spearman')),converged=fit$converged,iterations=fit$iter,embedding=embedding,scoring=scoring,fitting=fitting)
   positions[[as.character(k)]]<-list(initial=iso$t,one=onepos,final=final,warnings=warnings)
   traces[[as.character(k)]]<-fit$elbo_trace
  }
  tab<-do.call(rbind,rows);chosen<-tab$k[which.max(tab$score)]
  atomic_save(list(id=id,row=manifest[i,],group=g,truth=truth,input_hash=digest(X,algo='sha256'),candidates=tab,selected_k=chosen,positions=positions,traces=traces,source_hash=digest(file=file.path(out,'run.R'),algo='sha256'),package_version=as.character(packageVersion('MPCurver'))),path)
  cat(id,'default',tab$final_rho[tab$k==15],'selected',tab$final_rho[tab$k==chosen],'k',chosen,'\n');flush.console()
 }
}
