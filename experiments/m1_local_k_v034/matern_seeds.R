# Twenty prespecified new seeds; vary GP draws while retaining positions/noise.
source('experiments/m2_local_k_v034/common.R')
library(ggplot2)
out<-'experiments/m1_local_k_v034/matern_seeds';dir.create(out,recursive=TRUE,showWarnings=FALSE)
d<-readRDS(file.path(study,'data','rich_S4_r06.rds'));cols<-which(d$truth==1);truth<-d$latent[,1]
# Match the preceding generator exactly, including its display-grid draws.
tt<-c(truth,seq(0,1,length.out=601));z<-sqrt(5)*abs(outer(tt,tt,'-'))/.2
L<-t(chol((1+z+z^2/3)*exp(-z)+diag(1e-10,length(tt))))
seeds<-20260931L+0:19
writeLines(c('Twenty new GP seeds: 20260931 through 20260950, all retained.',
 'Matern nu=2.5, length-scale=0.2; same latent positions and observation noise as the preceding example.',
 'Each feature centered/scaled to unit sample variance; N300, D12, SNR4.',
 'PCA initialization, K50, fixed uniform prior, quantile discretization, relative tolerance 1e-6.',
 'This measures trajectory-draw sensitivity conditional on a fixed position/noise realization.'),file.path(out,'design.txt'))
rows<-list();fits<-list()
for(seed in seeds){
 set.seed(seed);draw<-L%*%matrix(rnorm(length(tt)*12),length(tt),12)
 signal<-scale(draw[1:300,]);X<-signal+d$unit_noise[,cols]/2
 raw<-as.numeric(MPCurver:::PCA_ordering(X)$t)
 init<-MPCurver:::.cavi_exact_k_init_from_ordering(X=X,ordering_vec=raw,K=50L,discretization='quantile',strict_K=TRUE)
 gamma0<-MPCurver:::.cavi_hard_gamma_from_cluster_rank(init$cluster_rank,300,50L)
 f<-MPCurver:::cavi(X=X,K=50L,responsibilities_init=gamma0,position_prior_init=rep(.02,50),position_prior='fixed',sigma2_init=init$sigma2,lambda_init=rep(1,12),rw_q=2L,ridge=0,discretization='quantile',max_iter=2000L,tol=1e-6,convergence='relative',verbose=FALSE)
 while(!f$converged && f$iter<10000L)f<-MPCurver:::do_cavi(f,iter=2000L,tol=1e-6)
 before<-(raw-min(raw))/diff(range(raw));after<-as.numeric(f$gamma%*%seq(0,1,length.out=50))
 if(cor(before,truth,method='spearman')<0){before<-1-before;after<-1-after}
 row<-data.frame(seed=seed,before=abs(cor(before,truth,method='spearman')),after=abs(cor(after,truth,method='spearman')),iterations=f$iter,converged=f$converged)
 rows[[as.character(seed)]]<-row;fits[[as.character(seed)]]<-list(X=X,signal=signal,truth=truth,before=before,after=after,fit=f,input_hash=digest(X,algo='sha256'))
 print(row);flush.console()
}
tab<-do.call(rbind,rows);stopifnot(all(tab$converged))
write.csv(tab,file.path(out,'summary.csv'),row.names=FALSE)
saveRDS(list(results=fits,settings=list(seeds=seeds,nu=2.5,length_scale=.2,source_hash=digest(file='experiments/m1_local_k_v034/matern_seeds.R',algo='sha256'))),file.path(out,'results.rds'))
theme_set(theme_minimal(base_size=12))
p<-ggplot(tab,aes(before,after))+geom_abline(slope=1,intercept=0,color='gray65')+geom_point(size=2)+geom_hline(yintercept=.95,linetype=2,color='gray65')+geom_vline(xintercept=.95,linetype=2,color='gray65')+coord_equal(xlim=c(0,1),ylim=c(0,1))+labs(title='PCA recovery across 20 new Matern trajectory seeds',subtitle='Matern 5/2 | length-scale = 0.2 | SNR = 4 | fixed uniform prior | K = 50',x='PCA initialization: absolute Spearman correlation',y='Converged MPCurve: absolute Spearman correlation',caption='Figure 1. Each point is one prespecified seed; all 20 are retained. Dashed lines mark recovery 0.95.\nTrue positions and observation noise are held fixed; only the 12 GP trajectories vary.')
ggsave(file.path(out,'recovery.png'),p,width=7.5,height=6.5,dpi=170,bg='white')
# Show all seeds, without selecting favorable examples.
pts<-do.call(rbind,lapply(seq_len(nrow(tab)),function(i){r<-fits[[as.character(tab$seed[i])]];data.frame(truth=rep(truth,2),estimate=c(r$before,r$after),stage=rep(c('PCA','MPCurve'),each=300),seed=sprintf('%d: %.2f -> %.2f',tab$seed[i],tab$before[i],tab$after[i]))}))
p2<-ggplot(pts,aes(truth,estimate,color=stage))+geom_point(size=.35,alpha=.5)+facet_wrap(~seed,ncol=5)+scale_color_manual(values=c(PCA='gray65',MPCurve='#0072B2'))+theme(legend.position='bottom')+labs(title='All 20 seeds: initial and final positions',x='True latent position',y='Estimated latent position',caption='Figure 2. PCA scores are scaled to [0,1]; MPCurve positions are posterior means. Labels show initial -> final recovery.')
ggsave(file.path(out,'all_seeds.png'),p2,width=14,height=10,dpi=170,bg='white')
print(summary(tab[,c('before','after')]));cat('PCA >=.9:',sum(tab$before>=.9),' PCA >=.95:',sum(tab$before>=.95),' final >=.95:',sum(tab$after>=.95),'\n')
