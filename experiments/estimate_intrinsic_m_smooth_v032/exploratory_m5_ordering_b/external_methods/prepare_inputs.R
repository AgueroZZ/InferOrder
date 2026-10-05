source('experiments/estimate_intrinsic_m_smooth_v032/auto_m_v034/common.R')
out <- file.path(study_dir,'exploratory_m5_ordering_b','external_methods')
dir.create(file.path(out,'inputs'),showWarnings=FALSE)
d <- load_fixed_dataset(manifest[manifest$id=='main_M5_S4_r001',])
cases <- data.frame(id=c('B_noise0','B_noise05','B_noise1','A_noise1','C_noise1','D_noise1','E_noise1'),
 ordering=c('B','B','B','A','C','D','E'),noise_scale=c(0,.5,1,1,1,1,1))
metadata <- list()
for(i in seq_len(nrow(cases))) {
 row <- cases[i,]; columns <- which(d$true_assign==row$ordering)
 X <- d$signal[,columns]+row$noise_scale*(d$X[,columns]-d$signal[,columns])
 truth <- d$latent_positions[,match(row$ordering,d$ordering_labels)]
 starts <- data.frame(truth=truth,pca=MPCurver:::PCA_ordering(X)$t,
  isomap10=MPCurver:::isomap_ordering(X,k=10L)$t,isomap15=MPCurver:::isomap_ordering(X,k=15L)$t)
 write.csv(X,file.path(out,'inputs',paste0(row$id,'_X.csv')),row.names=FALSE)
 write.csv(starts,file.path(out,'inputs',paste0(row$id,'_positions.csv')),row.names=FALSE)
 metadata[[row$id]] <- list(case=row,columns=columns,features=colnames(X),
  original_input_hash=d$input_hash,subset_hash=digest::digest(X,algo='sha256'))
}
write.csv(cases,file.path(out,'cases.csv'),row.names=FALSE)
saveRDS(metadata,file.path(out,'input_provenance.rds'))
