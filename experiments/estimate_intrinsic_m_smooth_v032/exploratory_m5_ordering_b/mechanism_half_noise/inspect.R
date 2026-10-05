source('experiments/estimate_intrinsic_m_smooth_v032/auto_m_v034/common.R')
base <- file.path(study_dir,'exploratory_m5_ordering_b')
out <- file.path(base,'mechanism_half_noise')
X <- as.matrix(read.csv(file.path(base,'external_methods/inputs/B_noise05_X.csv')))
pos <- read.csv(file.path(base,'external_methods/inputs/B_noise05_positions.csv'))
d <- load_fixed_dataset(manifest[manifest$id=='main_M5_S4_r001',])
anchor <- colnames(d$X)[intersect(which(d$true_assign=='B'),d$anchor_indices)]
writeLines(c(deparse(MPCurver:::cavi),deparse(MPCurver:::initialize_ordering_csmooth)),file.path(out,'frozen_cavi_source.R.txt'))
results <- list(); rows <- list()
for(i in c(0,1,2,5,10,20,50,2000)) {
 fit <- MPCurver::fit_mpcurve(X,method='PCA',intrinsic_dim=1,iter=i)$fit
 edf <- vapply(seq_len(ncol(X)),function(j)sum(diag(fit$posterior$cov[[j]])*colSums(fit$gamma)/fit$params$sigma2[j]),numeric(1))
 rows[[length(rows)+1L]] <- data.frame(budget=i,iteration=fit$iter,feature=colnames(X),anchor=colnames(X)==anchor,
  sigma2=fit$params$sigma2,lambda=fit$lambda_vec,conditional_edf=edf,
  rho=abs(cor(pos$truth,as.numeric(fit$gamma %*% seq(0,1,length.out=50)),method='spearman')))
 results[[as.character(i)]] <- fit
}
saveRDS(results,file.path(out,'native_prefix_fits.rds'))
write.csv(do.call(rbind,rows),file.path(out,'native_feature_diagnostics.csv'),row.names=FALSE)
print(do.call(rbind,rows)[do.call(rbind,rows)$budget %in% c(0,2000),],row.names=FALSE)
cat('Anchor:',anchor,'\n')
