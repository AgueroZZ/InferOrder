library(Rtsne)
library(umap)

set.seed(1)
source("./code/plot_ordering.R")
load("./data/loading_order/pancreas_factors.rdata")
load("./data/loading_order/pancreas.rdata")
cells <- subsample_cell_types(sample_info$celltype,n = 500)
Loadings <- fl_snmf_ldf$L[cells,]
celltype <- as.character(sample_info$celltype[cells])
names(celltype) <- rownames(Loadings)
batchtype <- as.character(sample_info$tech[cells])
names(batchtype) <- rownames(Loadings)

# tSNE
tsne <- Rtsne(Loadings, dims = 2, perplexity = 30, verbose = TRUE, check_duplicates = FALSE)
# plot tSNE
tsne_df <- data.frame(
  Dim1 = tsne$Y[,1],
  Dim2 = tsne$Y[,2],
  celltype = celltype,
  batchtype = batchtype
)
ggplot(tsne_df, aes(x = Dim1, y = Dim2, colour = celltype)) +
  geom_point(size = 2, alpha = 0.7) +
  theme_minimal() +
  labs(title = "t-SNE Visualization",
       x = "t-SNE 1", y = "t-SNE 2")
ggplot(tsne_df, aes(x = Dim1, y = Dim2, colour = batchtype)) +
  geom_point(size = 2, alpha = 0.7) +
  theme_minimal() +
  labs(title = "t-SNE Visualization",
       x = "t-SNE 1", y = "t-SNE 2")


# PCA
pca <- prcomp(Loadings, center = TRUE, scale. = FALSE)
# plot PCA
pca_df <- data.frame(
  PC1 = pca$x[,1],
  PC2 = pca$x[,2],
  celltype = celltype,
  batchtype = batchtype
)
ggplot(pca_df, aes(x = PC1, y = PC2, colour = celltype)) +
  geom_point(size = 2, alpha = 0.7) +
  theme_minimal() +
  labs(title = "PCA Visualization",
       x = "PC1", y = "PC2")
ggplot(pca_df, aes(x = PC1, y = PC2, colour = batchtype)) +
  geom_point(size = 2, alpha = 0.7) +
  theme_minimal() +
  labs(title = "PCA Visualization",
       x = "PC1", y = "PC2")


# Additive SmoothEM
source("./code/additive_EM.R")
K <- 20
Q1 <- make_random_walk_precision(K = K, d = ncol(Loadings), q = 2, lambda = 1, ridge = 0)
Q2 <- make_random_walk_precision(K = K, d = ncol(Loadings), q = 2, lambda = 1, ridge = 1e-10)
Q_list <- list(Q1, Q2)
set.seed(123)
init_params <- make_default_init_add(X = Loadings, Q_list = Q_list)
result_RW <- additive_smoothEM(
  data = Loadings,
  Q_list = Q_list,
  init_params = init_params,
  max_iter = 100,
  tol = 1e-4,
  relative_lambda = TRUE,
  verbose = TRUE
)



K1 <- K; K2 <- K
cluster_assign <- max.col(result_RW$gamma)
latent_grid <- expand.grid(i = seq_len(K1), j = seq_len(K2))
latent_grid$cluster <- seq_len(nrow(latent_grid))  # 1:(K1*K2)
obs_map <- data.frame(
  obs     = seq_len(nrow(Loadings)),
  cluster = cluster_assign
)
assign_df <- merge(obs_map, latent_grid, by = "cluster")
raw_cnt <- assign_df %>%
  count(i, j, name = "n_obs")
full_grid <- expand.grid(i = seq_len(K1), j = seq_len(K2))
cnt_full <- full_grid %>%
  left_join(raw_cnt, by = c("i", "j")) %>%
  replace_na(list(n_obs = 0))
ggplot(cnt_full, aes(x = j / K2, y = i / K1, fill = n_obs)) +
  geom_tile() +
  scale_fill_viridis_c(option = "D", name = "# obs") +
  labs(
    x     = "latent j / K2",
    y     = "latent i / K1",
    title = "Observations per Latent Grid Cell"
  ) +
  coord_equal() +
  theme_minimal()



library(ggplot2)
library(dplyr)
library(tidyr)
library(purrr)
library(patchwork)

plot_obs_label_overlay <- function(
    obs_label_vec,
    highlight_labels,
    assign_df,
    K1, K2,
    jitter_width = 0.2,
    jitter_height = 0.2
) {
  if (length(obs_label_vec) != nrow(assign_df)) {
    stop("Length of 'obs_label_vec' must match the number of rows in 'assign_df'")
  }

  assign_df$obs_label <- obs_label_vec

  df_plot <- assign_df %>%
    filter(obs_label %in% highlight_labels)

  n_colors <- length(highlight_labels)
  palette_fun <- scales::hue_pal()(n_colors)

  ggplot(df_plot, aes(x = j / K2, y = i / K1, color = obs_label)) +
    geom_jitter(width = jitter_width / K2, height = jitter_height / K1, size = 1.2, alpha = 0.8) +
    scale_color_manual(values = setNames(palette_fun, highlight_labels)) +
    coord_equal() +
    labs(
      title = "Highlighted Labels in Latent Grid",
      x = "j / K2", y = "i / K1",
      color = "Label"
    ) +
    theme_minimal() +
    theme(legend.position = "right")
}

plot_obs_label_density <- function(obs_label_vec, highlight_labels, assign_df, K1, K2) {
  if (length(obs_label_vec) != nrow(assign_df)) {
    stop("Length of 'obs_label_vec' must match the number of rows in 'assign_df'")
  }

  assign_df$obs_label <- obs_label_vec

  plot_list <- map(highlight_labels, function(lab) {
    # count how many obs of this label fall in each (i,j)
    df_lab <- assign_df %>%
      filter(obs_label == lab) %>%
      count(i, j, name = "n_label")

    # make full grid for missing tiles
    full_grid <- expand.grid(i = seq_len(K1), j = seq_len(K2))
    df_full <- full_grid %>%
      left_join(df_lab, by = c("i", "j")) %>%
      replace_na(list(n_label = 0))

    ggplot(df_full, aes(x = j / K2, y = i / K1, fill = n_label)) +
      geom_tile(color = NA) +
      scale_fill_viridis_c(option = "D", name = "# obs") +
      coord_equal() +
      labs(
        title = paste("Label:", lab),
        x = "j / K2", y = "i / K1"
      ) +
      theme_minimal()
  })

  wrap_plots(plotlist = plot_list, nrow = 2)
}


plot_obs_label_density(
  obs_label_vec    = celltype,
  highlight_labels = c("acinar", "ductal", "alpha", "beta"),
  assign_df        = assign_df,
  K1 = K, K2 = K
)

plot_obs_label_density(
  obs_label_vec    = batchtype,
  highlight_labels = c("inDrop1", "inDrop3", "smartseq2", "celseq2"),
  assign_df        = assign_df,
  K1 = K, K2 = K
)

plot_obs_label_overlay(
  obs_label_vec    = celltype,
  highlight_labels = unique(celltype),
  assign_df        = assign_df,
  K1 = K, K2 = K
)

plot_obs_label_overlay(
  obs_label_vec    = batchtype,
  highlight_labels = unique(batchtype),
  assign_df        = assign_df,
  K1 = K, K2 = K
)


