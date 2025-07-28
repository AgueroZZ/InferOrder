####################################################################
####################################################################
################### 1. Define functions ############################
####################################################################
####################################################################
####################################################################
#' Fit SmoothEM along a sequence of lambda values using warm‐start initialization
#'
#' This function runs your SmoothEM (e.g. RW2 prior) model repeatedly over a
#' user‐supplied vector of penalty parameters `lambdas`, using the previous
#' fit’s parameters to initialize the next.  Results for each λ are stored
#' in a named list for downstream inspection or plotting.
#'
#' @param data Numeric matrix of loadings (rows = observations, cols = features).
#' @param ordering_vec Numeric vector for initial ordering (e.g. PC1 scores).
#' @param K Integer; number of mixture components.
#' @param lambdas Numeric vector of λ values (should be sorted from largest to smallest).
#' @param q Integer; order of the random‐walk prior (default: 2).
#' @param ridge Numeric; small ridge term for precision‐matrix stability (default: 1e-10).
#' @param modelName Character; passed to `EM_algorithm` (default: "EEI").
#' @param max_iter Integer; maximum EM iterations (default: 500).
#' @param tol Numeric; convergence tolerance for EM (default: 1e-3).
#' @param relative_lambda Logical; treat λ as relative (default: TRUE).
#' @param verbose Logical; print EM progress (default: FALSE).
#'
#' @return Named list of length `length(lambdas)`.  Each element is the full
#'   output of `EM_algorithm`, with additional components:
#'   \describe{
#'     \item{clustering}{Hard cluster labels (max‐responsibility).}
#'     \item{order_idx}{Row indices ordered by `clustering`.}
#'   }
#'
#' @examples
#' lambdas <- c(5e12, 1e12, 5e11, 1e11, 5e10)
#' res_list <- fit_em_path(
#'   data         = Loadings,
#'   ordering_vec = PC1,
#'   K            = 50,
#'   lambdas      = lambdas,
#'   verbose      = TRUE
#' )
#'
fit_em_path <- function(data,
                        ordering_vec = NULL,
                        K,
                        lambdas,
                        q = 2,
                        ridge = 0,
                        modelName = "EEI",
                        max_iter = 200,
                        tol = 1e-3,
                        relative_lambda = TRUE,
                        warmstart = TRUE,
                        verbose = FALSE) {
  # Prepare an empty list to hold results, named by each lambda
  results <- vector("list", length(lambdas))
  names(results) <- paste0("lambda_", format(lambdas, scientific = TRUE))

  # If ordering_vec is not NULL
  if(!is.null(ordering_vec)){
    # Initialize EM parameters once using the supplied ordering vector
    init_params <- make_init(X = data, ordering_vec = ordering_vec, K = K)
  }
  else{
    # If no ordering_vec, use default initialization
    init_params <- make_default_init(X = data, K = K)
  }

  # Loop over each lambda value
  for (i in seq_along(lambdas)) {
    lam <- lambdas[i]

    # Build the random‐walk precision matrix (RW-q)
    Q_prior <- make_random_walk_precision(
      K      = K,
      d      = ncol(data),
      lambda = lam,
      q      = q,
      ridge  = ridge
    )

    # Run the EM algorithm with warm‐start
    res <- EM_algorithm(
      data            = data,
      Q_prior         = Q_prior,
      init_params     = init_params,
      max_iter        = max_iter,
      modelName       = modelName,
      tol             = tol,
      relative_lambda = relative_lambda,
      verbose         = verbose
    )

    # Compute hard clustering from the posterior responsibilities
    res$clustering <- apply(res$gamma, 1, which.max)
    # Store the ordering of rows by cluster label
    res$order_idx  <- order(res$clustering)

    # Save this result
    results[[i]] <- res

    if(warmstart) {
      # Update init_params for the next iteration (warm‐start)
      init_params <- res$params
    }
  }

  return(results)
}


#' Plot selected EM fits using ggplot2 in a grid layout
#'
#' This version uses gridExtra::grid.arrange() to lay out multiple ggplot
#' objects in a true grid.  par(mfrow) has no effect on grid graphics.
#'
#' @param result_lists Named list of EM results (from fit_em_path).
#' @param loadings     Data matrix of loadings (rows = observations).
#' @param indices      Integer vector or names selecting which fits to plot.
#'                     Default = 1 (the first fit).
#' @param ncol         Number of columns in the grid; if NULL, defaults to
#'                     ceiling(sqrt(number_of_plots)).
#' @param palette      Optional vector of fill colors, passed to plot_structure.
#' @param ...          Further args passed on to plot_structure().
#'
#' @return Invisibly NULL; draws a grid of ggplots.
#' @examples
#' # single plot
#' plot_em_path_results(result_lists, Loadings, indices = 1)
#'
#' # 3×3 grid for fits 1–9
#' plot_em_path_results(result_lists, Loadings, indices = 1:9, ncol = 3)
plot_em_path_results <- function(result_lists,
                                 loadings,
                                 indices = 1,
                                 ncol = NULL,
                                 palette = NULL,
                                 ...) {
  # ensure gridExtra is available
  if (!requireNamespace("gridExtra", quietly = TRUE)) {
    stop("Please install the 'gridExtra' package (install.packages('gridExtra'))")
  }

  # --- Determine which fits to plot ---
  sel_idx <- if (is.numeric(indices)) {
    indices
  } else {
    match(indices, names(result_lists))
  }
  sel_idx <- sel_idx[!is.na(sel_idx)]
  if (length(sel_idx) < 1) {
    stop("No valid indices selected in result_lists")
  }

  n_plots <- length(sel_idx)
  if (is.null(ncol)) {
    ncol <- ceiling(sqrt(n_plots))
  }

  # --- Generate ggplot objects ---
  plot_list <- vector("list", n_plots)
  for (j in seq_along(sel_idx)) {
    i   <- sel_idx[j]
    nm  <- names(result_lists)[i]
    res <- result_lists[[i]]

    # compute hard clusters and ordering
    res$clustering <- apply(res$gamma, 1, which.max)
    order_idx      <- order(res$clustering)
    sample_order   <- rownames(loadings)[order_idx]

    # build the ggplot object
    p <- plot_structure(
      loadings,
      plot_name = paste0(nm, "; Obj:", round(max(res$loglik_trace), 3)),
      order     = sample_order,
      palette   = palette,
      ...
    )
    plot_list[[j]] <- p
  }

  # --- Arrange & draw ---
  if (n_plots == 1) {
    print(plot_list[[1]])
  } else {
    gridExtra::grid.arrange(grobs = plot_list, ncol = ncol)
  }

  invisible(NULL)
}




####################################################################
####################################################################
################### 2. Load Data ###################################
####################################################################
####################################################################
####################################################################
library(tibble)
library(tidyr)
library(ggplot2)
library(Rtsne)
library(umap)
set.seed(1)
source("./code/plot_ordering.R")
load("./data/loading_order/pancreas_factors.rdata")
load("./data/loading_order/pancreas.rdata")
cells <- subsample_cell_types(sample_info$celltype,n = 500)
Loadings <- fl_snmf_ldf$L[cells,c(3,8,9,12,17,18,20,21)]
celltype <- as.character(sample_info$celltype[cells])
names(celltype) <- rownames(Loadings)
batchtype <- as.character(sample_info$tech[cells])
names(batchtype) <- rownames(Loadings)
Loadings <- Loadings[batchtype == "inDrop3",]
celltype <- celltype[batchtype == "inDrop3"]
batchtype <- batchtype[batchtype == "inDrop3"]




####################################################################
####################################################################
################### 3. Analysis ####################################
####################################################################
####################################################################
####################################################################
source("./code/prior_precision.R")
source("./code/general_EM.R")
PC1 <- prcomp(Loadings,center = TRUE, scale. = FALSE)$x[,1]
PC1_order <- order(PC1)
plot_structure(Loadings, order = rownames(Loadings)[PC1_order], plot_name = "PC")
lambdas_vec <- c(1e15, 1e10, 1e6, 100, 30, 5, 0.1, 0.001, 0.00001)
set.seed(1)
result_lists_warmstart <- fit_em_path(
  data         = Loadings,
  # ordering_vec = -PC1,
  K            = 100,
  max_iter = 500,
  lambdas      = lambdas_vec,
  ridge = 0,
  modelName = "EII",
  relative_lambda = FALSE,
  verbose      = TRUE
)
set.seed(1)
result_lists_no_warmstart <- fit_em_path(
  data         = Loadings,
  # ordering_vec = -PC1,
  K            = 100,
  max_iter = 500,
  lambdas      = lambdas_vec,
  warmstart     = FALSE,
  relative_lambda = FALSE,
  ridge = 0,
  modelName = "EII",
  verbose      = TRUE
)

lambdas <- as.numeric(sub("^lambda_", "", names(result_lists_warmstart)))
names(result_lists_warmstart) <- paste0("lambda_", format(signif(lambdas, 3), scientific = TRUE))
names(result_lists_no_warmstart) <- paste0("lambda_", format(signif(lambdas, 3), scientific = TRUE))

# 9‐panel (3×3) plot
plot_em_path_results(result_lists_warmstart, Loadings, indices = 1:9)
plot_em_path_results(result_lists_no_warmstart, Loadings, indices = 1:9)


# Check the length for all result
unlist(lapply(result_lists_warmstart, function(x) length(x$loglik_trace)))
unlist(lapply(result_lists_no_warmstart, function(x) length(x$loglik_trace)))



## Conjecture, the smoothed mean are comparable across different local optima
## However, the inferred locations tend to be different.
## Let's confirm if this is true

result1 <- result_lists_warmstart[[1]]
result2 <- result_lists_no_warmstart[[1]]

# obtain the mean
mu1 <- do.call(rbind, result1$params$mu)
mu2 <- do.call(rbind, result2$params$mu)

# look at the first coordinate
plot(mu1[,1], type = "o", col = "blue", ylim = range(c(mu1[,1], mu2[,1])),
     xlab = "Component", ylab = "Mean (1st coordinate)", main = "Comparison of Means")
points(mu2[,1], type = "o", col = "red")

# What if we sort the means?
plot(sort(mu1[,1]), type = "o", col = "blue", ylim = range(c(sort(mu1[,1]), sort(mu2[,1]))),
     xlab = "Component (sorted)", ylab = "Mean (1st coordinate)", main = "Comparison of Sorted Means")
points(sort(mu2[,1]), type = "o", col = "red")

# What if we sort the matrices by the second coordinate?
mu1_ordered <- mu1[order(mu1[,2]),]
mu2_ordered <- mu2[order(mu2[,2]),]

# plot again
plot(mu1_ordered[, 1], type = "o", col = "blue", ylim = range(c(mu1_ordered[,1], mu2_ordered[,1])),
     xlab = "Component (sorted by 2nd coord)", ylab = "Mean (1st coordinate)", main = "Comparison of Ordered Means")
points(mu2_ordered[,1], type = "o", col = "red")

