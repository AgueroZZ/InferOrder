# simulate mixture data with a given projection matrix
simulate_mixture <- function(n, K, d, seed=NULL, proj_mat = NULL) {
  if (!is.null(seed)) set.seed(seed)

  # True weights
  pi_true <- rep(1/K, K)

  # True means
  mu_list <- lapply(1:K, function(k) {
    # if proj_mat is provided, use it to generate means
    if (!is.null(proj_mat)) {
      return(as.numeric(proj_mat %*% runif(d, -5, 5)))
    }
    else {
      return(runif(d, -5, 5))
    }
  })

  # True covariances (positive definite)
  sigma_list <- lapply(1:K, function(k) {
    Sigma <- diag(0.1, d)
    return(Sigma)
  })

  # Sample assignments
  z <- sample(1:K, n, replace=TRUE, prob=pi_true)

  # Simulate data
  X <- matrix(NA, n, d)
  for (i in 1:n) {
    X[i, ] <- mvrnorm(1, mu_list[[z[i]]], sigma_list[[z[i]]])
  }

  list(
    X = X,
    z = z,
    true_params = list(pi=pi_true, mu=mu_list, sigma=sigma_list)
  )
}
