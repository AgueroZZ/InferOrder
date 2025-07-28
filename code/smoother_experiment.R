smoother_mgcv_generator <- function(K, pi = NULL, bs = "bs", num_basis = 30, locations = NULL, m = 2, sp = NULL, method = "ML") {
  if (is.null(locations)) {
    locations <- seq(0, 1, length.out = K)
  }
  if (is.null(pi)) {
    pi <- rep(1/K, K)
  } else if (length(pi) != K) {
    stop("pi must be a vector of length K")
  }

  smoother <- function(Z, X, params = NULL) {
    z_idx <- apply(Z, 1, which.max)
    t <- locations[z_idx]

    d <- ncol(X)
    U <- matrix(0, nrow = K, ncol = d)
    sigma <- numeric(d)

    for (i in 1:d) {
      df_current <- data.frame(x = X[, i], t = t)
      mod <- mgcv::gam(x ~ s(t, bs = bs, k = num_basis, m = m, sp = sp),
                       method = method,
                       family = gaussian(),
                       data = df_current)
      # mgcv::predict.gam returns a vector
      U[, i] <- predict(mod, newdata = data.frame(t = locations))
      sigma[i] <- var(mod$residuals)
    }

    list(pi = pi, mu = split(U, row(U)), sigma = sigma)
  }

  return(smoother)
}

smoother_iwp_generator <- function(K, pi = NULL, order = 2, locations = NULL, num_basis = 30, sd.prior = NULL){
  if (is.null(locations)) {
    locations <- seq(0, 1, length.out = K)
  }
  if (is.null(pi)) {
    pi <- rep(1/K, K)
  } else if (length(pi) != K) {
    stop("pi must be a vector of length K")
  }

  smoother <- function(Z, X, params = NULL) {
    z_idx <- apply(Z, 1, which.max)
    t <- locations[z_idx]

    d <- ncol(X)
    U <- matrix(0, nrow = K, ncol = d)
    sigma <- numeric(d)

    for (i in 1:d) {
      df_current <- data.frame(x = X[, i], t = t)
      mod <- BayesGP::model_fit(x ~ f(t, model = "IWP",
                                      sd.prior = sd.prior,
                                      order = order, k = num_basis),
                       family = "Gaussian",
                       data = df_current)
      U[, i] <- predict(mod, variable = "t", newdata = data.frame(t = locations))$mean
      sigma[i] <- BayesGP::post_table(mod)$median[3]^2
    }

    list(pi = pi, mu = split(U, row(U)), sigma = sigma)
  }

  return(smoother)
}


source("./code/simulate.R")
library(MASS)
library(mvtnorm)
library(Matrix)

n_mixture <- 20
palette_colors <- rainbow(n_mixture)
alpha_colors <- sapply(palette_colors, function(clr) adjustcolor(clr, alpha.f=0.3))
sim <- simulate_mixture(n=500, K = n_mixture, d=2, seed=123, proj_mat = matrix(c(1,-0.6,-0.6,1), nrow = 2, byrow = T))
plot(sim$X, col = alpha_colors[sim$z],
     pch = 19, cex = 0.5,
     xlab = "X1", ylab = "X2",
     main = "Simulated mixture of Gaussians in 2D")


source("./code/prior_precision.R")
Q_prior_RW1 <- make_random_walk_precision(K=50, d=2, lambda = 10)
set.seed(1234)
init_params <- make_default_init(sim$X, K=50)
result_RW1 <- stochastic_EM(
  X = sim$X,
  Q_prior = Q_prior_RW1,
  init_params = init_params,
  threshold = NULL,
  max_iter = 100,
  tol = 1e-5,
  verbose = TRUE
)
plot(sim$X, col=alpha_colors[sim$z], cex=1,
     xlab="X1", ylab="X2",
     pch=19, main="EM Fitted Means with RW1 Prior")
mu_matrix <- do.call(rbind, result_RW1$params$mu)
for (k in 1:(nrow(mu_matrix)-1)) {
  if (sqrt(sum((mu_matrix[k+1,] - mu_matrix[k,])^2)) > 1e-6) {
    arrows(mu_matrix[k,1], mu_matrix[k,2],
           mu_matrix[k+1,1], mu_matrix[k+1,2],
           col="orange", lwd=2, length=0.1)
  }
}
points(mu_matrix, pch=8, cex=1, lwd=1, col="orange")

plot(mu_matrix[,1] ~ seq(0,1, length.out = nrow(mu_matrix)),
     type = "b", col = "blue", pch = 19,
     xlab = "Index", ylab = "Mean",
     main = "Means from RW1 Prior",
     ylim = range(mu_matrix))
plot(mu_matrix[,2] ~ seq(0,1, length.out = nrow(mu_matrix)),
     type = "b", col = "blue", pch = 19,
     xlab = "Index", ylab = "Mean",
     main = "Means from RW1 Prior",
     ylim = range(mu_matrix))



smoother_mgcv <- smoother_mgcv_generator(K = 50, bs = "bs",
                                         num_basis = 20, m = c(4,3),
                                         sp = 100, method = "REML")
set.seed(123)
init_params <- make_default_init(sim$X, K=50)
result_mgcv <- stochastic_EM(
  X = sim$X,
  init_params = init_params,
  threshold = "hard",
  smoother = smoother_mgcv,
  max_iter = 100,
  tol = 1e-5,
  verbose = TRUE
)
plot(sim$X, col=alpha_colors[sim$z], cex=1,
     xlab="X1", ylab="X2",
     pch=19, main="EM Fitted Means with MGCV Smoothing (Fixed)")
mu_matrix <- do.call(rbind, result_mgcv$params$mu)
for (k in 1:(nrow(mu_matrix)-1)) {
  if (sqrt(sum((mu_matrix[k+1,] - mu_matrix[k,])^2)) > 1e-6) {
    arrows(mu_matrix[k,1], mu_matrix[k,2],
           mu_matrix[k+1,1], mu_matrix[k+1,2],
           col="orange", lwd=2, length=0.1)
  }
}
points(mu_matrix, pch=8, cex=1, lwd=1, col="orange")
plot(mu_matrix[,1] ~ seq(0,1, length.out = nrow(mu_matrix)),
     type = "b", col = "blue", pch = 19,
     xlab = "Index", ylab = "Mean",
     main = "Means from MGCV Smoother",
     ylim = range(mu_matrix))
plot(mu_matrix[,2] ~ seq(0,1, length.out = nrow(mu_matrix)),
     type = "b", col = "blue", pch = 19,
     xlab = "Index", ylab = "Mean",
     main = "Means from MGCV Smoother",
     ylim = range(mu_matrix))



smoother_iwp <- smoother_iwp_generator(K = 50, order = 2, sd.prior = 0.5,
                                         num_basis = 20)
set.seed(1234)
init_params <- make_default_init(sim$X, K=50)
result_iwp <- stochastic_EM(
  X = sim$X,
  init_params = init_params,
  threshold = "hard",
  smoother = smoother_iwp,
  max_iter = 100,
  tol = 1e-2,
  verbose = TRUE
)
plot(sim$X, col=alpha_colors[sim$z], cex=1,
     xlab="X1", ylab="X2",
     pch=19, main="EM Fitted Means with IWP prior (fully Bayes)")
mu_matrix <- do.call(rbind, result_iwp$params$mu)
for (k in 1:(nrow(mu_matrix)-1)) {
  if (sqrt(sum((mu_matrix[k+1,] - mu_matrix[k,])^2)) > 1e-6) {
    arrows(mu_matrix[k,1], mu_matrix[k,2],
           mu_matrix[k+1,1], mu_matrix[k+1,2],
           col="orange", lwd=2, length=0.1)
  }
}
points(mu_matrix, pch=8, cex=1, lwd=1, col="orange")

plot(mu_matrix[,1] ~ seq(0,1, length.out = nrow(mu_matrix)),
     type = "b", col = "blue", pch = 19,
     xlab = "Index", ylab = "Mean",
     main = "Means from IWP Smoother",
     ylim = range(mu_matrix))
plot(mu_matrix[,2] ~ seq(0,1, length.out = nrow(mu_matrix)),
     type = "b", col = "blue", pch = 19,
     xlab = "Index", ylab = "Mean",
     main = "Means from IWP Smoother",
     ylim = range(mu_matrix))



smoother_tf_generator <- function(K, pi = NULL, ord = 1, locations = NULL, lambda = NULL) {
  if (is.null(locations)) {
    locations <- seq(0, 1, length.out = K)
  }
  if (is.null(pi)) {
    pi <- rep(1 / K, K)
  } else if (length(pi) != K) {
    stop("pi must be a vector of length K")
  }

  smoother <- function(Z, X, params = NULL) {
    z_idx <- apply(Z, 1, which.max)
    t <- locations[z_idx]

    d <- ncol(X)
    U <- matrix(NA, nrow = K, ncol = d)
    sigma <- numeric(d)

    for (i in 1:d) {
      # Sort by pseudotime
      o <- order(t)
      t_sorted <- t[o]
      x_sorted <- X[o, i]

      # Genlasso trend filtering only works on integer grid
      # We map t_sorted to equally spaced grid
      n <- length(t_sorted)
      grid_idx <- 1:n

      tf_fit <- genlasso::trendfilter(x_sorted, ord = ord)

      # Choose lambda
      if (is.null(lambda)) {
        # Use cross-validation
        cv <- genlasso::cv.trendfilter(tf_fit)
        lambda_best <- cv$lambda.min
      } else {
        lambda_best <- lambda
      }

      # Fitted values on grid_idx
      beta_hat <- predict(tf_fit, lambda = lambda_best)$fit

      # Now, we want to evaluate at 'locations'
      # We'll linearly interpolate
      approx_fit <- approx(t_sorted, beta_hat, xout = locations, rule = 2)$y

      U[, i] <- approx_fit
      sigma[i] <- var(x_sorted - beta_hat)
    }

    list(pi = pi, mu = split(U, row(U)), sigma = sigma)
  }

  return(smoother)
}
smoother_tf <- smoother_tf_generator(K = 50, ord = 1, lambda = 1)
set.seed(123)
init_params <- make_default_init(sim$X, K=50)
result_tf <- stochastic_EM(
  X = sim$X,
  init_params = init_params,
  threshold = "hard",
  smoother = smoother_tf,
  max_iter = 100,
  tol = 1e-5,
  verbose = TRUE
)

plot(sim$X, col=alpha_colors[sim$z], cex=1,
     xlab="X1", ylab="X2",
     pch=19, main="EM Fitted Means with Trend Filtering")
mu_matrix <- do.call(rbind, result_tf$params$mu)
for (k in 1:(nrow(mu_matrix)-1)) {
  if (sqrt(sum((mu_matrix[k+1,] - mu_matrix[k,])^2)) > 1e-6) {
    arrows(mu_matrix[k,1], mu_matrix[k,2],
           mu_matrix[k+1,1], mu_matrix[k+1,2],
           col="orange", lwd=2, length=0.1)
  }
}
points(mu_matrix, pch=8, cex=1, lwd=1, col="orange")

plot(mu_matrix[,1] ~ seq(0,1, length.out = nrow(mu_matrix)),
     type = "b", col = "blue", pch = 19,
     xlab = "Index", ylab = "Mean",
     main = "Means from Trend Filtering Smoother",
     ylim = range(mu_matrix))
plot(mu_matrix[,2] ~ seq(0,1, length.out = nrow(mu_matrix)),
     type = "b", col = "blue", pch = 19,
     xlab = "Index", ylab = "Mean",
     main = "Means from Trend Filtering Smoother",
     ylim = range(mu_matrix))



library(fda)

smoother_monotone_generator <- function(
    K,
    pi = NULL,
    locations = NULL,
    rangeval = c(0, 1),
    nbasis = 15,
    lambda = 1e-2,
    iterlim = 5,
    dbglev = 0
) {
  # Default grid
  if (is.null(locations)) {
    locations <- seq(rangeval[1], rangeval[2], length.out = K)
  }
  if (is.null(pi)) {
    pi <- rep(1/K, K)
  } else if (length(pi) != K) {
    stop("pi must be a vector of length K")
  }

  # Return the smoothing function
  smoother <- function(Z, X, params = NULL) {
    z_idx <- apply(Z, 1, which.max)
    t <- locations[z_idx]

    n <- nrow(X)
    d <- ncol(X)

    U <- matrix(NA, nrow = K, ncol = d)
    sigma <- numeric(d)

    # Fit each dimension separately
    for (j in 1:d) {
      y <- X[, j]

      # Create basis
      wbasis <- create.bspline.basis(rangeval, nbasis)
      WfdParobj <- fdPar(wbasis, Lfdobj = 2, lambda = lambda)

      # Fit monotone smoother
      fit <- smooth.monotone(argvals = t, y = y, WfdParobj = WfdParobj, iterlim = iterlim, dbglev = dbglev)

      # Evaluate at grid locations
      monotone_part <- eval.monfd(locations, fit$Wfdobj)
      U[, j] <- fit$beta[1,] + fit$beta[2,] * monotone_part

      # Estimate residual variance
      fitted_t <- eval.monfd(t, fit$Wfdobj)
      fitted_y <- fit$beta[1,] + fit$beta[2,] * fitted_t
      sigma[j] <- var(y - fitted_y)
    }

    list(pi = pi, mu = split(U, row(U)), sigma = sigma)
  }

  return(smoother)
}

smoother_monotone <- smoother_monotone_generator(
  K = 50, rangeval = c(0, 1), nbasis = 20, lambda = 1e-1, iterlim = 5
)
set.seed(123)
init_params <- make_default_init(sim$X, K=50)
result_monotone <- stochastic_EM(
  X = sim$X,
  init_params = init_params,
  threshold = "hard",
  smoother = smoother_monotone,
  max_iter = 100,
  tol = 1e-5,
  verbose = TRUE
)

plot(sim$X, col=alpha_colors[sim$z], cex=1,
     xlab="X1", ylab="X2",
     pch=19, main="EM Fitted Means with Monotone Smoothing")
mu_matrix <- do.call(rbind, result_monotone$params$mu)
for (k in 1:(nrow(mu_matrix)-1)) {
  if (sqrt(sum((mu_matrix[k+1,] - mu_matrix[k,])^2)) > 1e-6) {
    arrows(mu_matrix[k,1], mu_matrix[k,2],
           mu_matrix[k+1,1], mu_matrix[k+1,2],
           col="orange", lwd=2, length=0.1)
  }
}
points(mu_matrix, pch=8, cex=1, lwd=1, col="orange")
plot(mu_matrix[,1] ~ seq(0,1, length.out = nrow(mu_matrix)),
     type = "b", col = "blue", pch = 19,
     xlab = "Index", ylab = "Mean",
     main = "Means from Monotone Smoother",
     ylim = range(mu_matrix))
plot(mu_matrix[,2] ~ seq(0,1, length.out = nrow(mu_matrix)),
     type = "b", col = "blue", pch = 19,
     xlab = "Index", ylab = "Mean",
     main = "Means from Monotone Smoother",
     ylim = range(mu_matrix))






