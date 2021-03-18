

#' Obtain the trajectory of LSA-LASSO via path algorithm
#'
#' @import genlasso
#' @export
lasso_path <- function(X, Y, lambda = NULL) {
  # dimension of the parameter
  dd <- ncol(X)
  nn <- nrow(X)
  # OLS estimator
  beta_tilde <- as.numeric(solve(t(X) %*% X) %*% t(X) %*% Y)

  # residuals
  ei <- as.numeric(Y - X %*% beta_tilde)

  # sandwich covariance matrix for OLS estimator
  Sigma_hat <- solve(t(X) %*% X) * mean(ei ^ 2)


  Omega <- solve(Sigma_hat)
  LL <- chol(Omega) # cholesky decomposition of the weight matrix

  gamma_tilde <- LL %*% beta_tilde

  # defining the objective function
  sol_results <- genlasso::genlasso(y = gamma_tilde, D = solve(LL), approx = TRUE, svd = TRUE, rtol = 1e-15, eps = 1e-15)

  ## summarize the solution path at the knots
  sol_path <- vector("list", 3)
  sol_path[[1]] <- c(2 * sol_results$lambda, 0)  ## multiply by 2 due to different parametrization

  sol_path[[2]] <- cbind(sol_results$beta, sol_results$bls)
  colnames(sol_path[[2]])[ncol(sol_path[[2]])] <- "0"
  colnames(sol_path[[2]]) <- round(2 * as.numeric(colnames(sol_path[[2]])), 2)
  sol_path[[2]] <- apply(sol_path[[2]], c(1, 2), function(num) ifelse(abs(num) <= 1e-10, 0, num))

  sol_path[[3]] <- apply(sol_path[[2]], 2, function(cc) as.numeric(solve(LL) %*% cc))
  sol_path[[3]] <- apply(sol_path[[3]], c(1, 2), function(num) ifelse(abs(num) <= 1e-10, 0, num))

  names(sol_path) <- c("lambda", "gamma", "beta")

  if (is.null(lambda)) {
    return(sol_path)
  } else if (is.numeric(lambda) & all(lambda >= 0)) { ## if a vector of lambda is given, return the parameter values at given lambdas
    res <- vector(mode = "list", length = 3)
    res[[1]] <- lambda
    res[[2]] <- res[[3]] <- matrix(NA, ncol = length(lambda), nrow = length(beta_tilde))
    for (i in 1:length(lambda)) {
      ll <- lambda[i] ## current tuning parameter

      if (ll >= max(sol_path$lambda)) {
        res[[2]][, i] <- res[[3]][, i] <- 0
      } else {
        lower_index <- min(which(sol_path$lambda <= ll))
        upper_index <- max(which(sol_path$lambda > ll))

        lower_lambda <- sol_path$lambda[lower_index]
        upper_lambda <- sol_path$lambda[upper_index]

        res[[2]][, i] <- (upper_lambda - ll) / (upper_lambda - lower_lambda) * sol_path[[2]][, lower_index] +
          (ll - lower_lambda) / (upper_lambda - lower_lambda) * sol_path[[2]][, upper_index]

        res[[3]][, i] <- (upper_lambda - ll) / (upper_lambda - lower_lambda) * sol_path[[3]][, lower_index] +
          (ll - lower_lambda) / (upper_lambda - lower_lambda) * sol_path[[3]][, upper_index]
      }
    }
    names(res) <- c("lambda", "gamma", "beta")
    colnames(res[[2]]) <- colnames(res[[3]]) <- round(lambda, 2)
    return(res)
  } else {
    stop("lambda has to be either null or a numeric vector with non-negative entries.")
  }

}



# Obtain the trajectory of LSA-LASSO via path algorithm
rlasso_path <- function(X, Y, lambda = NULL) {
  # dimension of the parameter
  dd <- ncol(X)
  nn <- nrow(X)
  # OLS estimator
  beta_tilde <- as.numeric(solve(t(X) %*% X) %*% t(X) %*% Y)

  # residuals
  ei <- as.numeric(Y - X %*% beta_tilde)

  # sandwich covariance matrix for OLS estimator
  Sigma_hat <- solve(t(X) %*% X) %*% (t(X) %*% diag(ei ^ 2) %*% X) %*%
    solve(t(X) %*% X)


  Omega <- solve(Sigma_hat)
  LL <- chol(Omega) # cholesky decomposition of the weight matrix

  gamma_tilde <- LL %*% beta_tilde


  # defining the objective function
  sol_results <- genlasso::genlasso(y = gamma_tilde, D = solve(LL), approx = TRUE, svd = TRUE, rtol = 1e-15, eps = 1e-15)

  ## summarize the solution path at the knots
  sol_path <- vector("list", 3)
  sol_path[[1]] <- c(2 * sol_results$lambda, 0)

  sol_path[[2]] <- cbind(sol_results$beta, sol_results$bls)
  colnames(sol_path[[2]])[ncol(sol_path[[2]])] <- "0"
  colnames(sol_path[[2]]) <- round(2 * as.numeric(colnames(sol_path[[2]])), 2)
  sol_path[[2]] <- apply(sol_path[[2]], c(1, 2), function(num) ifelse(abs(num) <= 1e-10, 0, num))

  sol_path[[3]] <- apply(sol_path[[2]], 2, function(cc) as.numeric(solve(LL) %*% cc))
  sol_path[[3]] <- apply(sol_path[[3]], c(1, 2), function(num) ifelse(abs(num) <= 1e-10, 0, num))

  names(sol_path) <- c("lambda", "gamma", "beta")

  if (is.null(lambda)) {
    return(sol_path)
  } else if (is.numeric(lambda) & all(lambda >= 0)) { ## if a vector of lambda is given, return the parameter values at given lambdas
    res <- vector(mode = "list", length = 3)
    res[[1]] <- lambda
    res[[2]] <- res[[3]] <- matrix(NA, ncol = length(lambda), nrow = length(beta_tilde))
    for (i in 1:length(lambda)) {
      ll <- lambda[i] ## current tuning parameter

      if (ll >= max(sol_path$lambda)) {
        res[[2]][, i] <- res[[3]][, i] <- 0
      } else {
        lower_index <- min(which(sol_path$lambda <= ll))
        upper_index <- max(which(sol_path$lambda > ll))

        lower_lambda <- sol_path$lambda[lower_index]
        upper_lambda <- sol_path$lambda[upper_index]

        res[[2]][, i] <- (upper_lambda - ll) / (upper_lambda - lower_lambda) * sol_path[[2]][, lower_index] +
          (ll - lower_lambda) / (upper_lambda - lower_lambda) * sol_path[[2]][, upper_index]

        res[[3]][, i] <- (upper_lambda - ll) / (upper_lambda - lower_lambda) * sol_path[[3]][, lower_index] +
          (ll - lower_lambda) / (upper_lambda - lower_lambda) * sol_path[[3]][, upper_index]
      }
    }
    names(res) <- c("lambda", "gamma", "beta")
    colnames(res[[2]]) <- colnames(res[[3]]) <- round(lambda, 2)
    return(res)
  } else {
    stop("lambda has to be either null or a numeric vector with non-negative entries.")
  }

}



#' Given the sample size and the number of folds, make folds for cross-validation
#'
#' @export

make_folds <- function(n, nfolds) {
  nn <- round(n / nfolds)
  sizes <- rep(nn, nfolds)
  sizes[nfolds] <- sizes[nfolds] + n - nn * nfolds
  b <- c(0, cumsum(sizes))
  ii <- sample(n)
  folds <- list()
  for (i in seq(nfolds))
    folds[[i]] <- ii[seq(b[i] + 1, b[i + 1])]
  folds
}




#' LASSO using path algorithm
#'
#' @param X design matrix
#' @param Y outcome vector
#' @param lambda a vector of tuning parameters. If left unspecified, using the
#' sequence of numbers from 0 to the smallest tuning parameter such that the
#' solution is the zero vector
#' @param nlambda the number of the sequence of the tuning parameters.
#'
#' @export

lasso <- function(X, Y, lambda = NULL, nlambda = 100) {
  ## obtain the solution path of LASSO
  pp <- ncol(X) ## number of the covariates
  nn <- nrow(X) ## number of observations
  if (nn <= pp) stop("This version of LASSO can only solve problems with n > p!")
  lasso_sol_path <- lasso_path(Y = Y, X = X)

  if (is.null(lambda)) lambda <- seq(0, lasso_sol_path$lambda[1], len = nlambda)

  nlam <- length(lambda)

  RES <- list(lambda = lambda,
              beta = vector(mode = "list", length = nlam))  ## make the list

  for (ii in 1:nlam) {
    lam <- lambda[ii]


    if (lam >= lasso_sol_path$lambda[1]) {
      RES$beta[[ii]] <- rep(0, pp)
    } else if (lam == 0){
      RES$beta[[ii]] <- lasso_sol_path$beta[, ncol(lasso_sol_path$beta)]
    } else {
      ## beta with the given lambda is the convex optimization of the two closest knots
      lam_up_id <- max(which(lasso_sol_path$lambda >= lam))
      lam_down_id <- min(which(lasso_sol_path$lambda < lam))
      lam_up <- lasso_sol_path$lambda[lam_up_id]
      lam_down <- lasso_sol_path$lambda[lam_down_id]
      beta_up <- lasso_sol_path$beta[, lam_up_id]
      beta_down <- lasso_sol_path$beta[, lam_down_id]
      RES$beta[[ii]] <- (lam - lam_down) / (lam_up - lam_down) * beta_up +
        (lam_up - lam) / (lam_up - lam_down) * beta_down
    }
  }
  return(RES)
}

#' LASSO using path algorithm, where an optimal tuning parameter is found by
#' cross validation
#'
#' @param X design matrix
#' @param Y outcome vector
#' @param lambda a vector of tuning parameters. If left unspecified, using the
#' sequence of numbers from 0 to the smallest tuning parameter such that the
#' solution is the zero vector
#' @param nlambda the number of the sequence of the tuning parameters.
#' @param nfolds number of folds in the cross-validation
#'
#' @export


cv.lasso <- function(X, Y, lambda = NULL, nlambda = 100, nfolds = 5){
  n <- nrow(X) ## sample size
  folds <- make_folds(n, nfolds)

  ## obtain the lambda vector
  if (is.null(lambda)) {
    lambda <- seq(0, lasso_path(X = X, Y = Y)$lambda[1], length.out = nlambda)
  }

  cv_err_mat <- matrix(nrow = nlambda, ncol = nfolds) ## residual sum of square



  for (jj in 1:nfolds) {
    ## get the training data and testing data
    testing_id <- folds[[jj]]
    training_id <- (1:n)[-testing_id]
    training_X <- X[training_id,]
    training_Y <- Y[training_id]
    testing_X <- X[testing_id,]
    testing_Y <- Y[testing_id]
    ## compute the rLASSO estimate using the training data
    fold_RES <- lasso(X = training_X, Y = training_Y, lambda = lambda)

    for (kk in 1:nlambda) {
      cv_err_mat[kk, jj] <- mean((testing_Y - testing_X %*% fold_RES$beta[[kk]]) ^ 2)
    }
  }
  cv_err <- rowMeans(cv_err_mat)

  opt_lambda <- lambda[which.min(cv_err)]

  opt_beta <- lasso(X = X, Y = Y, lambda = opt_lambda)$beta[[1]]

  RES <- list(lambda = lambda, cv_err = cv_err, opt_lambda = opt_lambda, opt_beta = opt_beta)
}

#' rLASSO using path algorithm
#'
#' @param X design matrix
#' @param Y outcome vector
#' @param lambda a vector of tuning parameters. If left unspecified, using the
#' sequence of numbers from 0 to the smallest tuning parameter such that the
#' solution is the zero vector
#' @param nlambda the number of the sequence of the tuning parameters.
#'
#' @export


rlasso <- function(X, Y, lambda = NULL, nlambda = 100) {
  ## obtain the solution path of rLASSO
  pp <- ncol(X) ## number of the covariates
  nn <- nrow(X) ## number of observations
  if (nn <= pp) stop("rLASSO can only solve problems with n > p!")
  rlasso_sol_path <- rlasso_path(Y = Y, X = X)

  if (is.null(lambda)) lambda <- seq(0, rlasso_sol_path$lambda[1], len = nlambda)

  nlam <- length(lambda)

  RES <- list(lambda = lambda,
              beta = vector(mode = "list", length = nlam)) ## make the list

  for (ii in 1:nlam) {
    lam <- lambda[ii]


    if (lam >= rlasso_sol_path$lambda[1]) {
      RES$beta[[ii]] <- rep(0, pp)
    } else if (lam == 0){
      RES$beta[[ii]] <- rlasso_sol_path$beta[, ncol(rlasso_sol_path$beta)]
    } else {
      ## beta with the given lambda is the convex optimization of the two closest knots
      lam_up_id <- max(which(rlasso_sol_path$lambda >= lam))
      lam_down_id <- min(which(rlasso_sol_path$lambda < lam))
      lam_up <- rlasso_sol_path$lambda[lam_up_id]
      lam_down <- rlasso_sol_path$lambda[lam_down_id]
      beta_up <- rlasso_sol_path$beta[, lam_up_id]
      beta_down <- rlasso_sol_path$beta[, lam_down_id]
      RES$beta[[ii]] <- (lam - lam_down) / (lam_up - lam_down) * beta_up +
        (lam_up - lam) / (lam_up - lam_down) * beta_down
    }
  }
  return(RES)
}


#' rLASSO using path algorithm, where an optimal tuning parameter is found by
#' cross validation
#'
#' @param X design matrix
#' @param Y outcome vector
#' @param lambda a vector of tuning parameters. If left unspecified, using the
#' sequence of numbers from 0 to the smallest tuning parameter such that the
#' solution is the zero vector
#' @param nlambda the number of the sequence of the tuning parameters.
#' @param nfolds number of folds in the cross-validation
#'
#' @export


cv.rlasso <- function(X, Y, lambda = NULL, nlambda = 100, nfolds = 5){
  n <- nrow(X) ## sample size
  folds <- make_folds(n, nfolds)

  ## obtain the lambda vector
  if (is.null(lambda)) {
    lambda <- seq(0, rlasso_path(X = X, Y = Y)$lambda[1], length.out = nlambda)
  }

  cv_err_mat <- matrix(nrow = nlambda, ncol = nfolds) ## residual sum of square



  for (jj in 1:nfolds) {
    ## get the training data and testing data
    testing_id <- folds[[jj]]
    training_id <- (1:n)[-testing_id]
    training_X <- X[training_id,]
    training_Y <- Y[training_id]
    testing_X <- X[testing_id,]
    testing_Y <- Y[testing_id]
    ## compute the rLASSO estimate using the training data
    fold_RES <- rlasso(X = training_X, Y = training_Y, lambda = lambda)

    for (kk in 1:nlambda) {
      cv_err_mat[kk, jj] <- mean((testing_Y - testing_X %*% fold_RES$beta[[kk]]) ^ 2)
    }
  }
  cv_err <- rowMeans(cv_err_mat)

  opt_lambda <- lambda[which.min(cv_err)]

  opt_beta <- rlasso(X = X, Y = Y, lambda = opt_lambda)$beta[[1]]

  RES <- list(lambda = lambda, cv_err = cv_err, opt_lambda = opt_lambda, opt_beta = opt_beta)
}


