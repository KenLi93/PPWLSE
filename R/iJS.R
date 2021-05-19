#' Imroved James-Stein estimator under heterogeneity
#' @export
iJS <- function(X, Y) {
  Y <- as.numeric(Y)
  X <- as.matrix(X)

  ## obtain the solution path of rLASSO
  pp <- ncol(X) ## number of the covariates
  nn <- nrow(X) ## number of observations


  if (nn <= pp) stop("James-Stein estimator can only solve problems with n > p!")

  ## compute the OLS estimator and sandwich variance
  beta_tilde <- as.numeric(solve(t(X) %*% X) %*% t(X) %*% Y)

  ## shinkage factor
  c_hat <- 1 - sum(diag(X %*% solve(t(X) %*% X) %*% t(X) %*%
                          diag(as.numeric((Y - X %*% beta_tilde) ^ 2)))) /
    as.numeric(t(beta_tilde) %*% t(X) %*% X %*% beta_tilde)

  beta_iJS <- as.numeric(c_hat * beta_tilde)

  return(beta_iJS)
}

#' Improved James-Stein positive-part
#' @export
iJSP <- function(X, Y) {
  Y <- as.numeric(Y)
  X <- as.matrix(X)

  ## obtain the solution path of rLASSO
  pp <- ncol(X) ## number of the covariates
  nn <- nrow(X) ## number of observations


  if (nn <= pp) stop("James-Stein estimator can only solve problems with n > p!")

  ## compute the OLS estimator and sandwich variance
  beta_tilde <- as.numeric(solve(t(X) %*% X) %*% t(X) %*% Y)

  ## shinkage factor
  c_hat <- 1 - sum(diag(X %*% solve(t(X) %*% X) %*% t(X) %*%
                          diag(as.numeric((Y - X %*% beta_tilde) ^ 2)))) /
    as.numeric(t(beta_tilde) %*% t(X) %*% X %*% beta_tilde)

  beta_iJSP <- as.numeric(max(c_hat, 0) * beta_tilde)

  return(beta_iJSP)
}

#' Improved Ridge regressions
#' @export
iridge <- function(X, Y, eps = 0.001) {
  Y <- as.numeric(Y)
  X <- as.matrix(X)

  ## obtain the solution path of rLASSO
  pp <- ncol(X) ## number of the covariates
  nn <- nrow(X) ## number of observations


  if (nn <= pp) stop("James-Stein estimator can only solve problems with n > p!")

  ## compute the OLS estimator and sandwich variance
  beta_tilde <- as.numeric(solve(t(X) %*% X) %*% t(X) %*% Y)


  WW <- solve(t(X) %*% X) ## covariance matrix of X
  II <- diag(nrow = pp) ## identify matrix of dimension p
  ee <- diag(as.numeric((Y - X %*% beta_tilde) ^ 2))
  pred_risk <- function(lambda) {
    beta_tilde %*% solve(II + lambda * WW) %*% t(X) %*% X %*%
      solve(II + lambda * WW) %*% beta_tilde -
      2 * beta_tilde %*% t(X) %*% X %*% solve(II + lambda * WW) %*%
      beta_tilde +
      2 * sum(diag(X %*% solve(II + lambda * WW) %*% solve(t(X) %*% X) %*%
                     t(X) %*% ee))
  }
  maxv <- max(eigen(-solve(WW))$values)
  opt_lam <- optim(0, pred_risk,
                   lower = min(maxv + eps, 0),
                   upper = nn,
                   method = "L-BFGS-B")$par
  beta_hat <- solve(II + opt_lam * WW) %*% beta_tilde
  return(beta_hat)
}

#' Improved Ridge regression positive part
#' @export
iridgep <- function(X, Y, eps = 0.001) {
  Y <- as.numeric(Y)
  X <- as.matrix(X)

  ## obtain the solution path of rLASSO
  pp <- ncol(X) ## number of the covariates
  nn <- nrow(X) ## number of observations


  if (nn <= pp) stop("James-Stein estimator can only solve problems with n > p!")

  ## compute the OLS estimator and sandwich variance
  beta_tilde <- as.numeric(solve(t(X) %*% X) %*% t(X) %*% Y)


  WW <- solve(t(X) %*% X) ## covariance matrix of X
  II <- diag(nrow = pp) ## identify matrix of dimension p
  ee <- diag(as.numeric((Y - X %*% beta_tilde) ^ 2))
  pred_risk <- function(lambda) {
    beta_tilde %*% solve(II + lambda * WW) %*% t(X) %*% X %*%
      solve(II + lambda * WW) %*% beta_tilde -
      2 * beta_tilde %*% t(X) %*% X %*% solve(II + lambda * WW) %*%
      beta_tilde +
      2 * sum(diag(X %*% solve(II + lambda * WW) %*% solve(t(X) %*% X) %*%
                     t(X) %*% ee))
  }
  maxv <- max(eigen(-solve(WW))$values)
  opt_lam <- max(optim(0, pred_risk,
                   lower = min(maxv + eps, 0),
                   upper = nn,
                   method = "L-BFGS-B")$par, 0)
  beta_hat <- solve(II + opt_lam * WW) %*% beta_tilde
  return(beta_hat)
}




#' ridge regression using generalized cross-validation
#' @export
ridge_gcv <- function(X, Y, eps = 0.001) {
  Y <- as.numeric(Y)
  X <- as.matrix(X)

  ## obtain the solution path of rLASSO
  pp <- ncol(X) ## number of the covariates
  nn <- nrow(X) ## number of observations


  if (nn <= pp) stop("James-Stein estimator can only solve problems with n > p!")

  ## compute the OLS estimator and sandwich variance
  beta_tilde <- as.numeric(solve(t(X) %*% X) %*% t(X) %*% Y)


  WW <- solve(t(X) %*% X) ## covariance matrix of X
  II <- diag(nrow = pp) ## identify matrix of dimension p
  In <- diag(nrow = nn)
  ee <- diag(as.numeric((Y - X %*% beta_tilde) ^ 2))
  pred_risk <- function(lambda) {
    nn * sum((Y - X %*% solve(t(X) %*% X + lambda * II) %*% t(X) %*% Y) ^ 2) /
      (sum(diag(In - X %*% solve(t(X) %*% X + lambda * II) %*% t(X)))) ^ 2
  }
  maxv <- max(eigen(-solve(WW))$values)
  opt_lam <- optim(0, pred_risk,
                       lower = min(maxv + eps, 0),
                       upper = nn,
                       method = "L-BFGS-B")$par
  beta_hat <- solve(II + opt_lam * WW) %*% beta_tilde
  return(beta_hat)
}





