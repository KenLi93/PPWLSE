#' Rotated James-Stein Estimator
#' @export
rJS <- function(X, Y, lambda = NULL, nlambda = 100) {
  ## obtain the solution path of rLASSO
  pp <- ncol(X) ## number of the covariates
  nn <- nrow(X) ## number of observations


  if (nn <= pp) stop("rRidge can only solve problems with n > p!")

  ## compute the OLS estimator and sandwich variance
  beta_tilde <- as.numeric(solve(t(X) %*% X) %*% t(X) %*% Y)
  Sigma_n <- solve(t(X) %*% X) %*% t(X) %*% diag(as.numeric((Y - X %*% beta_tilde)^2)) %*% X %*%
    solve(t(X) %*% X)

  if (is.null(lambda)) {
    lambda <- seq(0, sqrt(nn), length.out = nlambda)
  }

  nlam <- length(lambda)
  RES <- list(lambda = lambda,
              beta = vector(mode = "list", length = nlam)) ## make the list

  for (ii in 1:nlam) {
    lam <- lambda[ii]

    RES$beta[[ii]] <- solve(diag(nrow = pp, ncol = pp) + lam * Sigma_n %*% t(X) %*% (X)) %*% beta_tilde
  }
  return(RES)
}


#' Scaled OLS estimator
#' @export
scaledLS <- function(X, Y, lambda = NULL, nlambda = 100) {
  ## obtain the solution path of rLASSO
  pp <- ncol(X) ## number of the covariates
  nn <- nrow(X) ## number of observations


  if (nn <= pp) stop("rRidge can only solve problems with n > p!")

  ## compute the OLS estimator and sandwich variance
  beta_tilde <- as.numeric(solve(t(X) %*% X) %*% t(X) %*% Y)


  if (is.null(lambda)) {
    lambda <- seq(0, sqrt(nn), length.out = nlambda)
  }

  nlam <- length(lambda)
  RES <- list(lambda = lambda,
              beta = vector(mode = "list", length = nlam)) ## make the list

  for (ii in 1:nlam) {
    lam <- lambda[ii]

    RES$beta[[ii]] <- beta_tilde / (1 + lam)
  }
  return(RES)
}


#' James-Stein estimator for linear regression
#'
#' @export
JS <- function(X, Y) {
  ## obtain the solution path of rLASSO
  pp <- ncol(X) ## number of the covariates
  nn <- nrow(X) ## number of observations


  if (nn <= pp) stop("James-Stein estimator can only solve problems with n > p!")
  if (pp < 2) stop("James-Stein estimator must have pp >= 2")
  ## compute the OLS estimator and sandwich variance
  beta_tilde <- as.numeric(solve(t(X) %*% X) %*% t(X) %*% Y)
  sigma2 <- sum((Y - X %*% beta_tilde) ^ 2) / (nn - pp)

  beta_JS <- as.numeric((1 - (nn - pp) * (pp - 2) / ((nn - pp + 2) * beta_tilde %*% t(X) %*%
                                            X %*% beta_tilde))) * beta_tilde
  return(beta_JS)
}



#' Positive-Part James-Stein estimator for linear regression
#'
#' @export
PJS <- function(X, Y) {
  ## obtain the solution path of rLASSO
  pp <- ncol(X) ## number of the covariates
  nn <- nrow(X) ## number of observations


  if (nn <= pp) stop("James-Stein estimator can only solve problems with n > p!")

  ## compute the OLS estimator and sandwich variance
  beta_tilde <- as.numeric(solve(t(X) %*% X) %*% t(X) %*% Y)
  sigma2 <- sum((Y - X %*% beta_tilde) ^ 2) / (nn - pp)

  beta_PJS <- max(as.numeric((1 - (nn - pp) * (pp - 2) / ((nn - pp + 2) * beta_tilde %*% t(X) %*%
                                                       X %*% beta_tilde))), 0) * beta_tilde
  return(beta_PJS)
}
