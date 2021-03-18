
#' Rotated Ridge Estimator
#' @export
rridge <- function(X, Y, lambda = NULL, nlambda = 100) {
  ## obtain the solution path of rLASSO
  pp <- ncol(X) ## number of the covariates
  nn <- nrow(X) ## number of observations


  if (nn <= pp) stop("rRidge can only solve problems with n > p!")

  ## compute the OLS estimator and sandwich variance
  beta_tilde <- as.numeric(solve(t(X) %*% X) %*% t(X) %*% Y)
  Sigma_n <- solve(t(X) %*% X) %*% t(X) %*% diag(as.numeric((Y - X %*% beta_tilde)^2)) %*% X %*%
    solve(t(X) %*% X)

  if (is.null(lambda)) {
    ## find the smallest lambda such that max|\beta| < eps
    lambda <- seq(0, sqrt(nn), length.out = nlambda)
  }

  nlam <- length(lambda)
  RES <- list(lambda = lambda,
              beta = vector(mode = "list", length = nlam)) ## make the list

  for (ii in 1:nlam) {
    lam <- lambda[ii]

    RES$beta[[ii]] <- solve(diag(nrow = pp, ncol = pp) + lam * Sigma_n) %*% beta_tilde
  }
  return(RES)
}




#' Ridge Estimator
#'
#' @export
ridge <- function(X, Y, lambda = NULL, nlambda = 100) {
  ## obtain the solution path of rLASSO
  pp <- ncol(X) ## number of the covariates
  nn <- nrow(X) ## number of observations


  if (nn <= pp) stop("This version of Ridge can only solve problems with n > p!")

  ## compute the OLS estimator and sandwich variance
  beta_tilde <- as.numeric(solve(t(X) %*% X) %*% t(X) %*% Y)
  Sigma_n <- solve(t(X) %*% X)

  if (is.null(lambda)) {
    ## find the smallest lambda such that max|\beta| < eps
    lambda <- seq(0, sqrt(nn), length.out = nlambda)
  }

  nlam <- length(lambda)
  RES <- list(lambda = lambda,
              beta = vector(mode = "list", length = nlam)) ## make the list

  for (ii in 1:nlam) {
    lam <- lambda[ii]

    RES$beta[[ii]] <- solve(diag(nrow = pp, ncol = pp) + lam * Sigma_n) %*% beta_tilde
  }
  return(RES)
}
