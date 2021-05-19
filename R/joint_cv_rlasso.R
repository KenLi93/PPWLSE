#' Joint Cross-Validation
#' @export
cv_joint_shrinkage_est <- function(X, Y, lambda = NULL, nlambda = 100, nfolds = 5, nrep = 200){
  nn <- nrow(X) ## sample size


  ## obtain the lambda vector
  if (is.null(lambda)) {

    ## find the smallest lambda such that max|\beta| < eps
    rridge_lambda <- ridge_lambda <- rJS_lambda <- seq(0, nn, length.out = nlambda)
    lasso_lambda <- rlasso_lambda <-
     seq(0, max(rlasso_path(X = X, Y = Y)$lambda[1], lasso_path(X = X, Y = Y)$lambda[1]), length.out = nlambda)
  } else {
    rridge_lambda <- ridge_lambda <- rJS_lambda <- lasso_lambda <- rlasso_lambda <- lambda
  }

  rridge_cv_err_mat_full <- ridge_cv_err_mat_full <- rJS_cv_err_mat_full <-
    lasso_cv_err_mat_full <- rlasso_cv_err_mat_full <- matrix(NA, nrow = nlambda, ncol = nrep)

  for (ll in 1:nrep){
    rridge_cv_err_mat <- ridge_cv_err_mat <- rJS_cv_err_mat <-
      lasso_cv_err_mat <- rlasso_cv_err_mat <- matrix(nrow = nlambda, ncol = nfolds) ## residual sum of square

    folds <- make_folds(nn, nfolds)

    for (jj in 1:nfolds) {
      ## get the training data and testing data
      testing_id <- folds[[jj]]
      training_id <- (1:nn)[-testing_id]
      training_X <- X[training_id,]
      training_Y <- Y[training_id]
      testing_X <- X[testing_id,]
      testing_Y <- Y[testing_id]



      ## compute the rLASSO estimate using the training data
      fold_rridge_RES <- rridge(X = training_X, Y = training_Y, lambda = rridge_lambda)
      fold_ridge_RES <- ridge(X = training_X, Y = training_Y, lambda = ridge_lambda)
      fold_rJS_RES <- rJS(X = training_X, Y = training_Y, lambda = rJS_lambda)
      fold_rlasso_RES <- rlasso(X = training_X, Y = training_Y, lambda = lasso_lambda)
      fold_lasso_RES <- lasso(X = training_X, Y = training_Y, lambda = rlasso_lambda)

      for (kk in 1:nlambda) {
        ridge_cv_err_mat[kk, jj] <- mean((testing_Y - testing_X %*% fold_ridge_RES$beta[[kk]]) ^ 2)
        rridge_cv_err_mat[kk, jj] <- mean((testing_Y - testing_X %*% fold_rridge_RES$beta[[kk]]) ^ 2)
        rJS_cv_err_mat[kk, jj] <- mean((testing_Y - testing_X %*% fold_rJS_RES$beta[[kk]]) ^ 2)
        lasso_cv_err_mat[kk, jj] <- mean((testing_Y - testing_X %*% fold_lasso_RES$beta[[kk]]) ^ 2)
        rlasso_cv_err_mat[kk, jj] <- mean((testing_Y - testing_X %*% fold_rlasso_RES$beta[[kk]]) ^ 2)

      }

    }
    ridge_cv_err_mat_full[, ll] <- rowMeans(ridge_cv_err_mat)
    rridge_cv_err_mat_full[, ll] <- rowMeans(rridge_cv_err_mat)
    rJS_cv_err_mat_full[, ll] <- rowMeans(rJS_cv_err_mat)
    lasso_cv_err_mat_full[, ll] <- rowMeans(lasso_cv_err_mat)
    rlasso_cv_err_mat_full[, ll] <- rowMeans(rlasso_cv_err_mat)
  }
  ridge_cv_err <- rowMeans(ridge_cv_err_mat_full)
  rridge_cv_err <- rowMeans(rridge_cv_err_mat_full)
  rJS_cv_err <- rowMeans(rJS_cv_err_mat_full)
  lasso_cv_err <- rowMeans(lasso_cv_err_mat_full)
  rlasso_cv_err <- rowMeans(rlasso_cv_err_mat_full)

  ridge_opt_lambda <- ridge_lambda[which.min(ridge_cv_err)]
  rridge_opt_lambda <- rridge_lambda[which.min(rridge_cv_err)]
  rJS_opt_lambda <- rJS_lambda[which.min(rJS_cv_err)]
  lasso_opt_lambda <- lasso_lambda[which.min(lasso_cv_err)]
  rlasso_opt_lambda <- rlasso_lambda[which.min(rlasso_cv_err)]


  ridge_opt_beta <- ridge(X = X, Y = Y, lambda = ridge_opt_lambda)$beta[[1]]
  rridge_opt_beta <- rridge(X = X, Y = Y, lambda = rridge_opt_lambda)$beta[[1]]
  rJS_opt_beta <- rJS(X = X, Y = Y, lambda = rJS_opt_lambda)$beta[[1]]
  lasso_opt_beta <- lasso(X = X, Y = Y, lambda = lasso_opt_lambda)$beta[[1]]
  rlasso_opt_beta <- rlasso(X = X, Y = Y, lambda = rlasso_opt_lambda)$beta[[1]]

  RES <- list(
    ridge = list(lambda = ridge_lambda, cv_err = ridge_cv_err, opt_lambda = ridge_opt_lambda, opt_beta = ridge_opt_beta),
    rridge = list(lambda = rridge_lambda, cv_err = rridge_cv_err, opt_lambda = rridge_opt_lambda, opt_beta = rridge_opt_beta),
    rJS = list(lambda = rJS_lambda, cv_err = rJS_cv_err, opt_lambda = rJS_opt_lambda, opt_beta = rJS_opt_beta),
    lasso = list(lambda = lasso_lambda, cv_err = lasso_cv_err, opt_lambda = lasso_opt_lambda, opt_beta = lasso_opt_beta),
    rlasso = list(lambda = rlasso_lambda, cv_err = rlasso_cv_err, opt_lambda = rlasso_opt_lambda, opt_beta = rlasso_opt_beta)
  )
  RES
}

#' Joint Cross-Validation of Ridge and James-Stein
#' @export
## cv for ridge and JS
cv_ridge_JS <- function(X, Y, lambda = NULL, nlambda = 100, nfolds = 10, nrep = 200){
  nn <- nrow(X) ## sample size


  ## obtain the lambda vector
  if (is.null(lambda)) {

    ## find the smallest lambda such that max|\beta| < eps
    lambda <- seq(0, nn, length.out = nlambda)
  }

  ridge_cv_err_mat_full <- scaledLS_cv_err_mat_full <- matrix(NA, nrow = nlambda, ncol = nrep)

  for (ll in 1:nrep){
    ridge_cv_err_mat <- scaledLS_cv_err_mat <- matrix(nrow = nlambda, ncol = nfolds) ## residual sum of square

    folds <- make_folds(nn, nfolds)

    for (jj in 1:nfolds) {
      ## get the training data and testing data
      testing_id <- folds[[jj]]
      training_id <- (1:nn)[-testing_id]
      training_X <- X[training_id,]
      training_Y <- Y[training_id]
      testing_X <- X[testing_id,]
      testing_Y <- Y[testing_id]



      ## compute the rLASSO estimate using the training data
      fold_ridge_RES <- ridge(X = training_X, Y = training_Y, lambda = lambda)
      fold_scaledLS_RES <- scaledLS(X = training_X, Y = training_Y, lambda = lambda)


      for (kk in 1:nlambda) {
        ridge_cv_err_mat[kk, jj] <- mean((testing_Y - testing_X %*% fold_ridge_RES$beta[[kk]]) ^ 2)
        scaledLS_cv_err_mat[kk, jj] <- mean((testing_Y - testing_X %*% fold_scaledLS_RES$beta[[kk]]) ^ 2)
      }

    }
    ridge_cv_err_mat_full[, ll] <- rowMeans(ridge_cv_err_mat)
    scaledLS_cv_err_mat_full[, ll] <- rowMeans(scaledLS_cv_err_mat)
  }
  ridge_cv_err <- rowMeans(ridge_cv_err_mat_full)
  scaledLS_cv_err <- rowMeans(scaledLS_cv_err_mat_full)

  ridge_opt_lambda <- lambda[which.min(ridge_cv_err)]
  scaledLS_opt_lambda <- lambda[which.min(scaledLS_cv_err)]


  ridge_opt_beta <- ridge(X = X, Y = Y, lambda = ridge_opt_lambda)$beta[[1]]
  scaledLS_opt_beta <- scaledLS(X = X, Y = Y, lambda = scaledLS_opt_lambda)$beta[[1]]

  RES <- list(
    ridge = list(lambda = lambda, cv_err = ridge_cv_err, opt_lambda = ridge_opt_lambda, opt_beta = ridge_opt_beta),
    scaledLS = list(lambda = lambda, cv_err = scaledLS_cv_err, opt_lambda = scaledLS_opt_lambda, opt_beta = scaledLS_opt_beta)
  )
  RES
}



