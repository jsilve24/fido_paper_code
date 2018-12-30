require(mongrel)
source("src/fit_methods.R")
source("src/dataset_methods.R")
source("src/utils.R")

# main  function ---------------------------------------------------

#' Fit conjugate linear model
#' 
#' @param mdataset an mdataset object
#' @param n_samples number of samples to draw
#' @param decomposition which parameterization to use 
#'   ("eigen":default, "cholesky")
#' @param ret_mongrelfit return mongrelfit instead of mfit object (default:FALSE)
#' @param ... other parameters passed to function mongrel
#' @return mfit object (returns mongrelfit if ret_mongrelfit)
fit_CLM <- function(mdataset, n_samples=2000, ...) {
  eta.hat <- t(driver::alr(t(mdataset$Y+0.65)))
  start_time <- Sys.time()
  fit <- conjugateLinearModel(eta.hat, mdataset$X, mdataset$Theta, mdataset$Gamma, mdataset$Xi, mdataset$upsilon, n_samples=n_samples)
  end_time <- Sys.time()

  total_runtime <- end_time - start_time

  lambda_RMSE <- get_Lambda_RMSE(mdataset$Lambda_true, fit$Lambda)

  outside_percent <- get_95CI(mdataset$Lambda_true, fit$Lambda)

  metadata <- metadata(0, total_runtime, n_samples, lambda_RMSE, outside_percent)

  m <- mfit(N=mdataset$N, D=mdataset$D, Q=mdataset$Q, iter=as.integer(n_samples), 
            Lambda=fit$Lambda, 
            Sigma=fit$Sigma,
            mdataset=mdataset,
            metadata=metadata)
  # no Eta
  return(m)
}

