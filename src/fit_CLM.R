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

  # get Lambda MSE
  ref <- c(mdataset$Lambda_true)
  est_Lambda <- fit$Lambda
  dim(est_Lambda) <- c(nrow(fit$Lambda)*ncol(fit$Lambda),n_samples)
  lambda_MSE <- mean(apply(est_Lambda, 2, sample_SE, y=ref))

  # get count of Lambdas outside 95% CI
  intervals <- apply(fit$Lambda, c(1,2), function(x) quantile(x, probs=c(0.025, 0.975)))
  lower <- c(intervals[1,,])
  upper <- c(intervals[2,,])

  outside_95CI <- sum(apply(rbind(lower, upper, ref), 2, lambda_outside_bounds))/((mdataset$D-1)*mdataset$Q)

  metadata <- metadata(0, total_runtime, n_samples, lambda_MSE, outside_95CI)

  m <- mfit(N=mdataset$N, D=mdataset$D, Q=mdataset$Q, iter=as.integer(n_samples), 
            Lambda=fit$Lambda, 
            Sigma=fit$Sigma,
            mdataset=mdataset,
            metadata=metadata)
  # no Eta
  return(m)
}

