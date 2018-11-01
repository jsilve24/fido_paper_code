require(mongrel)
source("src/fit_methods.R")
source("src/dataset_methods.R")
source("src/utils.R")

# main  function ---------------------------------------------------

#' Fit multinomial model using mongrel
#' 
#' WARNING: Req that your current working directory is main folder
#'  "mongrel_paper_code"
#'  
#' WARNING: Currently uses a random seed... 
#' 
#' @param mdataset an mdataset object
#' @param n_samples number of samples to draw
#' @param decomposition which parameterization to use 
#'   ("eigen":default, "cholesky")
#' @param ret_mongrelfit return mongrelfit instead of mfit object (default:FALSE)
#' @param ... other parameters passed to function mongrel
#' @return mfit object (returns mongrelfit if ret_mongrelfit)
fit_mongrel <- function(mdataset, n_samples=2000, decomposition="eigen", ret_mongrelfit=FALSE, ...){

  start_time <- Sys.time()
  cat("Start time:",start_time,"\n")
  fit <- mongrel(Y=mdataset$Y, X=mdataset$X, upsilon=mdataset$upsilon, 
                 Theta=mdataset$Theta, Gamma=mdataset$Gamma, Xi=mdataset$Xi, 
                 init=random_mongrel_init(mdataset$Y), 
                 decomp_method=decomposition, n_samples=n_samples,  ...)
  end_time <- Sys.time()
  cat("End time:",end_time,"\n")
  if (ret_mongrelfit) return(fit)

  total_runtime <- end_time - start_time
  cat("Total runtime:",total_runtime,"\n")

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
  m$Eta <- fit$Eta
  return(m)
}

