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
  fit <- mongrel(Y=mdataset$Y, X=mdataset$X, upsilon=mdataset$upsilon, 
                 Theta=mdataset$Theta, Gamma=mdataset$Gamma, Xi=mdataset$Xi, 
                 init=random_mongrel_init(mdataset$Y), 
                 decomp_method=decomposition, n_samples=n_samples, ...)
  end_time <- Sys.time()
  if (ret_mongrelfit) return(fit)

  total_runtime <- end_time - start_time

  lambda_RMSE <- get_Lambda_RMSE(mdataset$Lambda_true, fit$Lambda)

  outside_percent <- get_95CI(mdataset$Lambda_true, fit$Lambda)

  metadata <- metadata(0, total_runtime, n_samples, lambda_RMSE, outside_percent)

  iterx <- ifelse(as.integer(n_samples)==0, 1L, as.integer(n_samples))
  m <- mfit(N=mdataset$N, D=mdataset$D, Q=mdataset$Q, iter=iterx, 
            Lambda=fit$Lambda, 
            Sigma=fit$Sigma,
            mdataset=mdataset,
            metadata=metadata)
  m$Eta <- fit$Eta
  m$Timer <- fit$Timer
  return(m)
}

