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
#' @param ... other parameters passed to function mongrel
#' @return mfit object
fit_mongrel <- function(mdataset, n_samples=2000, decomposition="eigen", ...){
  
  fit <- mongrel(Y=mdataset$Y, X=mdataset$X, upsilon=mdataset$upsilon, 
                 Theta=mdataset$Theta, Gamma=mdataset$Gamma, Xi=mdataset$Xi, 
                 init=random_mongrel_init(mdataset$Y), 
                 decomp_method=decomposition, n_samples=n_samples,  ...)

  m <- mfit(N=mdataset$N, D=mdataset$D, Q=mdataset$Q, iter=as.integer(n_samples), 
            Lambda=fit$Lambda, 
            Sigma=fit$Sigma, 
            mdataset=mdataset)
  
  return(m)
}


