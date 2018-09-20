require(mongrel)
require(driver)


# constructor for mfit objects ----------------------------------------

#' Constructor for S3 mfit (multinomial) objects
#' 
#' WARNING: Assumes everything is in ALR_D
#'
#' @inheritParams mdataset
#' @param iter number of samples (integer)
#' @param Lambda D-1 x Q x iter array
#' @param Sigma D-1 x D-1 x iter array
#' @param mdataset mdataset object used for fit
#' 
#' 
#' @return object of class mfit
#' @details internally verifys input
mfit <- function(N, D, Q, iter, Lambda, Sigma, mdataset){
  m <- new_mfit(N, D, Q, iter, Lambda, Sigma, mdataset)
  verify(m)
  return(m)
}

# Keep internal (don't call directly from other code)
new_mfit <- function(N, D, Q, iter, Lambda, Sigma, mdataset){
  m <- list(N=N, D=D, Q=Q, iter=iter, Lambda=Lambda, Sigma=Sigma, 
            mdataset=mdataset)
  class(m) <- c("mfit", "list")
  return(m)
}


# s3 methods --------------------------------------------------------------



#' Simple verification of passed mfit (multinomial) object
#' 
#' WARNING: Assumes everything is in ALR_D
#' 
#' @param m an object of class mfit
#' @param ... not used
#' @return throws error if any verification tests fail
verify.mfit <- function(m, ...){
  stopifnot(is.integer(m$N), is.integer(m$D), is.integer(m$Q), is.integer(m$iter))
  N <- m$N; D <- m$D; Q <- m$Q; iter <- m$iter
  check_dims(m$Lambda, c(D-1, Q, iter), "Lambda")
  check_dims(m$Sigma, c(D-1, D-1, iter), "Sigma")
  verify.mdataset(m$mdataset)
}
