require(mongrel)
require(driver)


#' Simulate mdataset from model (multinomial) 
#' Everything done in ALR_D
#' 
#' @param size N-vector of sequence depth of each sampel
#' @inheritParams mdataset
#' 
#' @details pulls D from nrow(Lambda_true)+1, pulls Q from ncol(Lambda_true), with 
#' no validation as to whether other passed objects have the same dimensions!
#' Uses Wishart distribution parameterization default to R and same as on wikipedia
#' - should be same as mongrel as well
#' Theta, Gamma, Xi, and Upsilon are just priors used added to the fitted dataset, 
#' they are not however used for the simulation - only Lambda_true and Sigma_true are
#' @examples 
#' N <- 10L
#' Q <- 3L
#' D <- 5L
#' size <- 1000:(1000+N)
#' X <- matrix(rnorm(Q*N), Q, N)
#' Lambda_true <- matrix(rnorm((D-1)*Q), D-1, Q)
#' Sigma_true <- Lambda_true %*% t(Lambda_true) #Lazy example
#' Theta <- Lambda_true # Lazy example
#' Gamma <- diag(Q)
#' Xi <- diag(D-1)
#' upsilon <- D
#' simulate_mdataset(N, size, X, Lambda_true, Sigma_true, Theta, Gamma, Xi, upsilon)
simulate_mdataset <- function(N, size, X, Lambda_true, Sigma_true, Theta, Gamma, Xi, upsilon){
  stopifnot(is.integer(N))
  D <- nrow(Lambda_true)+1L
  Q <- ncol(Lambda_true)
  E <- t(chol(Sigma_true)) %*% matrix(rnorm(N*(D-1)), D-1, N)
  eta <- Lambda_true%*%X + E
  pi <- driver::alrInv_array(eta, coords=1) # default is ALR_D
  Y <- matrix(0, D, N)
  for (j in 1:N){
    Y[,j] <- rmultinom(1, size=size[j], prob=pi[,j])
  }
  
  # Collect into mdataset object (with internal verify call)
  m <- mdataset(N=N, D=D, Q=Q, Y=Y, X=X, Lambda_true=Lambda_true, 
                Sigma_true=Sigma_true, Theta=Theta, Gamma=Gamma, Xi=Xi, 
                upsilon=upsilon)
  return(m)
}



# constructor for mdataset objects ----------------------------------------

#' Constructor for S3 mdataset (multinomial) objects
#' 
#' WARNING: Assumes everything is in ALR_D
#'
#' @param N number of samples (integer)
#' @param D number of categories (integer)
#' @param Q number of covariates (integer)
#' @param Y counts (D x N matrix)
#' @param X covariates (Q x N matrix)
#' @param Lambda_true true Lambda (D-1 x Q matrix)
#' @param Sigma_true true Sigma (D-1 x D-1 matrix)
#' @param Gamma prior (Q x Q matrix)
#' @param Xi prior (D-1 x D-1 matrix)
#' @param upsilon scalar > D-1
#'
#' @return object of class mdataset
#' @details internally verifys input
mdataset <- function(N, D, Q, Y, X, Lambda_true, Sigma_true, 
                     Theta, Gamma, Xi, upsilon){
  m <- new_mdataset(N, D, Q, Y, X, Lambda_true, Sigma_true, 
               Theta, Gamma, Xi, upsilon)
  verify(m)
  return(m)
}

# Keep internal (don't call directly from other code)
new_mdataset <- function(N, D, Q, Y, X, Lambda_true, Sigma_true, 
                        Theta, Gamma, Xi, upsilon){
  m <- list(N=N, D=D, Q=Q, Y=Y, X=X, Lambda_true=Lambda_true, Theta=Theta,  
            Sigma_true=Sigma_true, Gamma=Gamma, Xi=Xi, upsilon=upsilon)
  class(m) <- c("mdataset", "list")
  return(m)
}


# s3 methods --------------------------------------------------------------


#' Simple verification of passed mdataset (multinomial) object
#' 
#' WARNING: Assumes everything is in ALR_D
#' 
#' @param m an object of class mdataset
#' @param ... not used
#' @return throws error if any verification tests fail
verify.mdataset <- function(m, ...){
  stopifnot(is.integer(m$N), is.integer(m$D), is.integer(m$Q))
  N <- m$N; D <- m$D; Q <- m$Q;
  check_dims(m$Y, c(D, N), "Y")
  check_dims(m$X, c(Q, N), "X")
  check_dims(m$Lambda_true, c(D-1, Q), "Lambda_true")
  check_dims(m$Sigma_true, c(D-1, D-1), "Sigma_true")
  check_dims(m$Gamma, c(Q, Q), "Gamma")
  check_dims(m$Xi, c(D-1, D-1), "Xi")
  check_dims(m$upsilon, c(1), "upsilon")
  check_dims(m$Theta, c(D-1, Q), "Theta")
}
