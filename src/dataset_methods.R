require(mongrel)
require(driver)

#' Simulate hyperparameters, then data from the model (multinomial)
#' Returns an mdataset
#' 
#' @param N number of samples (integer)
#' @param D number of categories (integer)
#' @param Q number of covariates (integer)
#' @param sparse logical, if TRUE induces more zero counts in Y
#' 
#' @details Generates random X, Lambda_true, Sigma_true, Theta, Xi, Gamma, upsilon
#' @examples 
#' simulate_with_hyperparams(15L, 12L, 8L, TRUE)
simulate_with_hyperparams <- function(N=10L, D=10L, Q=5L, sparse=FALSE){
  stopifnot(is.integer(N), is.integer(D), is.integer(Q))
  
  X <- rbind(1, matrix(rnorm((Q-1)*N), Q-1, N))
  Lambda_true <- matrix(rnorm((D-1)*Q), D-1, Q)
  Sigma_true <- solve(rWishart(1, D+10, diag(D-1))[,,1])
  
  if(sparse) {
    # Add a bit of sparsity by changing mean of intercept
    Lambda_true[,1] <- seq(0, (D-1)*.3, length=D-1) + Lambda_true[,1]
    # account for increased uncertainty in intercept without acctually making
    # Theta an informed prior
    Gamma <- diag(c(10, rep(1, Q-1))) 
  }
  
  Theta <- matrix(0, D-1, Q)
  upsilon <- D+10
  Xi <- Sigma_true*(upsilon-(D-1)-1)
  Gamma <- diag(Q)
 
  sim_data <- simulate_mdataset(N, size=rpois(N, 5000), X, Lambda_true, Sigma_true, 
                                Theta, Gamma, Xi, upsilon)
  return(sim_data)
}

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
#' Q <- 4L
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
                Sigma_true=Sigma_true, Theta=Theta, Gamma=Gamma, Xi=Xi, upsilon=upsilon)
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



# constructor for metadata objects ----------------------------------------

#' Constructor for S3 metadata
#' 
#' @param warmup_runtime seconds of runtime used in warm-up or burn-in
#' @param total_runtime seconds of total execution time
#' @param mean_ess mean effective sample size over the Lambdas
#' @param lambda_RMSE root mean squared error of estimated and true Lambdas
#' @param outside_95CI proportion of true Lambda values outside 95% posterior CI
#'
#' @return object of class metadata
#' @details 
metadata <- function(warmup_runtime, total_runtime, mean_ess, lambda_RMSE, outside_95CI){
  m <- list(warmup_runtime=warmup_runtime, total_runtime=total_runtime,
            mean_ess=mean_ess, lambda_RMSE=lambda_RMSE, outside_95CI=outside_95CI)
  class(m) <- c("metadata", "list")
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


# internal
get_Lambda_RMSE <- function(Lambda_true, Lambda_estimated, RMSE=TRUE) {
  ref <- c(Lambda_true)
  estimates <- Lambda_estimated
  dim(estimates) <- c(dim(Lambda_estimated)[1]*dim(Lambda_estimated)[2],dim(Lambda_estimated)[3])
  MSE <- mean(apply(estimates, 2, function(x, y) { (x-y)**2 }, y=ref))
  if(RMSE) {
    return(sqrt(MSE))
  } else {
    return(MSE)
  }
}

lambda_outside_bounds <- function(x) {
  # x[1] is the 0.025 quantile
  # x[2] is the 0.975 quantile
  # x[3] is the value we're testing for inclusion
  if(x[3] < x[1] || x[3] > x[2]) {
    return(1)
  }
  return(0)
}

get_95CI <- function(Lambda_true, Lambda_estimated) {
  ref <- c(Lambda_true)
  intervals <- apply(Lambda_estimated, c(1,2), function(x) quantile(x, probs=c(0.025, 0.975)))
  lower <- c(intervals[1,,])
  upper <- c(intervals[2,,])
  return(sum(apply(rbind(lower, upper, ref), 2, lambda_outside_bounds))/(dim(Lambda_true)[1]*dim(Lambda_true)[2]))
}




