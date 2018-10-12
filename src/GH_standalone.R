require(driver)

#' @param eta is  vec(eta) 
matt_nll <- function(eta, Y, X, upsilon, Theta, Xi, Gamma){
  D <- nrow(Y)
  N <- ncol(Y)
  eta <- matrix(eta, D-1, N)
  A <- solve(diag(N) + t(X)%*%Gamma %*% X)
  K <- solve(Xi)
  delta <- (upsilon+N+D-2)/2
  E <- eta - Theta%*%X
  S <- diag(D-1) + K%*%E%*%A%*%t(E)
  nll <- delta*log(det(S))
  return(nll)
}

#' @param eta is  vec(eta) 
mult_nll <- function(eta, Y, X, upsilon, Theta, Xi, Gamma){
  D <- nrow(Y)
  N <- ncol(Y)
  nll <- 0
  eta <- matrix(eta, D-1, N)
  pi <- t(driver::alrInv(t(eta)))
  for (i in 1:N){
    nll <- nll - dmultinom(Y[,i], prob=pi[,i], log=TRUE) # NEGATIVE!
  }
  return(nll)
}

#' @param eta is  vec(eta) 
nll <- function(eta, Y, X, upsilon, Theta, Xi, Gamma){
  nll <- 0
  nll <- nll+ matt_nll(eta, Y, X, upsilon, Theta, Xi, Gamma)
  nll <- nll+mult_nll(eta, Y, X, upsilon, Theta, Xi, Gamma)
  return(nll)
}