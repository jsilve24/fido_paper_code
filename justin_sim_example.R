require(rstan)
require(mongrel)
source("src/fit_methods.R")
source("src/dataset_methods.R")
source("src/utils.R")
source("src/fit_stan.R")
source("src/plotting.R")
source("src/fit_mongrel.R")
source("src/GH_standalone.R")

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

# If you make changes to mongrel, you can load them with devtools 
# without having to install - example below. 
# devtools::load_all("~/Research/src/mongrel")

# helper functions --------------------------------------------------------

#' function to calculate hessian for model at a given eta 
#' @param eta (D-1) x N matrix 
#' @details uses hessMongrelCollapsed function of mongrel 
hessMC <- function(mdataset, eta){
  X <- mdataset$X
  A <- solve(diag(mdataset$N)+ t(X)%*%mdataset$Gamma%*%X)
  hessMongrelCollapsed(mdataset$Y, mdataset$upsilon, 
                       mdataset$Theta%*%X, solve(mdataset$Xi), 
                       A, eta)
}

lap_approx <- function(n, mu, Sigma){
  p <- length(mu)
  L <- t(chol(Sigma))
  Z <- matrix(rnorm(n*p), p, n)
  X <- L%*%Z
  X <- sweep(X, 1, mu, FUN=`+`)
}


# simulation --------------------------------------------------------------

N <- 7L
Q <- 4L
D <- 5L
X <- rbind(1, matrix(rnorm((Q-1)*N), Q-1, N))
Lambda_true <- matrix(rnorm((D-1)*Q), D-1, Q)
Sigma_true <- solve(rWishart(1, D+10, diag(D-1))[,,1])

# Add a bit of sparsity by changing mean of intercept
#Lambda_true[,1] <- seq(0, (D-1)*.3, length=D-1) + Lambda_true[,1]



Theta <- matrix(0, D-1, Q)
upsilon <- D+10
Xi <- Sigma_true*(upsilon-(D-1)-1)

# account for increased uncertainty in intercept without acctually making
# Theta an informed prior
Gamma <- diag(c(10, rep(1, Q-1))) 

sim_data <- simulate_mdataset(N, size=rpois(N, 5000), X, Lambda_true, Sigma_true, 
                              Theta, Gamma, Xi, upsilon)


# Calculate Sparsity 
sum(sim_data$Y==0)/prod(dim(sim_data$Y))


# analysis ----------------------------------------------------------------

fit.sc <- fit_mstan(sim_data, parameterization="collapsed", ret_stanfit=FALSE)
# fit.sc2 <- fit_mstan(sim_data, parameterization="collapsed", ret_stanfit=FALSE, 
#                     iter=5000)
fit.sco <- fit_mstan_optim(sim_data, hessian=TRUE)
fit.sco2 <- fit_mstan_optim(sim_data, hessian=TRUE, ret_stanfit = TRUE)
#fit.mongrel1 <- fit_mongrel(sim_data, ret_mongrelfit = TRUE)
fit.mongrel2 <- fit_mongrel(sim_data, decomposition = "eigen")






# Check stan hessian against our hessian calculated at Stan's MAP estimate
fit.sco2$hessian[22:28,22:28]
hess.mongrel2 <- hessMC(sim_data, fit.sco2$par$eta)
hess.mongrel2[22:28,22:28]

fit.foo <- fit.mongrel2
fit.foo$Eta <- lap_approx(2000, c(fit.sco2$par$eta), solve(-hess.mongrel2))
dim(fit.foo$Eta) <- c(D-1, Q, 2000)


# # Finite Differences
# nll_partial <- function(x) nll(x, sim_data$Y, sim_data$X, sim_data$upsilon, 
#                                sim_data$Theta, sim_data$Xi, sim_data$Gamma)
# hess <- numDeriv::hessian(nll_partial, c(fit.sco2$par$eta))
# #nlme::fdHess(c(fit.sco2$par$eta), nll_partial)
# 
# apply(fit.sco$Lambda, c(1,2), mean)
# apply(fit.mongrel2$Lambda, c(1,2), mean)
# 
# plotting ----------------------------------------------------------------

plot_lambda(list("sc1"=fit.sc, 
                 "sco"=fit.sco, 
                 "me2"=fit.mongrel2, 
                 "foo" = fit.foo),
            Lambda_true=Lambda_true)

