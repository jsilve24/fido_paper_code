require(rstan)
require(mongrel)
source("src/fit_methods.R")
source("src/dataset_methods.R")
source("src/utils.R")
source("src/fit_stan.R")
source("src/plotting.R")



# stupid simulation -------------------------------------------------------

N <- 10L
Q <- 4L
D <- 5L
size <- 1000:(1000+N)
X <- matrix(rnorm(Q*N), Q, N)
Lambda_true <- matrix(rnorm((D-1)*Q), D-1, Q)
Sigma_true <- Lambda_true %*% t(Lambda_true) #Lazy example
Theta <- Lambda_true # Lazy example
Gamma <- diag(Q)
Xi <- diag(D-1)
upsilon <- D
sim <- simulate_mdataset(N, size, X, Lambda_true, Sigma_true, Theta, Gamma, Xi, upsilon)

fit1 <- fit_mstan(sim, parameterization="collapsed",  ret_stanfit=TRUE)
fit2 <- fit_mstan(sim, parameterization = "uncollapsed", ret_stanfit=TRUE)



# Calculate mean ESS
summary(fit1)$summary[,"n_eff"] -> foo1
summary(fit2)$summary[,"n_eff"] -> foo2
mean(foo1[-length(foo1)]) # Drop l_prob which for some reason shows up at end 
mean(foo2[-length(foo2)]) # of this vector ... weird


fit1 <- fit_mstan(sim, parameterization="collapsed",  ret_stanfit=FALSE)
fit2 <- fit_mstan(sim, parameterization = "uncollapsed", ret_stanfit=FALSE)
mfits <- list("stan_collapsed"=fit1, "stan_uncollapsed"=fit2)
plot_lambda(mfits)

