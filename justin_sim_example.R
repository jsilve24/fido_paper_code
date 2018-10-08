require(rstan)
require(mongrel)
source("src/fit_methods.R")
source("src/dataset_methods.R")
source("src/utils.R")
source("src/fit_stan.R")
source("src/plotting.R")
source("src/fit_mongrel.R")

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)



# simulation --------------------------------------------------------------

N <- 50L
Q <- 4L
D <- 15L
X <- rbind(1, matrix(rnorm((Q-1)*N), Q-1, N))
Lambda_true <- matrix(rnorm((D-1)*Q), D-1, Q)
Sigma_true <- solve(rWishart(1, D+10, diag(D-1))[,,1])

# Add a bit of sparsity by changing mean of intercept
Lambda_true[,1] <- seq(0, (D-1)*.3, length=D-1) + Lambda_true[,1]



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
fit.mongrel1 <- fit_mongrel(sim_data)
fit.mongrel2 <- fit_mongrel(sim_data)
fit.mongrel3 <- fit_mongrel(sim_data)
fit.mongrel4 <- fit_mongrel(sim_data)

# plotting ----------------------------------------------------------------

plot_lambda(list("sc"=fit.sc, 
                 "me1"=fit.mongrel1, 
                 "me2"=fit.mongrel2, 
                 "me3"=fit.mongrel3, 
                 "me4"=fit.mongrel4), 
            Lambda_true = Lambda_true)


