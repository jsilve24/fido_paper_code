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
devtools::load_all("C:/Users/kim/Documents/mongrel")

# simulation --------------------------------------------------------------

sim_data <- simulate_with_hyperparams(N=10L, D=10L, Q=5L, sparse=TRUE)

# Calculate Sparsity 
sum(sim_data$Y==0)/prod(dim(sim_data$Y))

# analysis ----------------------------------------------------------------

fit.sc <- fit_mstan(sim_data, parameterization="collapsed", ret_stanfit=FALSE)
fit.su <- fit_mstan(sim_data, parameterization="uncollapsed", ret_stanfit=FALSE, iter = 3000)
#fit.sco <- fit_mstan_optim(sim_data, hessian=TRUE)
fit.mongrel <- fit_mongrel(sim_data, decomposition = "eigen")

# plotting ----------------------------------------------------------------

plot_lambda(list("sc" = fit.sc, 
                 "su" = fit.su, 
                 "me" = fit.mongrel),
            Lambda_true=sim_data$Lambda_true)

plot_eta(list("sc" = fit.sc, 
              "su" = fit.su, 
              "me" = fit.mongrel))
