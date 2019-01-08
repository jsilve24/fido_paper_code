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
devtools::load_all("/data/mukherjeelab/Mongrel/mongrel")
#devtools::load_all("C:/Users/kim/Documents/mongrel")

for(n in 1:1) {

  # simulation --------------------------------------------------------------

  sim_data <- simulate_with_hyperparams(N=as.integer(n*10), D=10L, Q=5L, sparse=FALSE)

  # Calculate Sparsity 
  print(paste("Percent zeros: ", sum(sim_data$Y==0)/prod(dim(sim_data$Y)), sep=""))

  # analysis ----------------------------------------------------------------

  fit.sc <- fit_mstan(sim_data, parameterization="collapsed", ret_stanfit=FALSE)
  fit.su <- fit_mstan(sim_data, parameterization="uncollapsed", ret_stanfit=FALSE)
  #fit.sco <- fit_mstan_optim(sim_data, hessian=TRUE)
  fit.mongrel.eigen <- fit_mongrel(sim_data, decomposition = "eigen")
  fit.mongrel.cholesky <- fit_mongrel(sim_data, decomposition = "cholesky")
  fit.mongrel.cholesky.partial <- fit_mongrel(sim_data, decomposition = "cholesky", calcPartialHess = TRUE)

  # plotting ----------------------------------------------------------------

  plot_lambda(list("sc" = fit.sc,
                   "su" = fit.su,
                   "me" = fit.mongrel.eigen, 
                   "mc" = fit.mongrel.cholesky,
                   "mc_partial" = fit.mongrel.cholesky.partial),
              Lambda_true=sim_data$Lambda_true,
              image_filename=paste("lambda_N",as.character(n*10),"_D10_Q5.png", sep=""))

#  plot_eta(list("sc" = fit.sc, 
#                "me" = fit.mongrel.eigen, 
#                "mc" = fit.mongrel.cholesky,
#                "mc_partial" = fit.mongrel.cholesky.partial),
#           image_filename=paste("eta_N",as.character(n*10),"_D10_Q5.png", sep=""))

}
