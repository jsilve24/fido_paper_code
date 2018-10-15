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

args = commandArgs(trailingOnly=TRUE)
if (length(args) < 5) {
        stop("Usage: Rscript simulate_efficiency.R {N} {D} {Q} {model=1,2,3} {log file}")
        # Rscript simulate_efficiency.R 10 10 5 2 output.log
	# model 1 : Stan (collapsed)
	#       2 : Mongrel (eigendecomposition)
	#       3 : Mongrel (Cholesky)
}

# need a try/catch here
N <- as.integer(args[1])
D <- as.integer(args[2])
Q <- as.integer(args[3])
model <- as.integer(args[4])
log_file <- args[5]

iter <- 2000L

rseed <- 2
set.seed(rseed)

# simulation --------------------------------------------------------------

sim_data <- simulate_with_hyperparams(N, D, Q, sparse=FALSE)

# Calculate Sparsity 
print(paste("Percent zeros: ", sum(sim_data$Y==0)/prod(dim(sim_data$Y)), sep=""))

# analysis ----------------------------------------------------------------

if(model == 1) {
  per_chain_it <- as.integer(iter/2)
  fit.sc <- fit_mstan(sim_data, parameterization="collapsed", ret_stanfit=FALSE, iter=per_chain_it)
  cat(paste("stan_collapsed,",fit.sc$metadata$mean_ess,",",fit.sc$metadata$warmup_runtime,",",
    fit.sc$metadata$total_runtime,",",N,",",D,",",Q,",",(4*per_chain_it),",",(2*per_chain_it),",",rseed,"\n",sep=""),file=log_file, append=TRUE)
}

if(model == 2) {
  fit.mongrel.eigen <- fit_mongrel(sim_data, decomposition = "eigen")
  cat(paste("mongrel_eigen,",fit.mongrel.eigen$metadata$mean_ess,",",fit.mongrel.eigen$metadata$warmup_runtime,",",
    fit.mongrel.eigen$metadata$total_runtime,",",N,",",D,",",Q,",",iter,",",0,",",rseed,"\n",sep=""),file=log_file, append=TRUE)
}

if(model == 3) {
  fit.mongrel.cholesky <- fit_mongrel(sim_data, decomposition = "cholesky")
  cat(paste("mongrel_cholesky,",fit.mongrel.cholesky$metadata$mean_ess,",",fit.mongrel.cholesky$metadata$warmup_runtime,",",
    fit.mongrel.cholesky$metadata$total_runtime,",",N,",",D,",",Q,",",iter,",",0,",",rseed,"\n",sep=""),file=log_file, append=TRUE)
}
