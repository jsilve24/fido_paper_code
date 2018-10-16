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
        stop("Usage: Rscript simulate_efficiency.R {N} {D} {Q} {model=1,2,3} {log file} {opt: step_size} {opt: max_iter} {opt: b1} {opt: eps_f}")
        # Rscript simulate_efficiency.R 10 10 5 2 output.log 0.003 10000 0.9 1e-10
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

step_size <- 0.003
max_iter <- 10000
b1 <- 0.9
eps_f <- 1e-10
if (length(args) >= 9) {
	step_size <- as.numeric(args[6])
	max_iter <- as.integer(args[7])
	b1 <- as.numeric(args[8])
	eps_f <- as.numeric(args[9])
}

iter <- 2000L

rseed <- 2
set.seed(rseed)

# simulation --------------------------------------------------------------

sim_data <- simulate_with_hyperparams(N, D, Q, sparse=FALSE)

# Calculate Sparsity 
percent_zero <- sum(sim_data$Y==0)/prod(dim(sim_data$Y))
print(paste("Percent zero counts: ",percent_zero,sep=""))

# analysis ----------------------------------------------------------------

if(model == 1) {
  per_chain_it <- as.integer(iter/2)
  fit.sc <- fit_mstan(sim_data, parameterization="collapsed", ret_stanfit=FALSE, iter=per_chain_it)
  cat(paste("stan_collapsed,",fit.sc$metadata$mean_ess,",",fit.sc$metadata$warmup_runtime,",",
    fit.sc$metadata$total_runtime,",",N,",",D,",",Q,",",(4*per_chain_it),",",(2*per_chain_it),",",percent_zero,",",
    fit.sc$metadata$lambda_MSE,",",fit.sc$metadata$outside_95CI,",",rseed,"\n",sep=""),file=log_file, append=TRUE)
}

if(model == 2) {
  fit.mongrel.eigen <- fit_mongrel(sim_data, decomposition="eigen", step_size=step_size, max_iter=max_iter, b1=b1, eps_f=eps_f)
  cat(paste("mongrel_eigen,",fit.mongrel.eigen$metadata$mean_ess,",",fit.mongrel.eigen$metadata$warmup_runtime,",",
    fit.mongrel.eigen$metadata$total_runtime,",",N,",",D,",",Q,",",iter,",",0,",",percent_zero,",",
    fit.mongrel.eigen$metadata$lambda_MSE,",",fit.mongrel.eigen$metadata$outside_95CI,",",rseed,"\n",sep=""),file=log_file, append=TRUE)
}

if(model == 3) {
  fit.mongrel.cholesky <- fit_mongrel(sim_data, decomposition="cholesky", step_size=step_size, max_iter=max_iter, b1=b1, eps_f=eps_f)
  cat(paste("mongrel_cholesky,",fit.mongrel.cholesky$metadata$mean_ess,",",fit.mongrel.cholesky$metadata$warmup_runtime,",",
    fit.mongrel.cholesky$metadata$total_runtime,",",N,",",D,",",Q,",",iter,",",0,",",percent_zero,",",
    fit.mongrel.cholesky$metadata$lambda_MSE,",",fit.mongrel.cholesky$metadata$outside_95CI,",",rseed,"\n",sep=""),file=log_file, append=TRUE)
}
