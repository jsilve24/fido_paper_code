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
if (length(args) < 6) {
        stop(paste("Usage: Rscript simulate_efficiency.R {N} {D} {Q} {random seed} {model='sc','su','me','mc','mcp','clm'}",
                   "{log file} {opt: step_size} {opt: max_iter} {opt: b1} {opt: eps_f}"))
        # Rscript simulate_efficiency.R 10 10 5 1 clm output.log 0.002 50000 0.99 1e-10
	# model 'clm' : conjugate linear model
	#       'sc'  : Stan (collapsed)
	#       'su'  : Stan (uncollapsed)
	#       'me'  : Mongrel (eigendecomposition)
	#       'mc'  : Mongrel (Cholesky)
	#       'mcp' : Mongrel (Cholesky), partial (multinomial) Hessian
}

# need a try/catch here
N <- as.integer(args[1])
D <- as.integer(args[2])
Q <- as.integer(args[3])
rseed <- as.integer(args[4])
model <- args[5]
log_file <- args[6]

step_size <- 0.003
max_iter <- 10000
b1 <- 0.9
eps_f <- 1e-10
calcPartialHess <- FALSE
if (length(args) >= 9) {
	step_size <- as.numeric(args[7])
	max_iter <- as.integer(args[8])
	b1 <- as.numeric(args[9])
	eps_f <- as.numeric(args[10])
}

iter <- 2000L

# data already generated but for reproducible sampling behavior from Stan, optimizer
# worth setting random seed?
set.seed(rseed)

writeLines(c("\nUsing parameter values",
             paste("\tN:",N),
             paste("\tD:",D),
             paste("\tQ:",Q),
             paste("\trseed:",rseed),
             paste("\tmodel:",model),
             paste("\tlog_file:",log_file),
             paste("\tstep_size:",step_size),
             paste("\tmax_iter:",max_iter),
             paste("\tb1:",b1),
             paste("\teps_f:",eps_f),
             paste("\tcalcPartialHess:",calcPartialHess),
             paste("\titer:",iter)
          ))

model_save_dir = "fitted_models"

# simulation --------------------------------------------------------------

destfile <- paste("simulated_data/N", N, "_D", D, "_Q", Q, "_R", rseed, ".RData", sep="")
if(!file.exists(destfile)){
  print(paste("Error: data file",destfile,"doesn't exist!"))
  quit()
}
load(destfile)

# Calculate Sparsity 
percent_zero <- sum(sim_data$Y==0)/prod(dim(sim_data$Y))
print(paste("Percent zero counts: ",percent_zero,sep=""))

# analysis ----------------------------------------------------------------

if(model == 'clm') {
  fit.clm <- conjugateLinearModel(sim_data$Y, sim_data$X, sim_data$Theta, sim_data$Gamma, sim_data$Xi, sim_data$upsilon, n_samples=iter)
  save(fit.clm, file=paste(model_save_dir,"/CLM_N",N,"_D",D,"_Q",Q,"_R",rseed,".RData",sep=""))
}

if(model == 'sc') {
  per_chain_it <- as.integer(iter/2)
  fit.sc <- fit_mstan(sim_data, parameterization="collapsed", ret_stanfit=FALSE, iter=per_chain_it)
  cat(paste("stan_collapsed,",fit.sc$metadata$mean_ess,",",fit.sc$metadata$warmup_runtime,",",
    fit.sc$metadata$total_runtime,",",N,",",D,",",Q,",",(4*per_chain_it),",",(2*per_chain_it),",",percent_zero,",",
    fit.sc$metadata$lambda_MSE,",",fit.sc$metadata$outside_95CI,",",rseed,"\n",sep=""),file=log_file, append=TRUE)
  save(fit.sc, file=paste(model_save_dir,"/SC_N",N,"_D",D,"_Q",Q,"_R",rseed,".RData",sep=""))
}

if(model == 'su') {
  per_chain_it <- as.integer(iter/2)
  fit.su <- fit_mstan(sim_data, parameterization="uncollapsed", ret_stanfit=FALSE, iter=per_chain_it)
  cat(paste("stan_uncollapsed,",fit.su$metadata$mean_ess,",",fit.su$metadata$warmup_runtime,",",
    fit.su$metadata$total_runtime,",",N,",",D,",",Q,",",(4*per_chain_it),",",(2*per_chain_it),",",percent_zero,",",
    fit.su$metadata$lambda_MSE,",",fit.su$metadata$outside_95CI,",",rseed,"\n",sep=""),file=log_file, append=TRUE)
  save(fit.su, file=paste(model_save_dir,"/SU_N",N,"_D",D,"_Q",Q,"_R",rseed,".RData",sep=""))
}

if(model == 'me') {
  fit.me <- fit_mongrel(sim_data, decomposition="eigen", step_size=step_size, max_iter=max_iter, b1=b1, eps_f=eps_f)
  cat(paste("mongrel_eigen,",fit.me$metadata$mean_ess,",",fit.me$metadata$warmup_runtime,",",
    fit.me$metadata$total_runtime,",",N,",",D,",",Q,",",iter,",",0,",",percent_zero,",",
    fit.me$metadata$lambda_MSE,",",fit.me$metadata$outside_95CI,",",rseed,"\n",sep=""),file=log_file, append=TRUE)
  save(fit.me, file=paste(model_save_dir,"/ME_N",N,"_D",D,"_Q",Q,"_R",rseed,".RData",sep=""))
}

if(model == 'mc') {
  fit.mc <- fit_mongrel(sim_data, decomposition="cholesky", step_size=step_size, max_iter=max_iter, b1=b1, eps_f=eps_f)
  cat(paste("mongrel_cholesky,",fit.mc$metadata$mean_ess,",",fit.mc$metadata$warmup_runtime,",",
    fit.mc$metadata$total_runtime,",",N,",",D,",",Q,",",iter,",",0,",",percent_zero,",",
    fit.mc$metadata$lambda_MSE,",",fit.mc$metadata$outside_95CI,",",rseed,"\n",sep=""),file=log_file, append=TRUE)
  save(fit.mc, file=paste(model_save_dir,"/MC_N",N,"_D",D,"_Q",Q,"_R",rseed,".RData",sep=""))
}

if(model == 'mcp') {
  fit.mcp <- fit_mongrel(sim_data, decomposition="cholesky", step_size=step_size, max_iter=max_iter, b1=b1, eps_f=eps_f, calcPartialHess=TRUE)
  cat(paste("mongrel_cholesky_partial,",fit.mcp$metadata$mean_ess,",",fit.mcp$metadata$warmup_runtime,",",
    fit.mcp$metadata$total_runtime,",",N,",",D,",",Q,",",iter,",",0,",",percent_zero,",",
    fit.mcp$metadata$lambda_MSE,",",fit.mcp$metadata$outside_95CI,",",rseed,"\n",sep=""),file=log_file, append=TRUE)
  save(fit.mcp, file=paste(model_save_dir,"/MCP_N",N,"_D",D,"_Q",Q,"_R",rseed,".RData",sep=""))
}
