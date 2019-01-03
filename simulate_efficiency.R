require(rstan)
require(mongrel)
source("src/fit_methods.R")
source("src/dataset_methods.R")
source("src/utils.R")
source("src/fit_mongrel.R")
source("src/fit_stan.R")
source("src/fit_CLM.R")
source("src/plotting.R")
source("src/GH_standalone.R")

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

# If you make changes to mongrel, you can load them with devtools 
# without having to install - example below. 
devtools::load_all("/data/mukherjeelab/Mongrel/mongrel")
#devtools::load_all("C:/Users/kim/Documents/mongrel")

args = commandArgs(trailingOnly=TRUE)
if (length(args) < 5) {
        stop(paste("Usage: Rscript simulate_efficiency.R {N} {D} {Q} {random seed} {model='sc','su','me','mc','mcp','clm','svbcm','svbcf','svbum','svbuf'}",
                   "{opt: iterations} {opt: optim_method} {opt: save_dir} {opt: file suffix} {opt: max_iter} {opt: eps_f} {opt: step_size} {opt: b1}"))
        # e.g. Rscript simulate_efficiency.R 20 30 5 1 mc 2000 adam fitted_models_temp _MKL_4 100000 1e-11 0.004 0.99
        # e.g. Rscript simulate_efficiency.R 20 30 5 1 mc 0 lbfgs fitted_models_temp _MKL_4 100000 1e-11
	# model 'clm'   : conjugate linear model
	#       'sc'    : Stan (collapsed)
	#       'su'    : Stan (uncollapsed)
	#       'me'    : Mongrel (eigendecomposition)
	#       'mc'    : Mongrel (Cholesky)
	#       'mcp'   : Mongrel (Cholesky), partial (multinomial) Hessian
	#	'svbcm'   : Stan (collapsed) - variational Bayes, meanfield
	#	'svbcf'   : Stan (collapsed) - variational Bayes, fullrank
	#	'svbum'   : Stan (uncollapsed) - variational Bayes, meanfield
	#	'svbuf'   : Stan (uncollapsed) - variational Bayes, fullrank
}

# required arguments
N <- as.integer(args[1])
D <- as.integer(args[2])
Q <- as.integer(args[3])
rseed <- as.integer(args[4])
model <- args[5]

# optional arguments
# 6: iterations
# 7: optimization method
# 8: save directory
# 9: save file suffix
# 10: optimization parameter, max iterations
# 11: optimization parameter, gradient stopping threshold
# 12: optimization parameter, step size (Adam only)
# 13: optimization parameter, momentum (Adam only

iter <- 2000L
if(length(args) >= 6) {
  iter <- as.integer(args[6])
}

optim_method <- "adam"
if(length(args) >= 7) {
  optim_method <- args[7]
}

model_save_dir <- "fitted_models"
if(length(args) >= 8) {
  model_save_dir <- args[8]
}

file_suffix <- NULL
if (length(args) >= 9) {
  file_suffix <- args[9]
}
if(iter < 1) {
  file_suffix <- paste(file_suffix, "MAP", sep="_")
}

max_iter <- 10000
eps_f <- 1e-10
if(length(args) >= 11) {
  max_iter <- as.integer(args[10])
  eps_f <- as.numeric(args[11])
}

step_size <- 0.003
b1 <- 0.9
if(length(args) >= 13) {
  step_size <- as.numeric(args[12])
  b1 <- as.numeric(args[13])
}

# data already generated but for reproducible sampling behavior from Stan, optimizer
# worth setting random seed?
set.seed(rseed)

writeLines(c("\nUsing parameter values",
             paste("\tN:",N),
             paste("\tD:",D),
             paste("\tQ:",Q),
             paste("\trseed:",rseed),
             paste("\titer:",iter),
             paste("\tmodel:",model),
             paste("\toptim_method:",optim_method),
             paste("\tsaving to:",model_save_dir),
             paste("\tmax_iter:",max_iter),
             paste("\teps_f:",eps_f),
             paste("\tstep_size:",step_size),
             paste("\tb1:",b1)
          ))

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
  fit.clm <- fit_CLM(sim_data, n_samples=iter)
  save(fit.clm, file=paste(model_save_dir,"/CLM_N",N,"_D",D,"_Q",Q,"_R",rseed,file_suffix,".RData",sep=""))

  if(FALSE) {
    # "sample" with a bunch of uncollapses; this should be the same as above but we'll use it to double-check
    Lambdas <- array(dim=c(D-1,Q,iter))
    set.seed(1)
    fit.u <- uncollapseMongrelCollapsed(eta.hat, sim_data$X, sim_data$Theta, sim_data$Gamma, sim_data$Xi, sim_data$upsilon, ret_mean=FALSE)
    Lambdas[,,1] <- fit.u$Lambda
    for (i in 2:iter) {
      set.seed(i)
      temp <- uncollapseMongrelCollapsed(eta.hat, sim_data$X, sim_data$Theta, sim_data$Gamma, sim_data$Xi, sim_data$upsilon, ret_mean=FALSE)
      Lambdas[,,i] <- temp$Lambda
    }
    fit.u$Lambda <- Lambdas
    save(fit.u, file=paste(model_save_dir,"/CLMU_N",N,"_D",D,"_Q",Q,"_R",rseed,file_suffix,".RData",sep=""))
  }
}

if(model == 'sc') {
  per_chain_it <- as.integer(iter/2)
  fit.sc <- fit_mstan(sim_data, ret_stanfit=FALSE, iter=per_chain_it)
  save(fit.sc, file=paste(model_save_dir,"/SC_N",N,"_D",D,"_Q",Q,"_R",rseed,file_suffix,".RData",sep=""))
}

if(model == 'su') {
  per_chain_it <- as.integer(iter/2)
  fit.su <- fit_mstan(sim_data, parameterization="uncollapsed", ret_stanfit=FALSE, iter=per_chain_it)
  save(fit.su, file=paste(model_save_dir,"/SU_N",N,"_D",D,"_Q",Q,"_R",rseed,file_suffix,".RData",sep=""))
}

if(model == 'svbcm') {
  fit.svbcm <- fit_mstan_vb(sim_data, ret_stanfit=FALSE, iter=iter)
  save(fit.svbcm, file=paste(model_save_dir,"/SVBCM_N",N,"_D",D,"_Q",Q,"_R",rseed,file_suffix,".RData",sep=""))
}

if(model == 'svbcf') {
  fit.svbcf <- fit_mstan_vb(sim_data, algorithm="fullrank", ret_stanfit=FALSE, iter=iter)
  save(fit.svbcf, file=paste(model_save_dir,"/SVBCF_N",N,"_D",D,"_Q",Q,"_R",rseed,file_suffix,".RData",sep=""))
}

if(model == 'svbum') {
  fit.svbum <- fit_mstan_vb(sim_data, parameterization="uncollapsed", ret_stanfit=FALSE, iter=iter)
  save(fit.svbum, file=paste(model_save_dir,"/SVBUM_N",N,"_D",D,"_Q",Q,"_R",rseed,file_suffix,".RData",sep=""))
}

if(model == 'svbuf') {
  fit.svbuf <- fit_mstan_vb(sim_data, parameterization="uncollapsed", algorithm="meanfield", ret_stanfit=FALSE, iter=iter)
  save(fit.svbuf, file=paste(model_save_dir,"/SVBUF_N",N,"_D",D,"_Q",Q,"_R",rseed,file_suffix,".RData",sep=""))
}

calcGradHess <- TRUE
if(iter == 0) {
  calcGradHess <- FALSE
}

if(model == 'me') {
  # double Mongrel iterations
  iter <- 2L*iter
  fit.me <- fit_mongrel(sim_data, decomposition="eigen", optim_method=optim_method, n_samples=iter, calcGradHess=calcGradHess, step_size=step_size, max_iter=max_iter, b1=b1, eps_f=eps_f)
  save(fit.me, file=paste(model_save_dir,"/ME_N",N,"_D",D,"_Q",Q,"_R",rseed,file_suffix,".RData",sep=""))
}

if(model == 'mc') {
  # double Mongrel iterations
  iter <- 2L*iter
  fit.mc <- fit_mongrel(sim_data, decomposition="cholesky", optim_method=optim_method, n_samples=iter, step_size=step_size, calcGradHess=calcGradHess, max_iter=max_iter, b1=b1, eps_f=eps_f)
  save(fit.mc, file=paste(model_save_dir,"/MC_N",N,"_D",D,"_Q",Q,"_R",rseed,file_suffix,".RData",sep="")) 
}

if(model == 'mcp') {
  fit.mcp <- fit_mongrel(sim_data, decomposition="cholesky", optim_method=optim_method, n_samples=iter, step_size=step_size, calcGradHess=calcGradHess, max_iter=max_iter, b1=b1, eps_f=eps_f, calcPartialHess=TRUE)
  save(fit.mcp, file=paste(model_save_dir,"/MCP_N",N,"_D",D,"_Q",Q,"_R",rseed,file_suffix,".RData",sep=""))
}
