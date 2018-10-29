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

N <- 100L
D <- 30L
Q <- 75L
rseed <- 2
model <- "me"

step_size <- 0.001
max_iter <- 50000
b1 <- 0.99
eps_f <- 1e-12
calcPartialHess <- FALSE

iter <- 2000L

# data already generated but for reproducible sampling behavior from Stan, optimizer
# worth setting random seed?
set.seed(rseed)

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

fit.me <- fit_mongrel(sim_data, decomposition="eigen", step_size=step_size, max_iter=max_iter, b1=b1, eps_f=eps_f)
