require(mongrel)
source("src/fit_methods.R")
source("src/dataset_methods.R")
source("src/utils.R")
source("src/fit_mongrel.R")

# If you make changes to mongrel, you can load them with devtools 
# without having to install - example below. 
devtools::load_all("/data/mukherjeelab/Mongrel/mongrel")
#devtools::load_all("C:/Users/kim/Documents/mongrel")

# need a try/catch here
N <- 100L
D <- 30L
Q <- 75L
model <- "me"

step_size <- 0.004
max_iter <- 50000
b1 <- 0.99
eps_f <- 1e-10
calcPartialHess <- FALSE

iter <- 2000L

warnings <- 0
errors <- 0

for(rseed in 1:100) {
  set.seed(rseed)
  sim_data <- simulate_with_hyperparams(N, D, Q, sparse=FALSE)
  # Calculate Sparsity 
  percent_zero <- sum(sim_data$Y==0)/prod(dim(sim_data$Y))
  cat(paste("Percent zero counts: ",percent_zero,"\n",sep=""))
  result = tryCatch({
    fit.me <- fit_mongrel(sim_data, decomposition="eigen", step_size=step_size, max_iter=max_iter, b1=b1, eps_f=eps_f)
  }, warning = function(w) {
    cat("warningstringXYZ")
  }, error = function(e) {
    cat("errorstringXYZ")
  }, finally = {})
}

cat("Warnings:",warnings,"\n")
cat("Errors:",errors,"\n")
