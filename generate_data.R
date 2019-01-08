source("src/fit_methods.R")
source("src/dataset_methods.R")
source("src/utils.R")

args = commandArgs(trailingOnly=TRUE)
if (length(args) < 5) {
  stop(paste("Usage: Rscript generate_data.R {N} {D} {Q} {output dir}"))
  # Rscript generate_data.R 10 10 5 1 simulated_data
  # Rscript generate_data.R 5000 1000 20 1 simulated_data
}

# need a try/catch here
N <- as.integer(args[1])
D <- as.integer(args[2])
Q <- as.integer(args[3])
rseed <- as.integer(args[4])
out_dir <- args[5]

set.seed(rseed)

# simulation --------------------------------------------------------------

# don't overwrite if already created
destfile <- paste(out_dir, "/N", N, "_D", D, "_Q", Q, "_R", rseed, ".RData", sep="")
if(!file.exists(destfile)){
	sim_data <- simulate_with_hyperparams(N, D, Q, sparse=FALSE)
	percent_zero <- sum(sim_data$Y==0)/prod(dim(sim_data$Y))
	print(paste("Percent zero counts: ",percent_zero,sep=""))
	save(sim_data, file=destfile)
}
