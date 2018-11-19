convert_to_seconds <- function(difftime_obj) {
	units <- attr(difftime_obj, "units")
	if (units == "secs") {
		return(as.double(difftime_obj))
	} else if (units == "mins") {
		return(60*as.double(difftime_obj))
	} else if (units == "hours") {
		return(60*60*as.double(difftime_obj))
	} else if (units == "days") {
		return(24*60*60*as.double(difftime_obj))
	} else if (units == "weeks") {
		return(7*24*60*60*as.double(difftime_obj))
	} else {
		return(-1)
	}
}

args = commandArgs(trailingOnly=TRUE)
if (length(args) < 3) {
        stop(paste("Usage: Rscript summarize_runs.R {which coord} {iterations} {output file}"))
        # Rscript summarize_runs.R N 2000 N_output.log
}

# need a try/catch here
coord <- args[1]
iter <- as.integer(args[2])
log_file <- args[3]

N_list <- NULL
D_list <- NULL
Q_list <- NULL
R_list <- c(1, 2, 3)

if(coord == "N") {
	N_list <- c(10, 20, 30, 50, 100, 250, 500, 750, 1000);
	D_list <- c(30);
	Q_list <- c(5);
} else if(coord == "D") {
	N_list <- c(100);
	D_list <- c(3, 5, 10, 25, 50, 75, 100, 250, 500);
	Q_list <- c(5);
} else if(coord == "Q") {
	N_list <- c(100);
	D_list <- c(30);
	Q_list <- c(2, 4, 10, 20, 50, 75, 100, 250, 500);
} else {
	quit()
}

model_list <- c("SU", "SC", "ME", "MC", "CLM")

per_chain_it <- as.integer(iter/2)
percent_zero <- -1

cat("model,ESS,burnin_runtime,sample_runtime,N,D,Q,total_iter,burnin_iter,percent_zero,lambda_MSE,percent_outside_95CI,random_seed\n", file=log_file)

for (m in model_list) {
	for (n in N_list) {
		for (d in D_list) {
			for (q in Q_list) {
				for (r in R_list) {
					cat("Checking",m,"N=",n,"D=",d,"Q=",q,"R=",r,"...\n")
					destfile <- paste("simulated_data/N",n,"_D",d,"_Q",q,"_R",r,".RData", sep="")
					if(file.exists(destfile)){
						# percent zero
						load(destfile)
						percent_zero <- sum(sim_data$Y==0)/prod(dim(sim_data$Y))
					} else {
						print(paste("Error: data file",destfile,"doesn't exist!"))
						percent_zero <- -1
					}
					destfile <- paste("fitted_models/",m,"_N",n,"_D",d,"_Q",q,"_R",r,".RData", sep="")
					if (file.exists(destfile)){
						load(destfile)
						if (m == "SU") {
							cat(paste("stan_uncollapsed,",fit.su$metadata$mean_ess,",",fit.su$metadata$warmup_runtime,",",
								fit.su$metadata$total_runtime,",",n,",",d,",",q,",",(4*per_chain_it),",",(2*per_chain_it),",",percent_zero,",",
								fit.su$metadata$lambda_MSE,",",fit.su$metadata$outside_95CI,",",r,"\n",sep=""),file=log_file,
								append=TRUE)
							rm(fit.su)
						} else if (m == "SC") {
							cat(paste("stan_collapsed,",fit.sc$metadata$mean_ess,",",fit.sc$metadata$warmup_runtime,",",
								fit.sc$metadata$total_runtime,",",n,",",d,",",q,",",(4*per_chain_it),",",(2*per_chain_it),",",percent_zero,",",
								fit.sc$metadata$lambda_MSE,",",fit.sc$metadata$outside_95CI,",",r,"\n",sep=""),file=log_file, append=TRUE)
							rm(fit.sc)
						} else if (m == "ME") {
							warmup_runtime <- fit.me$metadata$warmup_runtime
							total_runtime <- convert_to_seconds(fit.me$metadata$total_runtime)
							cat(paste("mongrel_eigen,",fit.me$metadata$mean_ess,",",warmup_runtime,",",total_runtime,",",n,",",d,",",q,",",iter*2,
								",",0,",",percent_zero,",",fit.me$metadata$lambda_MSE,",",fit.me$metadata$outside_95CI,",",r,"\n",sep=""),
								file=log_file, append=TRUE)
							rm(fit.me)
						} else if (m == "MC") {
							warmup_runtime <- fit.mc$metadata$warmup_runtime
							total_runtime <- convert_to_seconds(fit.mc$metadata$total_runtime)
							cat(paste("mongrel_cholesky,",fit.mc$metadata$mean_ess,",",warmup_runtime,",",total_runtime,",",n,",",d,",",q,",",iter*2,
								",",0,",",percent_zero,",",fit.mc$metadata$lambda_MSE,",",fit.mc$metadata$outside_95CI,",",r,"\n",sep=""),
								file=log_file, append=TRUE)
							rm(fit.mc)
						}
					} else {
						cat("Error: data file",destfile,"doesn't exist!\n")
					}
				}
			}
		}
	}
}
