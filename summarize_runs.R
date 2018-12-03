source("src/dataset_methods.R")

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
if (length(args) < 1) {
        stop(paste("Usage: Rscript summarize_runs.R {output file}"))
        # Rscript summarize_runs.R all_output.log
}

# need a try/catch here
log_file <- args[1]

# WARNING: NEED TO FIX IN fit_Stan(), CURRENTLY STORING sample_runtime IN MFIT'S total_runtime
# correctly referencing it as sample_runtime here but need to fix this for clarity/consistency
cat(paste("model,ESS,burnin_runtime,sample_runtime,Timer_HessianCalculation,Timer_LaplaceApproximation,Timer_Optimization,Timer_Overall,Timer_Uncollapse_Overall.Overall,",
    "sweep_param,sweep_value,total_iter,burnin_iter,percent_zero,lambda_RMSE,percent_outside_95CI,random_seed\n",sep=""), file=log_file)

write_Stan_out <- function(fit, model_name) {
  lambda_RMSE <- get_Lambda_MSE(Lambda_true, fit$Lambda)
  outside_CI <- get_95CI(Lambda_true, fit$Lambda)
  cat(paste(model_name,",",
    fit$metadata$mean_ess,",",
    fit$metadata$warmup_runtime,",",
    fit$metadata$total_runtime,",0,0,0,0,0,",
    coord,",",
    sweep_param_value,",",
    dim(fit$Lambda)[3],",",
    dim(fit$Lambda)[3]/2,",",
    percent_zero,",",
    lambda_RMSE,",",
    outside_CI,",",
    r,
    "\n",sep=""),file=log_file, append=TRUE)
}

write_Mongrel_out <- function(fit, model_name) {
  warmup_runtime <- 0
  total_runtime <- convert_to_seconds(fit$metadata$total_runtime)
  lambda_RMSE <- get_Lambda_MSE(Lambda_true, fit$Lambda)
  outside_CI <- get_95CI(Lambda_true, fit$Lambda)
  cat(paste(model_name,",",
    dim(fit$Lambda)[3],",",
    warmup_runtime,",",
    total_runtime,",",
    fit$Timer['HessianCalculation'][[1]],",",
    fit$Timer['LaplaceApproximation'][[1]],",",
    fit$Timer['Optimization'][[1]],",",
    fit$Timer['Overall'][[1]],",",
    fit$Timer['Uncollapse_Overall.Overall'][[1]],",",
    coord,",",
    sweep_param_value,",",
    dim(fit$Lambda)[3],",",
    0,",",
    percent_zero,",",
    lambda_RMSE,",",
    outside_CI,",",
    r,
    "\n",sep=""),
    file=log_file, append=TRUE)
}

N_list <- NULL
D_list <- NULL
Q_list <- NULL
R_list <- c(1, 2, 3)
coords <- c("N", "D", "Q")

for(coord in coords) {
  if(coord == "N") {
    N_list <- c(3, 5, 10, 20, 30, 50, 100, 250, 500, 750, 1000);
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
  }

  model_list <- c("SU", "SC", "ME", "MC", "CLM", "SVBCM", "SVBCF", "SVBUM", "SVBUF")

  percent_zero <- -1

  for (m in model_list) {
    for (n in N_list) {
      for (d in D_list) {
        for (q in Q_list) {
          for (r in R_list) {
            cat("Checking",m,"N=",n,"D=",d,"Q=",q,"R=",r,"...\n")
            if(coord == "N") {
              sweep_param_value <- n
            } else if(coord == "D") {
              sweep_param_value <- d
            } else {
              sweep_param_value <- q
            }
            destfile <- paste("simulated_data/N",n,"_D",d,"_Q",q,"_R",r,".RData", sep="")
            if(file.exists(destfile)){
              # percent zero
              load(destfile)
              percent_zero <- sum(sim_data$Y==0)/prod(dim(sim_data$Y))
            } else {
              print(paste("Error: data file",destfile,"doesn't exist!"))
              percent_zero <- -1
            }
            Lambda_true <- sim_data$Lambda_true
            destfile <- paste("fitted_models/",m,"_N",n,"_D",d,"_Q",q,"_R",r,".RData", sep="")
            if (file.exists(destfile)){
              load(destfile)
              if (m == "SU") {
                write_Stan_out(fit.su, "stan_collapsed")
                rm(fit.su)
              } else if (m == "SC") {
                write_Stan_out(fit.sc, "stan_uncollapsed")
                rm(fit.sc)
              } else if (m == "SVBCM") {
                write_Stan_out(fit.svbcm, "stan_collapsed_variationalbayes_meanfield")
                rm(fit.svbcm)
              } else if (m == "SVBCF") {
                write_Stan_out(fit.svbcf, "stan_collapsed_variationalbayes_fullrank")
                rm(fit.svbcf)
              } else if (m == "SVBUM") {
                write_Stan_out(fit.svbum, "stan_uncollapsed_variationalbayes_meanfield")
                rm(fit.svbum)
              } else if (m == "SVBUF") {
                write_Stan_out(fit.svbuf, "stan_uncollapsed_variationalbayes_fullrank")
                rm(fit.svbuf)
              } else if (m == "ME") {
                write_Mongrel_out(fit.me, "mongrel_eigen")
                rm(fit.me)
              } else if (m == "MC") {
                write_Mongrel_out(fit.mc, "mongrel_cholesky")
                rm(fit.mc)
              } else if (m == "CLM") {
                warmup_runtime <- 0
                total_runtime <- convert_to_seconds(fit.clm$metadata$total_runtime)
                lambda_RMSE <- get_Lambda_MSE(Lambda_true, fit.clm$Lambda)
                outside_CI <- get_95CI(Lambda_true, fit.clm$Lambda)
                cat(paste("conjugate_linear_model,",
                  fit.clm$metadata$mean_ess,",",
                  warmup_runtime,",",
                  total_runtime,",0,0,0,0,0,",
                  coord,",",
                  sweep_param_value,",",
                  dim(fit.clm$Lambda)[3],",",
                  0,",",
                  percent_zero,",",
                  lambda_RMSE,",",
                  outside_CI,",",
                  r,"\n",sep=""),
                  file=log_file, append=TRUE)
                rm(fit.clm)
              }
            } else {
              cat("Error: data file",destfile,"doesn't exist!\n")
            }
          }
        }
      }
    }
  }
}
