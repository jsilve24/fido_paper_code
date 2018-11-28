library(ggplot2)
source("src/dataset_methods.R")

Stan_iter <- 2000
use_log <- FALSE

# use Stan collapsed as baseline and leave replicates unaveraged

log_file <- "second_moment_data.log"
cat("model,sweep_param,sweep_value,sd_MSE\n", file=log_file)

varying <- c("N", "D", "Q")
for(vary in varying) {
  if (vary == 'N') {
    # Stan (collapsed) fails to finish at N=1000; omit
    N_list <- c(3, 5, 10, 20, 30, 50, 100, 250, 500, 750)
    D_list <- c(30)
    Q_list <- c(5)
  } else if(vary == 'D') {
    N_list <- c(100)
    D_list <- c(3, 5, 10, 25, 50, 75, 100, 250, 500)
    Q_list <- c(5)
  } else if(vary == 'Q') {
    N_list <- c(100)
    D_list <- c(30)
    Q_list <- c(2, 4, 10, 20, 50, 75, 100, 250, 500)
  }

  for (N in N_list) {
    for (D in D_list) {
      for (Q in Q_list) {
      	parameter_value <- N
      	if(vary == "D") {
      		parameter_value <- D
      	}
      	if(vary == "Q") {
      		parameter_value <- Q
      	}
        L.sd.full <- list()
        cat("Processing case", N, D, Q, "\n")
        for (R in 1:3) {
          load(paste("fitted_models/SC_N",N,"_D",D,"_Q",Q,"_R",R,".RData",sep=""))
          # get Stan std dev for each \lambda_{i,j}
          L.sd <- apply(fit.sc$Lambda, c(1,2), sd)
          L.sd <- c(L.sd) # column-major
          L.sd.full[[R]] <- L.sd
        }
        cat(paste(paste("stan_collapsed",vary,parameter_value,0,sep=","),"\n",sep=""), file=log_file, append=TRUE)
        for (R in 1:3) {
          # STAN (UNCOLLAPSED)
          tryCatch({
            load(paste("fitted_models/SU_N",N,"_D",D,"_Q",Q,"_R",R,".RData",sep=""))
            mse.su <- 0
            L.sd.su <- apply(fit.su$Lambda, c(1,2), sd)
            L.sd.su <- c(L.sd.su)
            for (i in 1:length(L.sd.full[[R]])) {
              mse.su <- mse.su + (L.sd.su[i] - L.sd.full[[R]][i])**2
            }
            mse.su <- mse.su / length(L.sd.full[[R]])
            cat(paste(paste("stan_uncollapsed",vary,parameter_value,mse.su,sep=","),"\n",sep=""), file=log_file, append=TRUE)
          }, error = function(e) { })

          # MONGREL CHOLESKY
          tryCatch({
            load(paste("fitted_models/MC_N",N,"_D",D,"_Q",Q,"_R",R,".RData",sep=""))
            mse.mc <- 0
            L.sd.mc <- apply(fit.mc$Lambda[,,1:Stan_iter], c(1,2), sd)
            L.sd.mc <- c(L.sd.mc)
            for (i in 1:length(L.sd.full[[R]])) {
              mse.mc <- mse.mc + (L.sd.mc[i] - L.sd.full[[R]][i])**2
            }
            mse.mc <- mse.mc / length(L.sd.full[[R]])
            cat(paste(paste("mongrel_cholesky",vary,parameter_value,mse.mc,sep=","),"\n",sep=""), file=log_file, append=TRUE)
          }, error = function(e) { })

          # CONJUGATE LINEAR MODEL
          tryCatch({
            load(paste("fitted_models/CLM_N",N,"_D",D,"_Q",Q,"_R",R,".RData",sep=""))
            mse.clm <- 0
            L.sd.clm <- apply(fit.clm$Lambda, c(1,2), sd)
            L.sd.clm <- c(L.sd.clm)
            for (i in 1:length(L.sd.full[[R]])) {
              mse.clm <- mse.clm + (L.sd.clm[i] - L.sd.full[[R]][i])**2
            }
            mse.clm <- mse.clm / length(L.sd.full[[R]])
            cat(paste(paste("conjugate_linear_model",vary,parameter_value,mse.clm,sep=","),"\n",sep=""), file=log_file, append=TRUE)
          }, error = function(e) { })
        }
      }
    }
  }
}
