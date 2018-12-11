library(ggplot2)
library(moments)
source("src/dataset_methods.R")

# use Stan collapsed as baseline and leave replicates unaveraged

args = commandArgs(trailingOnly=TRUE)
if (length(args) < 2) {
        stop(paste("Usage: Rscript summarize_3rd_moment.R {varying parameter} {output file}"))
        # Rscript summarize_3rd_moment.R N output.log
}

# need a try/catch here
vary <- args[1]
log_file <- args[2]

calc_skewness_manual <- function(Lambda) {
  d.1 <- dim(Lambda)[1]
  d.2 <- dim(Lambda)[2]
  d.3 <- dim(Lambda)[3]
  m.1 <- apply(Lambda, c(1,2), mean) # D-1 x Q
  s.1 <- apply(Lambda, c(1,2), sd) # D-1 x Q
  L.sk <- matrix(0, nrow=d.1, ncol=d.2)
  for(i in 1:d.1) {
    for(j in 1:d.2) {
      L.sk[i,j] <- mean(((Lambda[i,j,] - m.1[i,j])/s.1[i,j])**3)
    }
  }
  return(c(L.sk))
}

calc_skewness <- function(Lambda) {
  L.sk <- apply(Lambda, c(1,2), skewness)
  return(c(L.sk))
}

rmse_skewness <- function(est_Lambda_sk, ref_Lambda_sk) {
  rmse <- 0
  for (i in 1:length(ref_Lambda_sk)) {
    rmse <- rmse + (est_Lambda_sk[i] - ref_Lambda_sk[i])**2
  }
  rmse <- sqrt(rmse / length(ref_Lambda_sk))
  return(rmse)
}

cat("model,sweep_param,sweep_value,skew_RMSE\n", file=log_file)

stan_iter <- -1

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
      L.sk.full <- list()
      cat("Processing case", N, D, Q, "\n")
      for (R in 1:3) {
        load(paste("fitted_models/SC_N",N,"_D",D,"_Q",Q,"_R",R,".RData",sep=""))
        # get Stan std dev for each \lambda_{i,j}
        L.sk.full[[R]] <- calc_skewness(fit.sc$Lambda)
        cat("Baseline:",mean(L.sk.full[[R]]),"\n")
      }
      cat(paste(paste("stan_collapsed",vary,parameter_value,0,sep=","),"\n",sep=""), file=log_file, append=TRUE)
      for (R in 1:3) {
        # STAN (UNCOLLAPSED)
        tryCatch({
          load(paste("fitted_models/SU_N",N,"_D",D,"_Q",Q,"_R",R,".RData",sep=""))
          if(stan_iter < 0) { stan_iter <- dim(fit.su$Lambda)[3] }
          rmse <- rmse_skewness(calc_skewness(fit.su$Lambda), L.sk.full[[R]])
          cat(paste(paste("stan_uncollapsed",vary,parameter_value,rmse,sep=","),"\n",sep=""), file=log_file, append=TRUE)
        }, error = function(e) { })

        # STAN (COLLAPSED) VARIATIONAL BAYES MEAN FIELD
        tryCatch({
          load(paste("fitted_models/SVBCM_N",N,"_D",D,"_Q",Q,"_R",R,".RData",sep=""))
          rmse <- rmse_skewness(calc_skewness(fit.svbcm$Lambda), L.sk.full[[R]])
          cat(paste(paste("stan_collapsed_variationalbayes_meanfield",vary,parameter_value,rmse,sep=","),"\n",sep=""), file=log_file, append=TRUE)
        }, error = function(e) { })

        # STAN (COLLAPSED) VARIATIONAL BAYES FULL RANK
        tryCatch({
          load(paste("fitted_models/SVBCF_N",N,"_D",D,"_Q",Q,"_R",R,".RData",sep=""))
          rmse <- rmse_skewness(calc_skewness(fit.svbcf$Lambda), L.sk.full[[R]])
          cat(paste(paste("stan_collapsed_variationalbayes_fullrank",vary,parameter_value,rmse,sep=","),"\n",sep=""), file=log_file, append=TRUE)
        }, error = function(e) { })

        # STAN (UNCOLLAPSED) VARIATIONAL BAYES MEAN FIELD
        tryCatch({
          load(paste("fitted_models/SVBUM_N",N,"_D",D,"_Q",Q,"_R",R,".RData",sep=""))
          rmse <- rmse_skewness(calc_skewness(fit.svbum$Lambda), L.sk.full[[R]])
          cat(paste(paste("stan_uncollapsed_variationalbayes_meanfield",vary,parameter_value,rmse,sep=","),"\n",sep=""), file=log_file, append=TRUE)
        }, error = function(e) { })

        # STAN (UNCOLLAPSED) VARIATIONAL BAYES FULL RANK
        tryCatch({
          load(paste("fitted_models/SVBUF_N",N,"_D",D,"_Q",Q,"_R",R,".RData",sep=""))
          rmse <- rmse_skewness(calc_skewness(fit.svbuf$Lambda), L.sk.full[[R]])
          cat(paste(paste("stan_uncollapsed_variationalbayes_fullrank",vary,parameter_value,rmse,sep=","),"\n",sep=""), file=log_file, append=TRUE)
        }, error = function(e) { })

        # MONGREL CHOLESKY
        tryCatch({
          load(paste("fitted_models/MC_N",N,"_D",D,"_Q",Q,"_R",R,".RData",sep=""))
          L.sk <- calc_skewness(fit.mc$Lambda[,,1:stan_iter])
          rmse <- rmse_skewness(L.sk, L.sk.full[[R]])
          cat("MC (mean):",mean(L.sk),"\n")
          cat("MC (max):",max(L.sk),"\n")
          cat(paste(paste("mongrel_cholesky",vary,parameter_value,rmse,sep=","),"\n",sep=""), file=log_file, append=TRUE)
        }, error = function(e) { })

        # CONJUGATE LINEAR MODEL
        tryCatch({
          load(paste("fitted_models/CLM_N",N,"_D",D,"_Q",Q,"_R",R,".RData",sep=""))
          rmse <- rmse_skewness(calc_skewness(fit.clm$Lambda), L.sk.full[[R]])
          cat(paste(paste("conjugate_linear_model",vary,parameter_value,rmse,sep=","),"\n",sep=""), file=log_file, append=TRUE)
        }, error = function(e) { })
      }
    }
  }
}
