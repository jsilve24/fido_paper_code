library(ggplot2)
source("src/dataset_methods.R")

# use Stan collapsed as baseline and leave replicates unaveraged

args = commandArgs(trailingOnly=TRUE)
if (length(args) < 2) {
        stop(paste("Usage: Rscript summarize_2nd_moment.R {varying parameter} {output file}"))
        # Rscript summarize_2nd_moment.R N output.log
}

# need a try/catch here
vary <- args[1]
log_file <- args[2]

calc_msdod <- function(est_Lambda, ref_Lambda) {
  rmse <- 0
  L.sd <- apply(est_Lambda, c(1,2), sd)
  L.sd <- c(L.sd)
  for (i in 1:length(ref_Lambda)) {
    rmse <- rmse + (L.sd[i] - ref_Lambda[i])**2
  }
  rmse <- sqrt(rmse / length(ref_Lambda))
  return(rmse)
}

cat("model,sweep_param,sweep_value,sd_RMSE\n", file=log_file)

stan_iter <- -1

if(vary == "N") {
  N_list <- c(3, 5, 10, 20, 30, 50, 100, 250, 500, 750); # omit N=1000 since nothing to compare to
  D_list <- c(30);
  Q_list <- c(5);
} else if(vary == "D") {
  N_list <- c(100);
  D_list <- c(3, 5, 10, 25, 50, 75, 100, 250, 500);
  Q_list <- c(5);
} else if(vary == "Q") {
  N_list <- c(100);
  D_list <- c(30);
  Q_list <- c(2, 4, 10, 20, 50, 75, 100, 250, 500);
}

model_dir <- "fitted_models_2018-12-30"
MKL_suffix <- "_MKL_4"

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
      for (R in 1:3) {
        destfile <- paste(model_dir,"/SC_N",N,"_D",D,"_Q",Q,"_R",R,MKL_suffix,".RData",sep="")
	cat("Loading data file '",destfile,"'\n")
        load(destfile)
        # get Stan std dev for each \lambda_{i,j}
        L.sd <- apply(fit.sc$Lambda, c(1,2), sd)
        L.sd <- c(L.sd) # column-major
        L.sd.full[[R]] <- L.sd
      }
      cat(paste(paste("stan_collapsed",vary,parameter_value,0,sep=","),"\n",sep=""), file=log_file, append=TRUE)
      for (R in 1:3) {
        # STAN (UNCOLLAPSED)
        tryCatch({
          destfile <- paste(model_dir,"/SU_N",N,"_D",D,"_Q",Q,"_R",R,MKL_suffix,".RData",sep="")
          cat("Loading fit file '",destfile,"'\n")
          load(destfile)
          if(stan_iter < 0) { stan_iter <- dim(fit.su$Lambda)[3] }
          rmse <- calc_msdod(fit.su$Lambda, L.sd.full[[R]])
          cat(paste(paste("stan_uncollapsed",vary,parameter_value,rmse,sep=","),"\n",sep=""), file=log_file, append=TRUE)
        }, error = function(e) { })

        # STAN (COLLAPSED) VARIATIONAL BAYES MEAN FIELD
        tryCatch({
          destfile <- paste(model_dir,"/SVBCM_N",N,"_D",D,"_Q",Q,"_R",R,MKL_suffix,".RData",sep="")
          cat("Loading fit file '",destfile,"'\n")
          load(destfile)
          rmse <- calc_msdod(fit.svbcm$Lambda, L.sd.full[[R]])
          cat(paste(paste("stan_collapsed_variationalbayes_meanfield",vary,parameter_value,rmse,sep=","),"\n",sep=""), file=log_file, append=TRUE)
        }, error = function(e) { })

        # STAN (COLLAPSED) VARIATIONAL BAYES FULL RANK
        #tryCatch({
        #  destfile <- paste(model_dir,"/SVBCF_N",N,"_D",D,"_Q",Q,"_R",R,MKL_suffix,".RData",sep="")
        #  cat("Loading fit file '",destfile,"'\n")
        #  load(destfile)
        #  rmse <- calc_msdod(fit.svbcf$Lambda, L.sd.full[[R]])
        #  cat(paste(paste("stan_collapsed_variationalbayes_fullrank",vary,parameter_value,rmse,sep=","),"\n",sep=""), file=log_file, append=TRUE)
        #}, error = function(e) { })

        # STAN (UNCOLLAPSED) VARIATIONAL BAYES MEAN FIELD
        #tryCatch({
        #  destfile <- paste(model_dir,"/SVBUM_N",N,"_D",D,"_Q",Q,"_R",R,MKL_suffix,".RData",sep="")
        #  cat("Loading fit file '",destfile,"'\n")
        #  load(destfile)
        #  rmse <- calc_msdod(fit.svbum$Lambda, L.sd.full[[R]])
        #  cat(paste(paste("stan_uncollapsed_variationalbayes_meanfield",vary,parameter_value,rmse,sep=","),"\n",sep=""), file=log_file, append=TRUE)
        #}, error = function(e) { })

        # STAN (UNCOLLAPSED) VARIATIONAL BAYES FULL RANK
        #tryCatch({
        #  destfile <- paste(model_dir,"/SVBUF_N",N,"_D",D,"_Q",Q,"_R",R,MKL_suffix,".RData",sep="")
        #  cat("Loading fit file '",destfile,"'\n")
        #  load(destfile)
        #  rmse <- calc_msdod(fit.svbuf$Lambda, L.sd.full[[R]])
        #  cat(paste(paste("stan_uncollapsed_variationalbayes_fullrank",vary,parameter_value,rmse,sep=","),"\n",sep=""), file=log_file, append=TRUE)
        #}, error = function(e) { })

        # MONGREL CHOLESKY
        tryCatch({
          destfile <- paste(model_dir,"/MC_N",N,"_D",D,"_Q",Q,"_R",R,MKL_suffix,".RData",sep="")
          cat("Loading fit file '",destfile,"'\n")
          load(destfile)
          rmse <- calc_msdod(fit.mc$Lambda[,,1:stan_iter], L.sd.full[[R]])
          cat(paste(paste("mongrel_cholesky",vary,parameter_value,rmse,sep=","),"\n",sep=""), file=log_file, append=TRUE)
        }, error = function(e) { })

        # CONJUGATE LINEAR MODEL
        tryCatch({
          destfile <- paste(model_dir,"/CLM_N",N,"_D",D,"_Q",Q,"_R",R,MKL_suffix,".RData",sep="")
          cat("Loading fit file '",destfile,"'\n")
          load(destfile)
          rmse <- calc_msdod(fit.clm$Lambda, L.sd.full[[R]])
          cat(paste(paste("conjugate_linear_model",vary,parameter_value,rmse,sep=","),"\n",sep=""), file=log_file, append=TRUE)
        }, error = function(e) { })
      }
    }
  }
}
