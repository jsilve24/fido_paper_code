library(ggplot2)
source("src/dataset_methods.R")

Stan_iter <- 2000
use_log <- FALSE
combine_replicates <- FALSE
SU_as_baseline <- FALSE
# TRUE : use Stan (uncollapsed) as our baseline/reference measurement for "true" standard deviation
# FALSE : use Stan (collapsed) as our baseline; may want to use because Stan (uncollapsed) fails to finish in
#    many cases, meaning we have no baseline!

log_file <- "second_moment_data.log"
cat("model,sweep_param,sweep_value,sd_MSE\n", file=log_file)

varying <- c("N", "D", "Q")
for(vary in varying) {
  if (vary == 'N') {
    if(SU_as_baseline) {
      N_list <- c(10, 20, 30, 50, 100, 250, 500, 750, 1000)
    } else {
      # Stan (collapsed) failed to finish at N=1000; omit
      N_list <- c(10, 20, 30, 50, 100, 250, 500, 750)
    }
    D_list <- c(30)
    Q_list <- c(5)
  } else if(vary == 'D') {
    N_list <- c(100)
    if(SU_as_baseline) {
      # Stan (uncollapsed) failed to finish at D=250
      D_list <- c(3, 5, 10, 25, 50, 75, 100)
    } else {
      D_list <- c(3, 5, 10, 25, 50, 75, 100, 250, 500)
    }
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
        cat("Processing case", N, D, Q, "for Stan baseline\n")
        for (R in 1:3) {
          if(SU_as_baseline) {
            load(paste("fitted_models/SU_N",N,"_D",D,"_Q",Q,"_R",R,".RData",sep=""))
            # get Stan std dev for each \lambda_{i,j}
            L.sd <- apply(fit.su$Lambda, c(1,2), sd)
            L.sd <- c(L.sd) # column-major
            L.sd.full[[R]] <- L.sd
          } else {
            load(paste("fitted_models/SC_N",N,"_D",D,"_Q",Q,"_R",R,".RData",sep=""))
            # get Stan std dev for each \lambda_{i,j}
            L.sd <- apply(fit.sc$Lambda, c(1,2), sd)
            L.sd <- c(L.sd) # column-major
            L.sd.full[[R]] <- L.sd
          }
        }
        mse.mc <- 0
        mse.sc <- 0
        mse.clm <- 0
        for (R in 1:3) {
          fit.mc <- NULL; fit.sc <- NULL; fit.clm <- NULL
          # tryCatch in case some of these cases failed to produce output
          if(SU_as_baseline) {
            tryCatch({
              load(paste("fitted_models/SC_N",N,"_D",D,"_Q",Q,"_R",R,".RData",sep=""))
            }, error = function(e) { })
          }
          tryCatch({
            load(paste("fitted_models/MC_N",N,"_D",D,"_Q",Q,"_R",R,".RData",sep=""))
          }, error = function(e) { })
          tryCatch({
            load(paste("fitted_models/CLM_N",N,"_D",D,"_Q",Q,"_R",R,".RData",sep=""))
          }, error = function(e) { })
          cat("Processing case", N, D, Q, R, "for other models\n")

          if(SU_as_baseline && !is.null(fit.sc)) {
            # STAN (COLLAPSED)
            L.sd.sc <- apply(fit.sc$Lambda, c(1,2), sd)
            L.sd.sc <- c(L.sd.sc)
            for (i in 1:length(L.sd.full[[R]])) {
              mse.sc <- mse.sc + (L.sd.sc[i] - L.sd.full[[R]][i])**2
            }
          }

          if(!is.null(fit.mc)) {
            # get mean sq dev of \lambda_{i,j} from Stan baseline
            # MONGREL CHOLESKY
            # truncate Mongrel samples to same # as Stan
            L.sd.mc <- apply(fit.mc$Lambda[,,1:Stan_iter], c(1,2), sd)
            L.sd.mc <- c(L.sd.mc)
            for (i in 1:length(L.sd.full[[R]])) {
              mse.mc <- mse.mc + (L.sd.mc[i] - L.sd.full[[R]][i])**2
            }
          }

          if(!is.null(fit.clm)) {
            # CONJUGATE LINEAR MODEL
            L.sd.clm <- apply(fit.clm$Lambda, c(1,2), sd)
            L.sd.clm <- c(L.sd.clm)
            for (i in 1:length(L.sd.full[[R]])) {
              mse.clm <- mse.clm + (L.sd.clm[i] - L.sd.full[[R]][i])**2
            }
          }

          if(combine_replicates && R == 3) {
            if(SU_as_baseline) {
              mse.sc <- mse.sc / sum(unlist(lapply(L.sd.full, length)))
            }
            mse.mc <- mse.mc / sum(unlist(lapply(L.sd.full, length)))
            mse.clm <- mse.clm / sum(unlist(lapply(L.sd.full, length)))
          }

          if(!combine_replicates) {
            if(SU_as_baseline) {
              mse.sc <- mse.sc / length(L.sd.full[[R]])
            }
            mse.mc <- mse.mc / length(L.sd.full[[R]])
            mse.clm <- mse.clm / length(L.sd.full[[R]])

            if(SU_as_baseline) {
              cat(paste(paste("stan_uncollapsed",vary,parameter_value,0,sep=","),"\n",sep=""), file=log_file, append=TRUE)
              if(!is.null(fit.sc)) {
                cat(paste(paste("stan_collapsed",vary,parameter_value,mse.sc,sep=","),"\n",sep=""), file=log_file, append=TRUE)
              }
            } else {
              cat(paste(paste("stan_collapsed",vary,parameter_value,0,sep=","),"\n",sep=""), file=log_file, append=TRUE)
            }
            if(!is.null(fit.mc)) {
            	cat(paste(paste("mongrel_cholesky",vary,parameter_value,mse.mc,sep=","),"\n",sep=""), file=log_file, append=TRUE)
            }
            if(!is.null(fit.clm)) {
            	cat(paste(paste("conjugate_linear_model",vary,parameter_value,mse.clm,sep=","),"\n",sep=""), file=log_file, append=TRUE)
            }
          }

        }
        if(combine_replicates) {
          if(SU_as_baseline) {
          	cat(paste(paste("stan_uncollapsed",vary,parameter_value,0,sep=","),"\n",sep=""), file=log_file, append=TRUE)
            if(!is.null(fit.sc)) {
            	cat(paste(paste("stan_collapsed",vary,parameter_value,mse.sc,sep=","),"\n",sep=""), file=log_file, append=TRUE)
            }
          } else {
          	cat(paste(paste("stan_collapsed",vary,parameter_value,0,sep=","),"\n",sep=""), file=log_file, append=TRUE)
          }
          if(!is.null(fit.mc)) {
          	cat(paste(paste("mongrel_cholesky",vary,parameter_value,mse.mc,sep=","),"\n",sep=""), file=log_file, append=TRUE)
          }
          if(!is.null(fit.clm)) {
          	cat(paste(paste("conjugate_linear_model",vary,parameter_value,mse.clm,sep=","),"\n",sep=""), file=log_file, append=TRUE)
          }
        }
      }
    }
  }
}

#p <- ggplot(data, aes_string(x=vary, y="value", color="model")) + geom_point() + scale_x_log10()
#if(use_log) {
#  p <- p + scale_y_log10()
#}

#if(combine_replicates) {
#  ggsave(paste("mom2_varying-",vary,"_combined_reps.png",sep=""), plot=last_plot(), width=10, height=10, dpi=300)
#} else {
#  ggsave(paste("mom2_varying-",vary,".png",sep=""), plot=last_plot(), width=10, height=10, dpi=300)
#}
