library(ggplot2)
source("src/dataset_methods.R")

args = commandArgs(trailingOnly=TRUE)
if (length(args) < 1) {
	stop(paste("Usage: Rscript plot_2nd_moment.R {dimension}"))
        # Rscript plot_2nd_moment.R N
}
vary <- args[1]
use_log <- FALSE
combine_replicates <- FALSE

destfile <- paste("mom2_varying-",vary,".RData",sep="")
if(combine_replicates) {
	destfile <- paste("mom2_varying-",vary,"_combined_reps.RData",sep="")
}

if(!file.exists(destfile)){
	if (vary == 'N') {
		#N_list <- c(10, 20, 30, 50, 100, 250, 500, 750, 1000) # full
		N_list <- c(10, 20, 30, 50, 100, 250, 500, 750)
		D_list <- c(30)
		Q_list <- c(5)
	} else if(vary == 'D') {
		N_list <- c(100)
		#D_list <- c(3, 5, 10, 25, 50, 75, 100, 250, 500) # full
		D_list <- c(3, 5, 10, 25, 50, 75, 100)
		Q_list <- c(5)
	} else if(vary == 'Q') {
		N_list <- c(100)
		D_list <- c(30)
		#Q_list <- c(2, 4, 10, 20, 50, 75, 100, 250, 500) # full
		Q_list <- c(2, 4, 10, 20, 50, 75, 100, 250, 500)
	}

	no_models <- 3

	data <- data.frame(N=numeric(), D=numeric(), Q=numeric(), model=character(), value=numeric(), stringsAsFactors=FALSE)

	for (N in N_list) {
		for (D in D_list) {
			for (Q in Q_list) {
				L.sd.full <- list()
				cat(N, D, Q, "-- Stan baseline\n")
				for (R in 1:3) {
					load(paste("fitted_models/SU_N",N,"_D",D,"_Q",Q,"_R",R,".RData",sep=""))
					# get Stan std dev for each \lambda_{i,j}
					L.sd <- apply(fit.su$Lambda, c(1,2), sd)
					L.sd <- c(L.sd) # column-major
					L.sd.full[[R]] <- L.sd
				}
				mse.mc <- 0
				mse.sc <- 0
				mse.clm <- 0
				for (R in 1:3) {
					#load(paste("simulated_data/N",N,"_D",D,"_Q",Q,"_R",R,".RData",sep=""))
					fit.mc <- NULL; fit.sc <- NULL; fit.clm <- NULL
					tryCatch({
						load(paste("fitted_models/SC_N",N,"_D",D,"_Q",Q,"_R",R,".RData",sep=""))
					}, error = function(e) { })
					tryCatch({
						load(paste("fitted_models/MC_N",N,"_D",D,"_Q",Q,"_R",R,".RData",sep=""))
					}, error = function(e) { })
					tryCatch({
						load(paste("fitted_models/CLM_N",N,"_D",D,"_Q",Q,"_R",R,".RData",sep=""))
					}, error = function(e) { })
					cat(N, D, Q, R, "-- Mongrel, CLM, SC\n")

					if(!is.null(fit.sc)) {
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
						L.sd.mc <- apply(fit.mc$Lambda, c(1,2), sd)
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
						mse.sc <- mse.sc / sum(unlist(lapply(L.sd.full, length)))
						mse.mc <- mse.mc / sum(unlist(lapply(L.sd.full, length)))
						mse.clm <- mse.clm / sum(unlist(lapply(L.sd.full, length)))
					}

					if(!combine_replicates) {
						mse.sc <- mse.sc / length(L.sd.full[[R]])
						mse.mc <- mse.mc / length(L.sd.full[[R]])
						mse.clm <- mse.clm / length(L.sd.full[[R]])

						data[nrow(data)+1,] <- list(N, D, Q, "stan_uncollapsed", 0)
						if(!is.null(fit.sc)) {
							data[nrow(data)+1,] <- list(N, D, Q, "stan_collapsed", mse.sc)
						}
						if(!is.null(fit.mc)) {
							data[nrow(data)+1,] <- list(N, D, Q, "mongrel_cholesky", mse.mc)
						}
						if(!is.null(fit.clm)) {
							data[nrow(data)+1,] <- list(N, D, Q, "conjugate_linear_model", mse.clm)
						}
					}

				}
				if(combine_replicates) {
					data[nrow(data)+1,] <- list(N, D, Q, "stan_uncollapsed", 0)
					if(!is.null(fit.sc)) {
						data[nrow(data)+1,] <- list(N, D, Q, "stan_collapsed", mse.sc)
					}
					if(!is.null(fit.mc)) {
						data[nrow(data)+1,] <- list(N, D, Q, "mongrel_cholesky", mse.mc)
					}
					if(!is.null(fit.clm)) {
						data[nrow(data)+1,] <- list(N, D, Q, "conjugate_linear_model", mse.clm)
					}
				}
			}
		}
	}
	if(combine_replicates) {
		save(data, file=paste("mom2_varying-",vary,"_combined_reps.RData",sep=""))
	} else {
		save(data, file=paste("mom2_varying-",vary,".RData",sep=""))
	}
} else {
	load(destfile)
}

p <- ggplot(data, aes_string(x=vary, y="value", color="model")) + geom_point() + scale_x_log10()
if(use_log) {
	p <- p + scale_y_log10()
}

if(combine_replicates) {
	ggsave(paste("mom2_varying-",vary,"_combined_reps.png",sep=""), plot=last_plot(), width=10, height=10, dpi=300)
} else {
	ggsave(paste("mom2_varying-",vary,".png",sep=""), plot=last_plot(), width=10, height=10, dpi=300)
}
