library(ggplot2)
source("src/dataset_methods.R")

args = commandArgs(trailingOnly=TRUE)
if (length(args) < 1) {
	stop(paste("Usage: Rscript plot_2nd_moment.R {dimension}"))
        # Rscript plot_2nd_moment.R N
}
vary <- args[1]
use_log <- FALSE

destfile <- paste("mom2_varying-",vary,".RData",sep="")
if(!file.exists(destfile)){

	if (vary == 'N') {
		N_list <- c(10, 20, 30, 50, 100, 250, 500, 750) # omitting 1000
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

	no_models <- 3

	data <- data.frame(N=numeric(), D=numeric(), Q=numeric(), model=character(), R=numeric(), value=numeric(), stringsAsFactors=FALSE)

	for (N in N_list) {
		for (D in D_list) {
			for (Q in Q_list) {
				for (R in 1:3) {
					load(paste("simulated_data/N",N,"_D",D,"_Q",Q,"_R",R,".RData",sep=""))
					load(paste("fitted_models/SC_N",N,"_D",D,"_Q",Q,"_R",R,".RData",sep=""))
					load(paste("fitted_models/MC_N",N,"_D",D,"_Q",Q,"_R",R,".RData",sep=""))
					load(paste("fitted_models/CLM_N",N,"_D",D,"_Q",Q,"_R",R,".RData",sep=""))
					cat(N, D, Q, R, "\n")

					# get Stan std dev for each \lambda_{i,j}
					L.sd <- apply(fit.sc$Lambda, c(1,2), sd)
					L.sd <- c(L.sd) # column-major

					# get mean sq dev of MC and CLM \lambda_{i,j} from Stan baseline
					L.sd.mc <- apply(fit.mc$Lambda, c(1,2), sd)
					L.sd.mc <- c(L.sd.mc)
					L.sd.clm <- apply(fit.clm$Lambda, c(1,2), sd)
					L.sd.clm <- c(L.sd.clm)

					mse.mc <- 0
					for (i in 1:length(L.sd)) {
						mse.mc <- mse.mc + (L.sd.mc[i] - L.sd[i])**2
					}
					mse.mc <- mse.mc / length(L.sd)
					mse.clm <- 0
					for (i in 1:length(L.sd)) {
						mse.clm <- mse.clm + (L.sd.clm[i] - L.sd[i])**2
					}
					mse.clm <- mse.clm / length(L.sd)

					data[nrow(data)+1,] <- list(N, D, Q, "SC", R, 0)
					data[nrow(data)+1,] <- list(N, D, Q, "MC", R, mse.mc)
					data[nrow(data)+1,] <- list(N, D, Q, "CLM", R, mse.clm)
				}
			}
		}
	}
	save(data, file=paste("mom2_varying-",vary,".RData",sep=""))
} else {
	load(destfile)
}

p <- ggplot(data, aes_string(x=vary, y="value", color="model")) + geom_point() + scale_x_log10()
if(use_log) {
	p <- p + scale_y_log10()
}

ggsave(paste("mom2_varying-",vary,".png",sep=""), plot=last_plot(), width=10, height=10, dpi=300)
