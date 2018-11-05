source("src/dataset_methods.R")

N <- 100
#D_list <- c(3, 5, 10, 25, 50, 75, 100, 250) # 500 omitted b/c Mongrel eigen jobs don't finish
D_list <- c(3, 5, 10) # testing
Q <- 5

get_lambda_mse <- function(mfit) {
	est_Lambda <- mfit$Lambda
	dim(est_Lambda) <- c(nrow(mfit$Lambda)*ncol(mfit$Lambda),2000L)
	lambda_MSE <- mean(apply(est_Lambda, 2, sample_SE, y=ref))
	lambda_MSE
}

for (D in D_list) {
	var.sc <- 0
	var.me <- 0
	var.clm <- 0
	for (R in 1:3) {
		load(paste("simulated_data/N",N,"_D",D,"_Q",Q,"_R",R,".RData",sep=""))
		load(paste("fitted_models/SC_N",N,"_D",D,"_Q",Q,"_R",R,".RData",sep=""))
		load(paste("fitted_models/ME_N",N,"_D",D,"_Q",Q,"_R",R,".RData",sep=""))
		load(paste("fitted_models/CLM_N",N,"_D",D,"_Q",Q,"_R",R,".RData",sep=""))
		cat(N, D, Q, R, "\n")

		# get Lambda MSE
		ref <- c(sim_data$Lambda_true)
		cat("\tMSE SC:",get_lambda_mse(fit.sc),"\n")
		cat("\tMSE ME:",get_lambda_mse(fit.me),"\n")
		cat("\tMSE CLM:",get_lambda_mse(fit.clm),"\n\n")
	}
}

