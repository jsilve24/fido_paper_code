require(purrr)
require(driver)
require(dplyr)
require(ggstance)
require(ggplot2)

source("src/dataset_methods.R")
source("src/plotting.R")

load("fitted_models/CLM_N30_D30_Q5_R1.RData")
load("fitted_models/ME_N30_D30_Q5_R1.RData")

#output_dir <- "results_accuracy/run4/"
output_dir <- ""

N_ll <- 1
N_ul <- 5

D_ll <- 1
D_ul <- 29

Q_ll <- 1
Q_ul <- 5

print(paste("N range",N_ll,":",N_ul))
print(paste("D range",D_ll,":",D_ul))
print(paste("Q range",Q_ll,":",Q_ul))

do_plot <- TRUE

# should all be the same for all models
Lambda_true <- fit.me$mdataset$Lambda_true[D_ll:D_ul, Q_ll:Q_ul]

ME <- fit.me
CLM <- fit.clm

# truncate the Lambdas and etas; we'll just look at a subset
ME$Lambda <- ME$Lambda[D_ll:D_ul, Q_ll:Q_ul,]
CLM$Lambda <- CLM$Lambda[D_ll:D_ul, Q_ll:Q_ul,]

mfits <- list("fit.me"=ME, "fit.clm"=CLM)

# calculate MSE
ref <- c(Lambda_true)

get_MSE <- function(est_Lambda, trunc_d1, trunc_d2, d3, ref) {
  dim(est_Lambda) <- c(trunc_d1*trunc_d2,d3)
  return(mean(apply(est_Lambda, 2, sample_SE, y=ref)))
}

cat(paste("MSE ME:",get_MSE(ME$Lambda, nrow(ME$Lambda), ncol(ME$Lambda), ME$iter, ref),"\n"))
cat(paste("MSE CLM:",get_MSE(CLM$Lambda, nrow(CLM$Lambda), ncol(CLM$Lambda), ME$iter, ref),"\n"))

# get count of Lambdas outside 95% CI
get_outside_count <- function(est_Lambda, ref) {
  intervals <- apply(est_Lambda, c(1,2), function(x) quantile(x, probs=c(0.025, 0.975)))
  lower <- c(intervals[1,,])
  upper <- c(intervals[2,,])
  sum_outside <- sum(apply(rbind(lower, upper, ref), 2, lambda_outside_bounds))
  total_L <- nrow(est_Lambda)*ncol(est_Lambda)
  return(sum_outside)
}

cat(paste("Outside 95 CI count ME:",get_outside_count(ME$Lambda, ref),"\n"))
cat(paste("Outside 95 CI count CLM:",get_outside_count(CLM$Lambda, ref),"\n"))

if(do_plot) {
  plot_lambda(mfits, Lambda_true=Lambda_true, image_filename=paste(output_dir,"Lambda_subset_plot_N30_D30_Q5_R1.png",sep=""))
}








