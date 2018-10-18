require(purrr)
require(driver)
require(dplyr)
require(ggstance)
require(ggplot2)

source("src/dataset_methods.R")
source("src/plotting.R")

load("fit_MC_N100_D30_Q100.RData")
load("fit_ME_N100_D30_Q100.RData")
load("fit_SC_N100_D30_Q100.RData")

N_ll <- 1
N_ul <- 5 # max 100

D_ll <- 1
D_ul <- 5 # max 30

Q_ll <- 1
Q_ul <- 10 # max 100

print(paste("N range",N_ll,":",N_ul))
print(paste("D range",D_ll,":",D_ul))
print(paste("Q range",Q_ll,":",Q_ul))

do_plot <- TRUE

# should all be the same and by quick crappy visual inspection they are
#Lambda_true <- fit.mongrel.cholesky$mdataset$Lambda_true[D_ll:D_ul, Q_ll:Q_ul]
#Lambda_true <- fit.mongrel.eigen$mdataset$Lambda_true[D_ll:D_ul, Q_ll:Q_ul]
Lambda_true <- fit.sc$mdataset$Lambda_true[D_ll:D_ul, Q_ll:Q_ul]

MC <- fit.mongrel.cholesky
ME <- fit.mongrel.eigen
SC <- fit.sc

# truncate the Lambdas and etas; we'll just look at a subset
MC$Lambda <- MC$Lambda[D_ll:D_ul, Q_ll:Q_ul,]
ME$Lambda <- ME$Lambda[D_ll:D_ul, Q_ll:Q_ul,]
SC$Lambda <- SC$Lambda[D_ll:D_ul, Q_ll:Q_ul,]

MC$Eta <- MC$Eta[D_ll:D_ul, N_ll:N_ul,]
ME$Eta <- ME$Eta[D_ll:D_ul, N_ll:N_ul,]
SC$Eta <- SC$Eta[D_ll:D_ul, N_ll:N_ul,]

mfits <- list("fit.mc"=MC, "fit.me"=ME, "fit.sc"=SC)

# calculate MSE
ref <- c(Lambda_true)

get_MSE <- function(est_Lambda, trunc_d1, trunc_d2, d3, ref) {
  dim(est_Lambda) <- c(trunc_d1*trunc_d2,d3)
  return(mean(apply(est_Lambda, 2, sample_SE, y=ref)))
}

print(paste("MSE MC:",get_MSE(MC$Lambda, nrow(MC$Lambda), ncol(MC$Lambda), MC$iter, ref)))
print(paste("MSE ME:",get_MSE(ME$Lambda, nrow(ME$Lambda), ncol(ME$Lambda), ME$iter, ref)))
print(paste("MSE SC:",get_MSE(SC$Lambda, nrow(SC$Lambda), ncol(SC$Lambda), SC$iter, ref)))

# get count of Lambdas outside 95% CI
get_outside_count <- function(est_Lambda, ref) {
  intervals <- apply(est_Lambda, c(1,2), function(x) quantile(x, probs=c(0.025, 0.975)))
  lower <- c(intervals[1,,])
  upper <- c(intervals[2,,])
  sum_outside <- sum(apply(rbind(lower, upper, ref), 2, lambda_outside_bounds))
  total_L <- nrow(est_Lambda)*ncol(est_Lambda)
  return(sum_outside)
}

print(paste("Outside 95 CI count MC:",get_outside_count(MC$Lambda, ref)))
print(paste("Outside 95 CI count ME:",get_outside_count(ME$Lambda, ref)))
print(paste("Outside 95 CI count SC:",get_outside_count(SC$Lambda, ref)))

if(do_plot) {
  plot_lambda(mfits, Lambda_true=Lambda_true, image_filename="Lambda_subset_plot_N100_D30_Q100.png")
  plot_eta(mfits, Eta_true=NULL, image_filename="Eta_subset_plot_N100_D30_Q100.png")
}








