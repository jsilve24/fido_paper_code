require(ggplot2)

source("src/plotting.R")

# expects logs in current directory

exclude_su_mcp <- FALSE
file_suffix <- ""

datN <- read.csv("N_output.log")
datD <- read.csv("D_output.log")
datQ <- read.csv("Q_output.log")

# plot zeros

print("Plotting zeros...")
zero_N <- datN[datN$model == "mongrel_cholesky",]
zero_D <- datD[datD$model == "mongrel_cholesky",]
zero_Q <- datQ[datQ$model == "mongrel_cholesky",]

plot_sparsity(zero_N, "N", image_filename="results_accuracy/zeros/zeros_N_plot.png")
plot_sparsity(zero_D, "D", image_filename="results_accuracy/zeros/zeros_D_plot.png")
plot_sparsity(zero_Q, "Q", image_filename="results_accuracy/zeros/zeros_Q_plot.png")

rm(zero_N); rm(zero_D); rm(zero_Q)

#datN <- datN[datN$model != "stan_uncollapsed",]

# plot seconds per effective sample size

print("Plotting SpES...")
plot_SpES(datN, "N", image_filename="results_efficiency/SpES/SpES_N_plot.png", log_y=TRUE)
plot_SpES(datD, "D", image_filename="results_efficiency/SpES/SpES_D_plot.png", log_y=TRUE)
plot_SpES(datQ, "Q", image_filename="results_efficiency/SpES/SpES_Q_plot.png", log_y=TRUE)

# plot seconds per effective sample size

MSE_N <- datN[datN$model != "mongrel_eigen",]
MSE_D <- datD[datD$model != "mongrel_eigen",]
MSE_Q <- datQ[datQ$model != "mongrel_eigen",]

print("Plotting MSE...")
plot_accuracy(MSE_N, "N", image_filename="results_accuracy/lambda_MSE/MSE_N.png", fit_line=FALSE)
plot_accuracy(MSE_D, "D", image_filename="results_accuracy/lambda_MSE/MSE_D.png", fit_line=FALSE)
plot_accuracy(MSE_Q, "Q", image_filename="results_accuracy/lambda_MSE/MSE_Q.png", fit_line=FALSE)

# plot 95% credible interval

#print("Plotting 95 CI...")
#plot_accuracy(MSE_N, "N", use_95CI=TRUE, image_filename="results_accuracy/CI_95/95CI_N.png")
#plot_accuracy(MSE_D, "D", use_95CI=TRUE, image_filename="results_accuracy/CI_95/95CI_D.png")
#plot_accuracy(MSE_Q, "Q", use_95CI=TRUE, image_filename="results_accuracy/CI_95/95CI_Q.png")

rm(MSE_N); rm(MSE_D); rm(MSE_Q)


