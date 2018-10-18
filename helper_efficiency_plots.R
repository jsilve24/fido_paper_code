require(ggplot2)

source("src/plotting.R")

#log_dir <- "results_efficiency/default_optim_params/"
log_dir <- "results_efficiency/run3/"
#log_dir <- ""

datN <- read.csv(paste(log_dir,"run_N-varying.log",sep=""))
datD <- read.csv(paste(log_dir,"run_D-varying.log",sep=""))
datQ <- read.csv(paste(log_dir,"run_Q-varying.log",sep=""))

# plot seconds per effective sample size

print("Plotting SpES...")
plot_SpES(datN, "N", image_filename=paste(log_dir,"N_SpES_plot.png",sep=""))
plot_SpES(datD, "D", image_filename=paste(log_dir,"D_SpES_plot.png",sep=""))
plot_SpES(datQ, "Q", image_filename=paste(log_dir,"Q_SpES_plot.png",sep=""))

print("Plotting MSE...")
plot_accuracy(datN, "N", image_filename=paste(log_dir,"N_MSE_plot.png",sep=""))
plot_accuracy(datD, "D", image_filename=paste(log_dir,"D_MSE_plot.png",sep=""))
plot_accuracy(datQ, "Q", image_filename=paste(log_dir,"Q_MSE_plot.png",sep=""))

print("Plotting 95 CI...")
plot_accuracy(datN, "N", use_95CI=TRUE, image_filename=paste(log_dir,"N_95CI_plot.png",sep=""))
plot_accuracy(datD, "D", use_95CI=TRUE, image_filename=paste(log_dir,"D_95CI_plot.png",sep=""))
plot_accuracy(datQ, "Q", use_95CI=TRUE, image_filename=paste(log_dir,"Q_95CI_plot.png",sep=""))

print("Plotting zeros...")
plot_sparsity(datN, "N", image_filename=paste(log_dir,"N_zeros_plot.png",sep=""))
plot_sparsity(datD, "D", image_filename=paste(log_dir,"D_zeros_plot.png",sep=""))
plot_sparsity(datQ, "Q", image_filename=paste(log_dir,"Q_zeros_plot.png",sep=""))







