require(ggplot2)

source("src/plotting.R")

#log_dir <- "results_efficiency/default_optim_params/"
log_dir <- ""

datN <- read.csv(paste(log_dir,"run_N-varying.log",sep=""))
datD <- read.csv(paste(log_dir,"run_D-varying.log",sep=""))
datQ <- read.csv(paste(log_dir,"run_Q-varying.log",sep=""))

plot_SpES(datN, "N", image_filename=paste(log_dir,"N_SpES_plot.png",sep=""))
plot_SpES(datD, "D", image_filename=paste(log_dir,"D_SpES_plot.png",sep=""))
plot_SpES(datQ, "Q", image_filename=paste(log_dir,"Q_SpES_plot.png",sep=""))
