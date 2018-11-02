require(ggplot2)

source("src/plotting.R")

log_dir <- "results_efficiency/run4/"
#log_dir <- ""

#log_suffix <- "_2018-10-25"
log_suffix <- "-output"
exclude_su_mcp <- FALSE
file_suffix <- ""

datN <- read.csv(paste(log_dir,"N",log_suffix,".log",sep=""))
datD <- read.csv(paste(log_dir,"D",log_suffix,".log",sep=""))
datQ <- read.csv(paste(log_dir,"Q",log_suffix,".log",sep=""))

if(exclude_su_mcp) {
  datN <- datN[datN$model != "stan_uncollapsed",]
  datD <- datD[datD$model != "stan_uncollapsed",]
  datQ <- datQ[datQ$model != "stan_uncollapsed",]
} else {
  file_suffix <- "_all"
}

# take out the trash
datN <- datN[datN$model != "mongrel_cholesky_partial",]
datD <- datD[datD$model != "mongrel_cholesky_partial",]
datQ <- datQ[datQ$model != "mongrel_cholesky_partial",]

# plot seconds per effective sample size

#print("Plotting SpES...")
#plot_SpES(datN, "N", image_filename=paste(log_dir,"SpES_N_plot_log",file_suffix,".png",sep=""))
#plot_SpES(datD, "D", image_filename=paste(log_dir,"SpES_D_plot_log",file_suffix,".png",sep=""))
#plot_SpES(datQ, "Q", image_filename=paste(log_dir,"SpES_Q_plot_log",file_suffix,".png",sep=""))

#plot_SpES(datN, "N", image_filename=paste(log_dir,"SpES_N_plot",file_suffix,".png",sep=""), log_y=FALSE)
#plot_SpES(datD, "D", image_filename=paste(log_dir,"SpES_D_plot",file_suffix,".png",sep=""), log_y=FALSE)
#plot_SpES(datQ, "Q", image_filename=paste(log_dir,"SpES_Q_plot",file_suffix,".png",sep=""), log_y=FALSE)

print("Plotting MSE...")
plot_accuracy(datN, "N", image_filename=paste(log_dir,"MSE_N_plot",file_suffix,".png",sep=""), fit_line=FALSE)
plot_accuracy(datD, "D", image_filename=paste(log_dir,"MSE_D_plot",file_suffix,".png",sep=""), fit_line=FALSE)
plot_accuracy(datQ, "Q", image_filename=paste(log_dir,"MSE_Q_plot",file_suffix,".png",sep=""), fit_line=FALSE)

#print("Plotting 95 CI...")
#plot_accuracy(datN, "N", use_95CI=TRUE, image_filename=paste(log_dir,"95CI_N_plot",file_suffix,".png",sep=""))
#plot_accuracy(datD, "D", use_95CI=TRUE, image_filename=paste(log_dir,"95CI_D_plot",file_suffix,".png",sep=""))
#plot_accuracy(datQ, "Q", use_95CI=TRUE, image_filename=paste(log_dir,"95CI_Q_plot",file_suffix,".png",sep=""))

#print("Plotting zeros...")
#plot_sparsity(datN, "N", image_filename=paste(log_dir,"zeros_N_plot.png",sep=""))
#plot_sparsity(datD, "D", image_filename=paste(log_dir,"zeros_D_plot.png",sep=""))
#plot_sparsity(datQ, "Q", image_filename=paste(log_dir,"zeros_Q_plot.png",sep=""))







