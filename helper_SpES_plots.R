require(ggplot2)

source("src/plotting.R")

datN <- read.csv("results_efficiency/default_optim_params/run_N-varying.log")
datD <- read.csv("results_efficiency/default_optim_params/run_D-varying.log")
datQ <- read.csv("results_efficiency/default_optim_params/run_Q-varying.log")

plot_SpES(datN, "N", image_filename="results_efficiency/default_optim_params/N_SpES_plot.png")
plot_SpES(datD, "D", image_filename="results_efficiency/default_optim_params/D_SpES_plot.png")
plot_SpES(datQ, "Q", image_filename="results_efficiency/default_optim_params/Q_SpES_plot.png")
