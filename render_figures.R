require(ggplot2)
library(dplyr)
library(grid)
library(gridExtra)

use_legend <- TRUE

figure_no <- 1

apply_common_settings <- function(p, use_legend=FALSE, strip_ylabels=FALSE, horiz=FALSE) {
  p <- p + theme_minimal() +
    theme(panel.border = element_rect(color="black", fill=NA, size=0.5)) + 
    xlab("") +
    ylab("")
  if(horiz) {
    p <- p + facet_grid(. ~ sweep_param, scales="free")
  } else {
    p <- p + facet_grid(sweep_param ~ ., scales="free")
  }
  p <- p + theme(panel.spacing = unit(2, "lines")) +
    theme(aspect.ratio = 1) +
    theme(plot.margin=unit(c(0,0,0,0),"pt"))
  if(strip_ylabels) {
    p <- p + theme(strip.text.y = element_blank())
  } else {
    if(horiz) {
      p <- p + theme(strip.text.x = element_text(angle=0, face="bold", size=14))
    } else {
      p <- p + theme(strip.text.y = element_text(angle=0, face="bold", size=14))
    }
  }
  if(use_legend) {
    p <- p + theme(legend.position="bottom") +
      theme(legend.position="bottom",legend.direction="vertical") + 
      theme(legend.title = element_blank()) +
      theme(legend.text = element_text(size=7))
  } else {
    p <- p + theme(legend.position="none")
  }
  return(p)
}

if(figure_no == 1) {
  dat <- read.csv("all_output.txt")

  # just grab zeros from one model
  zeros <- filter(dat, model == "Mongrel (Cholesky)")

  p1 <- zeros %>% 
    ggplot(aes(x=sweep_value, y=percent_zero)) +
    geom_point() +
    scale_x_log10() +
    ylim(0, 1)
  p1 <- apply_common_settings(p1, use_legend=use_legend, strip_ylabels=TRUE)

  p2 <- dat %>% 
    mutate(SpES=1/(ESS/(sample_runtime+burnin_runtime))) %>% 
    ggplot(aes(x=sweep_value, y=SpES, color=model)) +
    geom_point() + 
    geom_smooth(method="lm") +
    scale_x_log10() +
    scale_y_log10(
      breaks = scales::trans_breaks("log10", function(x) 10^x),
      labels = scales::trans_format("log10", scales::math_format(10^.x))
    )
  p2 <- apply_common_settings(p2, use_legend=use_legend, strip_ylabels=TRUE)

  p3 <- dat %>% 
    ggplot(aes(x=sweep_value, y=lambda_MSE, color=model)) +
    geom_point() +
    scale_x_log10()
  p3 <- apply_common_settings(p3, use_legend=use_legend)

  p <- grid.arrange(p1, p2, p3, nrow=1)
  p
#  ggsave("figure_drafts/F1.png", plot=p, dpi=300, width=11, height=10, units="in")
}

if(figure_no == 2) {
  data <- read.csv("second_moment_data.log")

  p <- ggplot(data, aes(x=sweep_value, y=sd_MSE, color=model)) +
    geom_point() +
    scale_x_log10()
  p <- apply_common_settings(p, use_legend=FALSE, strip_ylabels=FALSE, horiz=TRUE)
  p
  ggsave("figure_drafts/F2.png", plot=p, dpi=300, width=11, height=4, units="in")

  # combine replicates?
}