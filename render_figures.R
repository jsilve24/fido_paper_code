require(ggplot2)
library(dplyr)
library(grid)
library(gridExtra)

dat <- read.csv("all_output.txt")

# just grab zeros from one model
zeros <- filter(dat, model == "Mongrel (Cholesky)")

p1 <- zeros %>% 
  ggplot(aes(x=sweep_value, y=percent_zero)) +
  geom_point() +
  scale_x_log10() +
  ylim(0, 1) +
  theme_minimal() +
  xlab("Sweep parameter") +
  ylab("Percent zero counts")
p1 <- p1 + facet_grid(sweep_param ~ ., scales="free")
p1 <- p1 + theme(panel.spacing = unit(2, "lines")) +
  theme(legend.position="bottom") +
  theme(plot.margin=unit(c(0,0.5,2,0),"cm"))

p2 <- dat %>% 
  mutate(SpES=1/(ESS/(sample_runtime+burnin_runtime))) %>% 
  ggplot(aes(x=sweep_value, y=SpES, color=model)) +
  geom_point() + 
  geom_smooth(method="lm") +
  scale_x_log10() +
  scale_y_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) + 
  theme_minimal() +
  xlab("Sweep parameter") +
  ylab("Seconds per effective sample size")
p2 <- p2 + facet_grid(sweep_param ~ ., scales="free")
p2 <- p2 + theme(panel.spacing = unit(2, "lines")) +
  theme(legend.position="bottom") + 
  theme(legend.title = element_text(size=8)) +
  theme(legend.text = element_text(size=8)) +
  theme(plot.margin=unit(c(0,0.5,0,0.5),"cm"))

p3 <- dat %>% 
  ggplot(aes(x=sweep_value, y=lambda_MSE, color=model)) +
  geom_point() +
  #geom_smooth(method="lm") +
  scale_x_log10() +
  theme_minimal() +
  xlab("Sweep parameter") +
  ylab("MSE elements of Lambda")
p3 <- p3 + facet_grid(sweep_param ~ ., scales="free")
p3 <- p3 + theme(panel.spacing = unit(2, "lines")) +
  theme(legend.position="bottom") +
  theme(legend.title = element_text(size=8)) +
  theme(legend.text = element_text(size=8)) +
  theme(plot.margin=unit(c(0,0.5,0,0.5),"cm"))

grid.arrange(p1, p2, p3, nrow=1)


