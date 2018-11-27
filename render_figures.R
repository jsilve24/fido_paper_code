library(dplyr)
library(grid)
library(gridExtra)
library(driver)
library(ggplot2)
library(purrr)
library(ggstance)

apply_common_settings <- function(base_p, horiz=FALSE, use_legend=FALSE) {
  # apply grid
  if(horiz) {
    p <- base_p + facet_grid(. ~ sweep_param, scales="free_x")
  } else {
    p <- base_p + facet_grid(sweep_param ~ ., scales="free_y")
  }

  # spacing
  p <- p + theme(panel.spacing = unit(0.5, "lines")) +
    theme(aspect.ratio = 1) +
    theme(plot.margin=unit(c(0,0,0,0),"pt"))

  # apply theme and border
  p <- p + theme_minimal() +
    theme(panel.border = element_rect(color="black", fill=NA, size=0.5))

  # remove legend and labels
  p <- p + theme(legend.position="none") +
    xlab("") +
    ylab("") +
    theme(strip.text.x = element_blank()) +
    theme(strip.text.y = element_blank())

  return(p)
}

render_F1_C1 <- function(dat) {
  # just grab zeros from one model, they're (verifiably) all the same
  zeros <- filter(dat, model == "mongrel_cholesky")
  p <- zeros %>% 
    ggplot(aes(x=sweep_value, y=percent_zero)) +
    geom_point() +
    scale_x_log10() +
    ylim(0, 1)
  p <- apply_common_settings(p)
  return(p)
}

render_F1_C2 <- function(dat, use_CLM=FALSE) {
  # dat includes: SU, SC, ME, MC, CLM
  dat_filtered <- filter(dat, !(model %in% c("mongrel_eigen")))
  if(!use_CLM) {
    dat_filtered <- filter(dat_filtered, !(model %in% c("conjugate_linear_model")))
  }
  p <- dat_filtered %>% 
    mutate(SpES=1/(ESS/(sample_runtime+burnin_runtime))) %>% 
    ggplot(aes(x=sweep_value, y=SpES, color=model)) +
    geom_point() + 
    geom_smooth(method="lm") +
    scale_x_log10() +
    scale_y_log10(
      breaks = scales::trans_breaks("log10", function(x) 10^x),
      labels = scales::trans_format("log10", scales::math_format(10^.x))
    )
  p <- apply_common_settings(p)
  return(p)
}

render_F1_C3 <- function(dat, use_CLM=FALSE) {
  dat_filtered <- filter(dat, !(model %in% c("mongrel_eigen")))
  if(!use_CLM) {
    dat_filtered <- filter(dat_filtered, !(model %in% c("conjugate_linear_model")))
  }
  p <- dat_filtered %>% 
    ggplot(aes(x=sweep_value, y=lambda_MSE, color=model)) +
    geom_point() +
    scale_x_log10()
  p <- apply_common_settings(p)
  return(p)
}

render_F1 <- function() {
  dat <- read.csv("all_output.log")
  c1 <- render_F1_C1(dat)

  # render once including conjugate linear model
  c2 <- render_F1_C2(dat, use_CLM=TRUE)
  c3 <- render_F1_C3(dat, use_CLM=TRUE)
  p <- grid.arrange(c1, c2, c3, nrow=1)
  ggsave("figure_drafts/F1_wCLM.png", plot=p, width=10, height=8, units="in")  

  # render once without
  c2 <- render_F1_C2(dat)
  c3 <- render_F1_C3(dat)
  p <- grid.arrange(c1, c2, c3, nrow=1)
  ggsave("figure_drafts/F1.png", plot=p, width=10, height=8, units="in")  
}

render_F2 <- function(use_stan_collapsed=TRUE) {
  # fix log name
  if(use_stan_collapsed) {
    dat <- read.csv("second_moment_data_SCbaseline.log")
  } else {
    dat <- read.csv("second_moment_data_SUbaseline.log")
  }

  # render once including the conjugate linear model
  p <- ggplot(dat, aes(x=sweep_value, y=sd_MSE, color=model)) +
    geom_point() +
    scale_x_log10()
  p <- apply_common_settings(p, use_legend=FALSE, horiz=TRUE)
  if(use_stan_collapsed) {
    ggsave("figure_drafts/F2_CLM_SCbaseline.png", plot=p, width=10, height=4, units="in")
  } else {
    ggsave("figure_drafts/F2_CLM_SUbaseline.png", plot=p, width=10, height=4, units="in")
  }
  
  # render once without CLM
  dat_filtered <- filter(dat, !(model %in% c("conjugate_linear_model")))
  p <- ggplot(dat_filtered, aes(x=sweep_value, y=sd_MSE, color=model)) +
    geom_point() +
    scale_x_log10()
  p <- apply_common_settings(p, use_legend=FALSE, horiz=TRUE)
  if(use_stan_collapsed) {
    ggsave("figure_drafts/F2_SCbaseline.png", plot=p, width=10, height=4, units="in")
  } else {
    ggsave("figure_drafts/F2_SUbaseline.png", plot=p, width=10, height=4, units="in")
  }
  rm(dat)
  rm(dat_filtered)
}

subset_Lambda <- function(mfits, lower, upper) {
  Lambda_true <- mfits$mongrel_cholesky$mdataset$Lambda_true[ll:ul,ll:ul]
  mfits$mongrel_cholesky$Lambda <- mfits$mongrel_cholesky$Lambda[ll:ul,ll:ul,]
  mfits$stan_collapsed$Lambda <- mfits$stan_collapsed$Lambda[ll:ul,ll:ul,]
  mfits$stan_uncollapsed$Lambda <- mfits$stan_uncollapsed$Lambda[ll:ul,ll:ul,]
  mfits$conjugate_linear_model$Lambda <- mfits$conjugate_linear_model$Lambda[ll:ul,ll:ul,]
  return(list("Lambda_true"=Lambda_true, "mfits"=mfits))
}

plot_posterior_intervals <- function(mfits, ll, ul, image_filename=NULL) {
  Lambda_true <- mfits$mongrel_cholesky$mdataset$Lambda_true[ll:ul,ll:ul]
  mfits$mongrel_cholesky$Lambda <- mfits$mongrel_cholesky$Lambda[ll:ul,ll:ul,]
  mfits$stan_collapsed$Lambda <- mfits$stan_collapsed$Lambda[ll:ul,ll:ul,]
  mfits$stan_uncollapsed$Lambda <- mfits$stan_uncollapsed$Lambda[ll:ul,ll:ul,]
  mfits$conjugate_linear_model$Lambda <- mfits$conjugate_linear_model$Lambda[ll:ul,ll:ul,]

  lt <- gather_array(Lambda_true, value, coord, covariate)
  p <- mfits %>% 
    map("Lambda") %>% 
    map(gather_array, value, coord, covariate, iter) %>%
    bind_rows(.id = "Model") %>% 
    group_by(Model, coord, covariate) %>% 
    summarise_posterior(value) %>% 
    ungroup() %>% 
    ggplot(aes(x=mean, y=coord)) +
    geom_segment(data=lt, aes(x=value, xend=value, y=coord-.3, yend=coord+.3)) +
    geom_pointrangeh(aes(xmin=p2.5, xmax=p97.5, color=Model), 
                     position=position_dodge2v(height=.3)) +
    facet_grid(coord ~ covariate, scales="free_y") +
    ylab("Coordinate") +
    theme_minimal() +
    theme(axis.title.x=element_blank()) +
    theme(axis.title.y=element_blank()) +
    theme(axis.text.y=element_blank()) +
    theme(strip.text.x=element_blank()) + 
    theme(strip.text.y=element_blank()) + 
    theme(panel.spacing=unit(0.5, "lines")) +
    theme(legend.position="none") +
    theme(panel.border = element_rect(color="black", fill=NA, size=0.5))
  if(!is.null(image_filename)) {
    ggsave(image_filename, plot=p, width=8, height=8, dpi=300)
  }
}

render_SF1 <- function() {
 load("fitted_models/MC_N100_D30_Q250_R1.RData")
 load("fitted_models/SC_N100_D30_Q250_R1.RData")
 load("fitted_models/SU_N100_D30_Q250_R1.RData")
 load("fitted_models/CLM_N100_D30_Q250_R1.RData")
 mfits <- list("mongrel_cholesky"=fit.mc,
               "stan_collapsed"=fit.sc,
               "stan_uncollapsed"=fit.su,
               "conjugate_linear_model"=fit.clm)
 plot_posterior_intervals(mfits, 1, 5, image_filename="figure_drafts/SF1_fail_case.png")

  load("fitted_models/MC_N30_D30_Q5_R1.RData")
  load("fitted_models/SC_N30_D30_Q5_R1.RData")
  load("fitted_models/SU_N30_D30_Q5_R1.RData")
  load("fitted_models/CLM_N30_D30_Q5_R1.RData")
  mfits <- list("mongrel_cholesky"=fit.mc,
                "stan_collapsed"=fit.sc,
                "stan_uncollapsed"=fit.su,
                "conjugate_linear_model"=fit.clm)
  plot_posterior_intervals(mfits, 1, 5, image_filename="figure_drafts/SF1_good_case.png")
}

render_F1()
render_F2()
render_F2(use_stan_collapsed=FALSE)
render_SF1()

