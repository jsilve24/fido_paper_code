library(dplyr)
# library(grid)
# library(gridExtra)
library(driver)
library(ggplot2)
library(purrr)
library(ggstance)
library(ggpubr)
library(RColorBrewer)

build_color_palette <- function() {
  my_colors <- brewer.pal(8, "Set1")
  names(my_colors) <- c("mongrel_cholesky",
                        "stan_collapsed",
                        "stan_uncollapsed",
                        "conjugate_linear_model",
                        "stan_collapsed_variationalbayes_meanfield",
                        "stan_collapsed_variationalbayes_fullrank",
                        "stan_uncollapsed_variationalbayes_meanfield",
                        "stan_uncollapsed_variationalbayes_fullrank")
  return(scale_colour_manual(name="model", values=my_colors))
}

apply_common_settings <- function(base_p, horiz=FALSE, use_legend=FALSE, share_y=FALSE, use_manual_color=TRUE) {
  # apply grid
  if(horiz) {
    p <- base_p + facet_grid(. ~ sweep_param, scales="free_x")
  } else {
    if(share_y) {
      p <- base_p + facet_grid(sweep_param ~ .)
    } else {
      p <- base_p + facet_grid(sweep_param ~ ., scales="free_y")
    }
  }

  # spacing
  p <- p + theme(panel.spacing = unit(0.5, "lines")) +
    theme(plot.margin=unit(c(0,0,0,0),"pt")) +
    theme(aspect.ratio = 1)

  # apply theme and border
  p <- p + theme_minimal() +
    theme(panel.border = element_rect(color="black", fill=NA, size=0.5))

  if(!use_legend) {
    p <- p + theme(legend.position="none")
  }

  # remove legend and labels
  p <- p + xlab("")
  if(!use_legend) {
    p <- p + theme(strip.text.x = element_blank()) +
        theme(strip.text.y = element_blank())
  }

  if(use_manual_color) {
    col_scale <- build_color_palette()
    p <- p + col_scale
  } else {
    p <- p + scale_colour_brewer(palette = "Set1")
  }

  # later, maybe?
  # theme(aspect.ratio = 1)

  return(p)
}

render_F1_C1 <- function(dat, use_legend=FALSE) {
  # just grab zeros from one model, they're (verifiably) all the same
  zeros <- filter(dat, model == "mongrel_cholesky")
  p <- zeros %>% 
    ggplot(aes(x=sweep_value, y=percent_zero)) +
    geom_point() +
    geom_smooth(method="loess", se=FALSE, color="black") +
    scale_x_log10() +
    ylim(0, 1) +
    ylab("percent zeros")
  p <- apply_common_settings(p, use_legend=use_legend)
  return(p)
}

render_F1_C2 <- function(dat, use_CLM=FALSE, use_legend=FALSE) {
  dat_filtered <- filter(dat, !(model %in% c("mongrel_eigen")))
  if(!use_CLM) {
    dat_filtered <- filter(dat_filtered, !(model %in% c("conjugate_linear_model")))
  }
  p <- dat_filtered %>% 
    mutate(SpES=1/(ESS/(sample_runtime+burnin_runtime))) %>% 
    ggplot(aes(x=sweep_value, y=SpES, color=model)) +
    geom_point() + 
#    geom_smooth(method="loess", se=FALSE) +
    scale_x_log10() +
    scale_y_log10(
      breaks = scales::trans_breaks("log10", function(x) 10^x),
      labels = scales::trans_format("log10", scales::math_format(10^.x))
    ) +
    ylab("second per effective sample size")
  p <- apply_common_settings(p, use_legend=use_legend)
  return(p)
}

render_F1_C3 <- function(dat, use_CLM=FALSE, use_legend=FALSE) {
  dat_filtered <- filter(dat, !(model %in% c("mongrel_eigen")))
  if(!use_CLM) {
    dat_filtered <- filter(dat_filtered, !(model %in% c("conjugate_linear_model")))
  }
  p <- dat_filtered %>% 
    ggplot(aes(x=sweep_value, y=lambda_RMSE, color=model)) +
    geom_point() +
#    geom_smooth(method="loess", se=FALSE) +
    scale_x_log10() +
    ylab(expression(paste("MSE ",Lambda,sep="")[ij]))
  p <- apply_common_settings(p, use_legend=use_legend)
  return(p)
}

render_F1_C4 <- function(use_CLM=FALSE, use_legend=FALSE, Q_only=FALSE) {
  dat <- read.csv("second_moment_data.log")
  if(!use_CLM) {
    dat <- filter(dat, !(model %in% c("conjugate_linear_model")))
  }
  if(Q_only) {
    dat <- filter(dat, sweep_param=="Q")
  }
  
  p <- ggplot(dat, aes(x=sweep_value, y=sd_RMSE, color=model)) +
    geom_point() +
#    geom_smooth(method="loess", se=FALSE) +
    scale_x_log10() +
    ylab("MSE of deviation")

  p <- apply_common_settings(p, use_legend=use_legend, share_y=TRUE)
  return(p)
}

render_F1 <- function(use_legend=FALSE) {
  cat("Rendering Figure 1\n")
  dat <- read.csv("first_moment_data.log")
  c1 <- render_F1_C1(dat, use_legend=FALSE)

  # render without CLM
  c2 <- render_F1_C2(dat, use_legend=FALSE)
  c3 <- render_F1_C3(dat, use_legend=FALSE)
  c4 <- render_F1_C4(use_legend=use_legend)
  if(use_legend) {
    p <- ggarrange(c1, c2, c3, c4, ncol=4, nrow=1, widths = c(1, 1, 1, 2))
    ggsave("figure_drafts/legends/F1.png", plot=p, width=16, height=8, units="in")  
  } else {
    p <- ggarrange(c1, c2, c3, c4, ncol=4, nrow=1, widths = c(1, 1, 1, 1))
    ggsave("figure_drafts/no_legends/F1.png", plot=p, width=10, height=6, units="in")  
  }
}

render_S1 <- function(use_legend=FALSE) {
  cat("Rendering Supplemental Figure 1\n")
  p <- render_F1_C4(use_CLM=TRUE, use_legend=use_legend, Q_only=TRUE)
  if(use_legend) {
    ggsave("figure_drafts/legends/S1.png", plot=p, width=6.5, height=3, units="in")
  } else {
    ggsave("figure_drafts/no_legends/S1.png", plot=p, width=4, height=4, units="in")
  }
}

subset_Lambda <- function(mfits, lower, upper) {
  Lambda_true <- mfits$mongrel_cholesky$mdataset$Lambda_true[ll:ul,ll:ul]
  mfits$mongrel_cholesky$Lambda <- mfits$mongrel_cholesky$Lambda[ll:ul,ll:ul,]
  mfits$stan_collapsed$Lambda <- mfits$stan_collapsed$Lambda[ll:ul,ll:ul,]
  mfits$stan_uncollapsed$Lambda <- mfits$stan_uncollapsed$Lambda[ll:ul,ll:ul,]
  mfits$conjugate_linear_model$Lambda <- mfits$conjugate_linear_model$Lambda[ll:ul,ll:ul,]
  return(list("Lambda_true"=Lambda_true, "mfits"=mfits))
}

plot_posterior_intervals <- function(mfits, ll, ul, use_legend=FALSE, image_filename=NULL) {
  Lambda_true <- mfits$mongrel_cholesky$mdataset$Lambda_true[ll:ul,ll:ul]
  mfits$mongrel_cholesky$Lambda <- mfits$mongrel_cholesky$Lambda[ll:ul,ll:ul,]
  mfits$stan_collapsed$Lambda <- mfits$stan_collapsed$Lambda[ll:ul,ll:ul,]
  mfits$stan_uncollapsed$Lambda <- mfits$stan_uncollapsed$Lambda[ll:ul,ll:ul,]
  mfits$conjugate_linear_model$Lambda <- mfits$conjugate_linear_model$Lambda[ll:ul,ll:ul,]
  mfits$stan_collapsed_variationalbayes_meanfield$Lambda <- mfits$stan_collapsed_variationalbayes_meanfield$Lambda[ll:ul,ll:ul,]

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
    theme(panel.border = element_rect(color="black", fill=NA, size=0.5))
  if(!use_legend) {
    p <- p + theme(legend.position="none")
  }
  col_scale <- build_color_palette()
  p <- p + col_scale
  if(!is.null(image_filename)) {
    if(use_legend) {
      ggsave(paste("figure_drafts/legends/",image_filename,sep=""), plot=p, width=10, height=8, dpi=300)
    } else {
      ggsave(paste("figure_drafts/no_legends/",image_filename,sep=""), plot=p, width=8, height=8, dpi=300)
    }
  }
}

render_S2 <- function(use_legend=FALSE) {
  cat("Rendering Supplemental Figure 2a\n")
  load("fitted_models/MC_N30_D30_Q5_R1.RData")
  load("fitted_models/SC_N30_D30_Q5_R1.RData")
  load("fitted_models/SU_N30_D30_Q5_R1.RData")
  load("fitted_models/CLM_N30_D30_Q5_R1.RData")
  load("fitted_models/SVBCM_N30_D30_Q5_R1.RData")
  mfits <- list("mongrel_cholesky"=fit.mc,
                "stan_collapsed"=fit.sc,
                "stan_uncollapsed"=fit.su,
                "conjugate_linear_model"=fit.clm,
                "stan_collapsed_variationalbayes_meanfield"=fit.svbcm)
  plot_posterior_intervals(mfits, 1, 5, use_legend=use_legend, image_filename="S2_good_case.png")

  cat("Rendering Supplemental Figure 2b\n")
  load("fitted_models/MC_N100_D30_Q250_R1.RData")
  load("fitted_models/SC_N100_D30_Q250_R1.RData")
  load("fitted_models/SU_N100_D30_Q250_R1.RData")
  load("fitted_models/CLM_N100_D30_Q250_R1.RData")
  load("fitted_models/SVBCM_N100_D30_Q250_R1.RData")
  mfits <- list("mongrel_cholesky"=fit.mc,
               "stan_collapsed"=fit.sc,
               "stan_uncollapsed"=fit.su,
               "conjugate_linear_model"=fit.clm,
               "stan_collapsed_variationalbayes_meanfield"=fit.svbcm)
  plot_posterior_intervals(mfits, 1, 5, use_legend=use_legend, image_filename="S2_fail_case.png")
}

extract_opt_cols <- function(df, old_colname, new_colname) {
  d <- select(df, c("sweep_param",
                      "sweep_value",
                      old_colname)) %>%
    mutate(opt_part=new_colname)
  colnames(d)[which(colnames(d)==old_colname)] <- "value"
  return(d)
}

# proportion of runtime taken up by optimization 
render_S3 <- function(use_legend=FALSE) {
  cat("Rendering Supplemental Figure 3\n")
  dat <- read.csv("first_moment_data.log")
  dat <- filter(dat, model %in% c("mongrel_cholesky"))
  dat <- select(dat, c("sweep_param",
                       "sweep_value",
                       "Timer_HessianCalculation",
                       "Timer_LaplaceApproximation",
                       "Timer_Optimization",
                       "Timer_Overall",
                       "Timer_Uncollapse_Overall.Overall"))
  d1 <- extract_opt_cols(dat, "Timer_HessianCalculation", "HessianCalculation")
  d2 <- extract_opt_cols(dat, "Timer_LaplaceApproximation", "LaplaceApproximation")
  d3 <- extract_opt_cols(dat, "Timer_Optimization", "Optimization")
  d4 <- extract_opt_cols(dat, "Timer_Overall", "Overall")
  d5 <- extract_opt_cols(dat, "Timer_Uncollapse_Overall.Overall", "Uncollapse_Overall.Overall")
  dat <- rbind(d1, d2, d3, d4, d5)
  p <- dat %>%
    ggplot(aes(x=sweep_value, y=value, color=opt_part)) +
    geom_point() + 
    geom_smooth(method="loess", se=FALSE) + 
    scale_x_log10() +
    scale_y_log10(
      breaks = scales::trans_breaks("log10", function(x) 10^x),
      labels = scales::trans_format("log10", scales::math_format(10^.x))
    )
#    facet_grid(. ~ sweep_param, scales="free_x")
  p <- apply_common_settings(p, horiz=TRUE, use_legend=use_legend, use_manual_color=FALSE)
  if(use_legend) {
    ggsave("figure_drafts/legends/S3.png", plot=p, width=8, height=2.5, dpi=300)
  } else {
    ggsave("figure_drafts/no_legends/S3.png", plot=p, width=8, height=3, dpi=300)
  }
}

use_l <- TRUE
render_F1(use_legend=use_l)
render_S1(use_legend=use_l)
render_S2(use_legend=use_l)
render_S3(use_legend=use_l)

use_l <- FALSE
render_F1(use_legend=use_l)
render_S1(use_legend=use_l)
render_S2(use_legend=use_l)
render_S3(use_legend=use_l)

