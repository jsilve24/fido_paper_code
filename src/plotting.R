require(purrr)
require(driver)
require(dplyr)
require(ggstance)


#' Plot Lambdas from multiple model outputs on same plot
#' 
#' @param mfits named list of mfit objects
#' @param Lambda_true if TRUE extracts Lambda_true from mfits and plots
#' @param image_filename if not NULL save the plot out to a file of this name
#' 
#' @details 
#' Covariates are the horizontal facets, other things are labeled
plot_lambda <- function(mfits, Lambda_true=NULL, image_filename=NULL){
  if (!is.null(Lambda_true)){ 
    lt <- gather_array(Lambda_true, value, coord, covariate)
  }
  
  p <- mfits %>% 
    map("Lambda") %>% 
    map(gather_array, value, coord, covariate, iter) %>%
    bind_rows(.id = "Model") %>% 
    group_by(Model, coord, covariate) %>% 
    summarise_posterior(value) %>% 
    ungroup() %>% 
    ggplot(aes(x=mean, y=coord))
  if (!is.null(Lambda_true)){
     p <- p + geom_segment(data=lt, 
                           aes(x=value, xend=value, 
                               y=coord-.3, yend=coord+.3))
  }
  p <- p +
    geom_pointrangeh(aes(xmin=p2.5, xmax=p97.5, color=Model), 
                   position=position_dodge2v(height=.3)) +
    facet_grid(~covariate) +
    ylab("Coordinate") +
    theme_minimal() +
    theme(axis.title.x=element_blank()) +
    scale_color_brewer(palette="Set1")
  if(!is.null(image_filename)) {
    ggsave(image_filename, plot=last_plot(), width=10, height=10, dpi=600)
  } else {
    show(p)
  }
}




#' Plot Eta from multiple model outputs on same plot
#' 
#' @param mfits named list of mfit objects
#' @param Eta_true if TRUE extracts Eta_true from mfits and plots
#' @param image_filename if not NULL save the plot out to a file of this name
#' 
#' @details 
#' Covariates are the horizontal facets, other things are labeled
plot_eta <- function(mfits, Eta_true=NULL, image_filename=NULL){
  if (!is.null(Eta_true)){ 
    lt <- gather_array(Eta_true, value, coord, covariate)
  }
  
  p <- mfits %>% 
    map("Eta") %>% 
    map(gather_array, value, coord, covariate, iter) %>%
    bind_rows(.id = "Model") %>% 
    group_by(Model, coord, covariate) %>% 
    summarise_posterior(value) %>% 
    ungroup() %>% 
    ggplot(aes(x=mean, y=coord))
  if (!is.null(Eta_true)){
    p <- p + geom_segment(data=lt, 
                          aes(x=value, xend=value, 
                              y=coord-.3, yend=coord+.3))
  }
  p <- p +
    geom_pointrangeh(aes(xmin=p2.5, xmax=p97.5, color=Model), 
                     position=position_dodge2v(height=.3)) +
    facet_grid(~covariate, scales="free") +
    ylab("Coordinate") +
    theme_minimal() +
    theme(axis.title.x=element_blank()) +
    scale_color_brewer(palette="Set1")
  if(!is.null(image_filename)) {
    ggsave(image_filename, plot=last_plot(), width=10, height=10, dpi=600)
  } else {
    show(p)
  }
}


#' Plot runtime as seconds per effective sample (SpES)
#' 
#' @param dat table of results from simulation runs
#' @param varying_param x-axis parameter, increasing in dimension: 
#' "N", "D", or "Q"
#' @param image_filename if not NULL, filename (including extension) 
#' to save plot out to
#' 
#' @details Columns of the data table are expected to be: model, ESS,
#' burnin_runtime, sample_runtime, N, D, Q, total_iter, burnin_iter,
#' percent_zero, lambda_RMSE, percent_outside_95CI, random_seed
plot_SpES <- function(dat, varying_param, image_filename=NULL, log_y=TRUE){
  p <- dat %>% 
    mutate(SpES=1/(ESS/(sample_runtime+burnin_runtime))) %>% 
    ggplot(aes_string(x=varying_param, y="SpES", color="model")) +
    geom_point() + 
    geom_smooth(method="lm") +
    scale_x_log10()
  if (log_y) {
    p <- p + scale_y_log10()
  }
  if(!is.null(image_filename)) {
    ggsave(image_filename, plot=last_plot(), width=10, height=10, dpi=600)
  } else {
    show(p)
  }
}


#' Plot accuracy over increasing dimension
#' 
#' @param dat table of results from simulation runs
#' @param varying_param x-axis parameter, increasing in dimension: 
#' "N", "D", or "Q"
#' @param use_95CI if TRUE, use percent of lambdas outside 95% 
#' CI as accuracy measure
#' @param image_filename if not NULL, filename (including extension) 
#' to save plot out to
#' 
#' @details Columns of the data table are expected to be: model, ESS,
#' burnin_runtime, sample_runtime, N, D, Q, total_iter, burnin_iter,
#' percent_zero, lambda_RMSE, percent_outside_95CI, random_seed
plot_accuracy <- function(dat, varying_param, use_95CI=FALSE, image_filename=NULL, fit_line=TRUE){
  accuracy_measure <- "lambda_RMSE"
  if(use_95CI) {
    accuracy_measure <- "percent_outside_95CI"
  }
  p <- dat %>% 
    ggplot(aes_string(x=varying_param, y=accuracy_measure, color="model")) +
    geom_point()
  if (fit_line) {
    P <- p + geom_smooth(method="lm")
  }
  p <- p + scale_x_log10()
  if(use_95CI) {
    p <- p + ylim(0, 1)
  }
  if(!is.null(image_filename)) {
    ggsave(image_filename, plot=last_plot(), width=10, height=10, dpi=600)
  } else {
    show(p)
  }
}


#' Plot percent zero counts over increasing dimension
#' 
#' @param dat table of results from simulation runs
#' @param varying_param x-axis parameter, increasing in dimension: 
#' "N", "D", or "Q"
#' @param image_filename if not NULL, filename (including extension) 
#' to save plot out to
#' 
#' @details Columns of the data table are expected to be: model, ESS,
#' burnin_runtime, sample_runtime, N, D, Q, total_iter, burnin_iter,
#' percent_zero, lambda_RMSE, percent_outside_95CI, random_seed
plot_sparsity <- function(dat, varying_param, image_filename=NULL){
  p <- dat %>% 
    ggplot(aes_string(x=varying_param, y="percent_zero")) +
    geom_point() +
    scale_x_log10() +
    ylim(0, 1)
  if(!is.null(image_filename)) {
    ggsave(image_filename, plot=last_plot(), width=10, height=10, dpi=600)
  } else {
    show(p)
  }
}





