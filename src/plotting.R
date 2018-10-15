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
  p +
    geom_pointrangeh(aes(xmin=p2.5, xmax=p97.5, color=Model), 
                   position=position_dodge2v(height=.3)) +
    facet_grid(~covariate) +
    ylab("Coordinate") +
    theme_minimal() +
    theme(axis.title.x=element_blank()) +
    scale_color_brewer(palette="Set1")
  if(!is.null(image_filename)) {
    ggsave(image_filename, plot=last_plot(), width=10, height=10, dpi=1200)
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
  p +
    geom_pointrangeh(aes(xmin=p2.5, xmax=p97.5, color=Model), 
                     position=position_dodge2v(height=.3)) +
    facet_grid(~covariate, scales="free") +
    ylab("Coordinate") +
    theme_minimal() +
    theme(axis.title.x=element_blank()) +
    scale_color_brewer(palette="Set1")
  if(!is.null(image_filename)) {
    ggsave(image_filename, plot=last_plot(), width=10, height=10, dpi=1200)
  }
}
