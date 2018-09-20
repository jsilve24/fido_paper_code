require(purrr)
require(driver)
require(dplyr)
require(ggstance)


#' Plot Lambdas from multiple model outputs on same plot
#' 
#' @param mfits named list of mfit objects
#' 
#' @details 
#' Covariates are the horizontal facets, other things are labeled
plot_lambda <- function(mfits){
  mfits %>% 
    map("Lambda") %>% 
    map(gather_array, value, coord, covariate, iter) %>%
    bind_rows(.id = "Model") %>% 
    group_by(Model, coord, covariate) %>% 
    summarise_posterior(value) %>% 
    ungroup() %>% 
    ggplot(aes(x=mean, y=coord, xmin=p2.5, xmax=p97.5, color=Model))+
    geom_pointrangeh(position=position_dodge2v(height=.3)) +
    facet_grid(~covariate) +
    ylab("Coordinate") +
    theme_minimal() +
    theme(axis.title.x=element_blank()) +
    scale_color_brewer(palette="Set1")
}
