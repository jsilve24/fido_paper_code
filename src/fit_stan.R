require(rstan)
require(mongrel)
source("src/fit_methods.R")
source("src/dataset_methods.R")
source("src/utils.R")


# Options for stan - recommended
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)



# main multinomial function ---------------------------------------------------

#' Fit multinomial model using stan
#' 
#' WARNING: Req that your current working directory is main folder
#'  "mongrel_paper_code"
#'  
#' WARNING: Currently uses a random seed... 
#' 
#' @param mdataset an mdataset object
#' @param chains number of chains to run
#' @param iter number of samples from each chain (note: includes warmup)
#' @param ... other parameters passed to function stan
#' @return mfit object
fit_mstan <- function(mdataset, chains=4, iter=2000, ...){
  init <- list()
  for (i in 1:chains){
    init[[i]] <- list(eta=mongrel::random_mongrel_init(mdataset$Y))
  }
  
  fit <- stan("src/rel_lm_collapsed.stan", data=mdataset, chains=chains, 
              init=init, iter=iter, pars=c("B", "Sigma"),  ...)
  pars <- rstan::extract(fit, c("B", "Sigma"))
  rm(fit) # free up memory
  
  m <- mfit(N=mdataset$N, D=mdataset$D, Q=mdataset$Q, iter=dim(pars$B)[1], 
            Lambda=aperm(pars$B, c(2,3,1)), 
            Sigma=aperm(pars$Sigma, c(2,3,1)), 
            mdataset=mdataset)
  
  return(m)
}


