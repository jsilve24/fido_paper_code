library(stray)
library(tidyverse)
library(phyloseq)
library(driver)
library(Biostrings)
set.seed(59292)



# helper functions --------------------------------------------------------

SE_from_dist <- function(dist, sigma = 1, rho = median(as.matrix(dist(t(X)))), 
                         jitter = 1e-10){
  dist <- as.matrix(dist)
  G <- sigma^2 * exp(-dist^2/(2 * rho^2)) + jitter * diag(ncol(dist))
  return(G)
}

# main analysis -----------------------------------------------------------



data(mallard_family)
setwd("~/Research/mongrel/stray_paper_code/output/")


# Just take hourly samples
mallard_family <- prune_samples((sample_data(mallard_family)$time > "2015-11-20 15:00:00 UTC") &
                                  (sample_data(mallard_family)$time < "2015-11-25 16:00:00 UTC"), mallard_family)

# Order samples - to make plotting easy later
o <- order(sample_data(mallard_family)$time)
otu_table(mallard_family) <- otu_table(mallard_family)[o,]
sample_data(mallard_family) <- sample_data(mallard_family)[o,]

# Extract Data / dimensions from Phyloseq object
Y <- t(as(otu_table(mallard_family), "matrix"))
rownames(Y) <- taxa_names(mallard_family)
D <- ntaxa(mallard_family)
N <- nrow(sample_data(mallard_family))

# X in hours
X <- as.numeric(sample_data(mallard_family)$time)
X <- t((X-min(X)) / 3600)

# extract vessel numbers
X <- rbind(X, sample_data(mallard_family)$Vessel)

# Specify Priors
Gamma.time <- function(X) SE(X[1,,drop=F], sigma=5, rho=10) # Create partial function 
Gamma.vessel <- function(X) {
    Roh <- driver::onehot(as.factor(X[2,]))
    return(tcrossprod(Roh))
}
Gamma <- function(X){
  Gt <- Gamma.time(X)
  Gv <- Gamma.vessel(X)
  return(Gt*Gv)
}
Theta <- function(X) matrix(0, D-1, ncol(X))
upsilon <- D-1+3
stdis <- stringDist(refseq(mallard_family))
Xi <- SE_from_dist(stdis)
GG <- cbind(diag(nrow(Xi)-1), -1)
Xi <- GG %*% Xi %*% t(GG)
Xi <- cov2cor(Xi)

# Now fit the model
fit <- stray::basset(Y, X, upsilon, Theta, Gamma, Xi)
print(fit$Timer)

# Convert to CLR
fit.clr <- to_clr(fit)

# Smooth over observed time points
X_predict <- 1:(max(X))
TT <- length(X_predict)
X_predict <- rep(t(1:(max(X))), 4)
X_predict <- rbind(X_predict, rep(1:4, each=TT))
predicted <- predict(fit.clr, X_predict, jitter=1) 

# Now visualize
family_names <- as(tax_table(mallard_family)[,"Family"], "vector")
Y_clr_tidy <- clr_array(Y+0.65, parts = 1) %>% 
  gather_array(mean, coord, sample) %>% 
  mutate(time = X[1,sample], 
         vessel=X[2,sample],
         coord = paste0("clr_", family_names[coord]))

predicted_tidy <- gather_array(predicted, val, coord, sample, iter) %>% 
  mutate(time = X_predict[1,sample], 
         vessel = X_predict[2,sample]) %>% 
  filter(!is.na(val)) %>% 
  group_by(time, vessel, coord) %>% 
  summarise_posterior(val, na.rm=TRUE) %>% 
  ungroup() %>% 
  mutate(coord = paste0("clr_", family_names[coord]))

p <- ggplot(predicted_tidy, aes(x = time, y=mean, group=vessel)) +
  #geom_point(data = Y_clr_tidy, alpha=0.5, aes(color=factor(vessel))) +
  geom_ribbon(aes(ymin=p2.5, ymax=p97.5), fill="darkgrey", alpha=0.5) +
  geom_ribbon(aes(ymin=p25, ymax=p75), fill="darkgrey", alpha=0.9)+
  geom_line(color="blue") +
  facet_wrap(~coord, scales="free_y", nrow=2) +
  theme_minimal()+
  theme(axis.title.y = element_blank(), 
        axis.text.x = element_text(angle=45)) +
  xlab("Hour")
ggsave("basset.pdf", plot=p, width=10, height=6, units="in")
