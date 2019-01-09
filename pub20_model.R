#devtools::install_github("jsilve24/mongrel")
library(driver)
library(tidyverse)
library(phyloseq)
require(rstan)
require(mongrel)

setwd("~/Research/mongrel/mongrel_paper_code/") # Sorry for hardcoding this!


source("src/fit_methods.R")
source("src/dataset_methods.R")
source("src/utils.R")
source("src/fit_stan.R")
source("src/plotting.R")
source("src/fit_mongrel.R")
source("src/GH_standalone.R")


set.seed(4)


options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

 
# load_data ---------------------------------------------------------------


ps_raw <- readRDS("data/phyloseq.rds")
ps <- subset_samples(ps_raw, X.SampleID != "ZH8")
ps <- subset_samples(ps, X.SampleID != "ZH18")
ps <- subset_samples(ps, X.SampleID != "ZH20")


scfa <- read_csv("data/20PUB_blSCFA_4Mongrel.csv")
scfa <- as.data.frame(scfa)
rownames(scfa) <- scfa[,1]
scfa <- scfa[,-1]
scfa <- as.matrix(scfa)

# Preprocess Data ---------------------------------------------------------

# Amalgamate to Genus Level and filter, taxa must be present with at 
# least three counts in at least 3 samples
ps.genus <- tax_glom(ps, "Genus")
ps.genus <- filter_taxa(ps.genus, function(x) sum(x>=3)>=3, TRUE)
Y <- t(as(otu_table(ps.genus), "matrix"))
X <- t(scfa)
X <- X[,colnames(Y)]
X <- rbind("intercept"=1, X)

# Set Dimensions
N <- ncol(X)
D <- nrow(Y)
Q <- nrow(X)

# Set Priors
upsilon <- D+3
Theta <- matrix(0, D-1, Q)
#OldGamma <- diag(Q)
#Gamma <- OldGamma
#Gamma <- matrix(c(5,0,0,0,0,2,-0.1,-0.2,0,-0.1,2,.3,0,-0.2,.3,2), nrow=4, ncol=4)
#Gamma <- matrix(c(5,0,0,0,0,5,-0.2,-0.4,0,-0.2,5,.4,0,-0.4,.4,5), nrow=4, ncol=4)
Gamma <- matrix(c(5,0,0,0,0,2,0.2,0.2,0,0.2,2,.4,0,0.2,.4,2), nrow=4, ncol=4)
G <- cbind(diag(D-1), -1) ## alr log-constrast matrix
Xi <- 0.5*G%*%diag(D)%*%t(G) ## default is iid on base scale
Xi <- (upsilon-D)*Xi



# Efficiency Comparisons --------------------------------------------------

Lambda_true <- matrix(-1, D-1, Q)  # just to use mdataset without validation
Sigma_true <- matrix(-1, D-1, D-1) #  "  " 

md <- mdataset(N, D, Q, Y, X, Lambda_true, Sigma_true, Theta, Gamma, Xi, upsilon)

fit.mongrel <- fit_mongrel(md, decomposition="cholesky",ncores=4)
fit.stancollapsed <- fit_mstan(md, parameterization="collapsed")
fit.stanuncollapsed <- fit_mstan(md, parameterization="uncollapsed")
fit.vb <- fit_mstan_vb(md, parameterization = "uncollapsed")

#plot_lambda(list(fit.mongrel, fit.stancollapsed))
cov.names <- c("Intercept", "Acetate", 
               "Propionate", "Butyrate")
p <- bind_rows(gather_array(fit.stancollapsed$Lambda), 
          gather_array(fit.mongrel$Lambda), .id="Implementation") %>% 
  rename(Coordinate=dim_1, Covariate=dim_2, iter=dim_3) %>% 
  mutate(Implementation=ifelse(Implementation==1, "HMC Collapsed", "LA Collapsed")) %>% 
  mutate(Covariate = cov.names[Covariate]) %>% 
  mutate(Covariate = factor(Covariate, levels = cov.names)) %>% 
  mutate(Coordinate = factor(Coordinate)) %>% 
  group_by(Implementation, Coordinate, Covariate) %>% 
  summarise_posterior(var) %>% 
  ungroup() %>% 
  ggplot(aes(x=Coordinate, y=mean, color=Implementation)) + 
  geom_linerange(aes(ymin=p2.5, ymax=p97.5), position=position_dodge(width=.3), alpha=0.7) + 
  geom_point(alpha=0.7) +
  facet_grid(~Covariate) +
  coord_flip() +
  theme_minimal() +
  scale_color_brewer(palette = "Set1") + 
  theme(axis.title.x = element_blank(), 
        legend.position="bottom", 
        axis.text.y = element_blank()) 
ggsave("Real_Data_Comparision.pdf", plot=p, height=7, width=5, units="in")  


 
# Analyze Mongrel Model output --------------------------------------------

# Fit Model
fit <- mongrel(Y, X, upsilon, Theta, Gamma, Xi)

p <- ppc(fit)
ggsave("ppc_real_data.pdf", plot=p, height=5, width=5, units="in")


# change names 
tax <- as(tax_table(ps.genus)[,c("Genus", "Family")], "matrix")
tmp.names <- rownames(tax)
tax <- paste(tax[,"Family"], tax[,"Genus"], sep="_")
names(tax) <- tmp.names; rm(tmp.names)
names_categories(fit) <- unname(tax[names_categories(fit)])

# View results in CLR coordinates
fit.clr <- mongrel_to_clr(fit)
summary.clr <- summary(fit.clr, pars="Lambda")

# Focus on just the elements of lambda with posterior 95% credible regions
# not covering zero
focus <- summary.clr$Lambda %>% 
  filter(sign(p2.5) == sign(p97.5),
         covariate != "intercept") %>% 
  arrange(abs(mean))
print(focus, n=50)
write.table(focus, "output/5_2_2_4summary_nonzero_clr.tsv")

# Create plot of taxa that have non-zero-covering regions
focus.coord <- unique(focus$coord)
focus.covariates <- rownames(X)[-1]
p <- plot(fit.clr, par="Lambda", focus.coord=focus.coord, 
          focus.cov = focus.covariates) +
  geom_vline(xintercept = 0, alpha=0.3, color="red")
ggsave("output/5_2_2_4posterior_nonzero_clr.pdf", plot=p, height=5, width=8, units="in")

# Check what prior looked like for these... 
prior.clr <- sample_prior(fit.clr)
p <- plot(prior.clr, par="Lambda", focus.coord=focus.coord, 
          focus.cov = focus.covariates) +
  geom_vline(xintercept = 0, alpha=0.3, color="red")
ggsave("output/5_2_2_4prior_nonzero_clr.pdf", plot=p, height=7, width=10, units="in")

# check posterior predictive fit 
p <- ppc(fit)
ggsave("output/5_2_2_4ppc.pdf", plot=p, height=5, width=8, units="in")

