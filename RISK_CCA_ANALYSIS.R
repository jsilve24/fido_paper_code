library(tidyverse)
library(phyloseq)
library(MicrobeDS)
library(stray)

set.seed(899)

setwd("~/Research/mongrel/stray_paper_code/") # Sorry for hardcoding this!


source("src/fit_methods.R")
source("src/dataset_methods.R")
source("src/utils.R")
source("src/fit_stan.R")
source("src/plotting.R")
source("src/fit_mongrel.R")
source("src/GH_standalone.R")



data("RISK_CCFA")
# drop any super low abundant taxa and samples
dat <- RISK_CCFA %>% 
  subset_samples(disease_stat!="missing", 
                 immunosup!="missing") %>% 
  #subset_samples(diseasesubtype %!in% c("UC", "control")) %>%
  subset_samples(diagnosis %in% c("no", "CD")) %>% 
  subset_samples(steroids=="false") %>% 
  subset_samples(antibiotics=="false") %>% 
  subset_samples(biologics=="false") %>% 
  subset_samples(biopsy_location=="Terminal ileum") %>% 
  tax_glom("Family") %>% 
  prune_samples(sample_sums(.) >= 5000,.) %>%
  filter_taxa(function(x) sum(x > 3) > 0.10*length(x), TRUE)

sum(otu_table(dat)==0)/prod(dim(otu_table(dat)))

#Create Design Matrix and OTU Table
sample_dat <- as.data.frame(as(sample_data(dat),"matrix")) %>% 
  mutate(age = as.numeric(as.character(age)),
         diagnosis = relevel(diagnosis, ref="no"), 
         disease_stat = relevel(disease_stat, ref="non-inflamed"))
X <- t(model.matrix(~diagnosis + disease_stat+age, data=sample_dat))
Y <- otu_table(dat)


# Pick Priors
# Set Dimensions
N <- ncol(X)
D <- nrow(Y)
Q <- nrow(X)

# Set Priors
upsilon <- D+3
Theta <- matrix(0, D-1, Q)
Gamma <- diag(Q)
G <- cbind(diag(D-1), -1) ## alr log-constrast matrix
Xi <- 0.5*G%*%diag(D)%*%t(G) ## default is iid on base scale
Xi <- (upsilon-D)*Xi


# View Posteriors -------------------------------------------------------------

# Compile mongrel against MKL for increased numerical stability in this analysis
fit <- mongrel(Y, X)

p <- ppc(fit)
ggsave("output/ppc_real_data.pdf", plot=p, height=5, width=5, units="in")

tax <- as(tax_table(dat)[,c("Class", "Family")], "matrix")
tmp.names <- rownames(tax)
tax <- paste(tax[,"Class"], tax[,"Family"], sep="_")
names(tax) <- tmp.names; rm(tmp.names)
names_categories(fit) <- unname(tax[names_categories(fit)])

cov.names <- c("Intercept", "CD", "Inflamed","Age")
names_covariates(fit) <- cov.names


fit.clr <- mongrel_to_clr(fit)

# # Small values - more stably convert LAMBDA to CLR
# A <- create_alr_base(ncategories(fit), ncategories(fit), inv=TRUE)
# B <- create_clr_base(ncategories(fit))
# B <- B%*%A
# for (i in 1:fit.clr$iter){
#   fit.clr$Lambda[,,i] <- B%*%fit$Lambda[,,i]
#   
# }
summary.clr <- summary(fit.clr, pars="Lambda")

focus <- summary.clr$Lambda %>% 
  filter(sign(p2.5) == sign(p97.5),
         covariate == "CD") %>% 
  arrange(abs(mean))
print(focus, n=50)
write.table(focus, "output/real_summary_nonzero_clr.tsv")

# Create plot of taxa that have non-zero-covering regions
focus.coord <- unique(focus$coord)
focus.covariates <- rownames(X)
p <- plot(fit.clr, par="Lambda", focus.coord=focus.coord, 
          focus.cov = c("CD")) +
  geom_vline(xintercept = 0, alpha=0.3, color="red")
ggsave("output/real_posterior_nonzero_clr.pdf", plot=p, height=4, width=6, units="in")

# Check what prior looked like for these... 
prior.clr <- sample_prior(fit.clr)
p <- plot(prior.clr, par="Lambda", focus.coord=focus.coord, 
          focus.cov = c("CD")) +
  geom_vline(xintercept = 0, alpha=0.3, color="red")
ggsave("output/real_prior_nonzero_clr.pdf", plot=p, height=7, width=10, units="in")



# Investigate Differences -------------------------------------------------

grab_taxa_counts <- function(phyloseq, family){
  bar <- (as(tax_table(phyloseq), "matrix")[,"Family"]== family)
  foo <- prune_taxa(bar, phyloseq) 
  baz <- prune_taxa(!bar, phyloseq)
  data.frame(diagnosis= sample_data(foo)$diagnosis, 
             counts = c(otu_table(foo)))
}


list("Veillonellaceae" = grab_taxa_counts(dat, "Veillonellaceae"), 
     "Peptostreptococcaceae"= grab_taxa_counts(dat, "Peptostreptococcaceae")) %>% 
  bind_rows(.id="Taxa") %>% 
  mutate(diagnosis = ifelse(diagnosis=="CD", "CD", "Healthy" ), 
         counts=counts+1) %>% 
  ggplot(aes(x = diagnosis, y = counts)) +
  geom_boxplot(position = position_dodge()) +
  scale_y_log10() +
  ylab("Counts") +
  theme_minimal() +
  facet_grid(~Taxa) +
  theme(legend.title=element_blank(), 
        axis.title.x = element_blank())
ggsave("output/compare_to_prior_analysis.pdf", height=5, width=4, units="in")



 # Smaller for Stan --------------------------------------------------------


# drop any super low abundant taxa and samples
ids <- sample(unique(as.character(sample_data(dat)$host_subject_id)), 75)
dat <- subset_samples(dat, host_subject_id %in% ids)


sum(otu_table(dat)==0)/prod(dim(otu_table(dat)))

#Create Design Matrix and OTU Table
sample_dat <- as.data.frame(as(sample_data(dat),"matrix")) %>% 
  mutate(age = as.numeric(as.character(age)),
         diagnosis = relevel(diagnosis, ref="no"), 
         disease_stat = relevel(disease_stat, ref="non-inflamed"))
X <- t(model.matrix(~diagnosis + disease_stat+age, data=sample_dat))
Y <- otu_table(dat)


# Pick Priors
# Set Dimensions
N <- ncol(X)
D <- nrow(Y)
Q <- nrow(X)

# Set Priors
upsilon <- D+3
Theta <- matrix(0, D-1, Q)
Gamma <- diag(Q)
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
cov.names <- c("Intercept", "CD", "Inflamed","Age")
p <- bind_rows(gather_array(fit.stancollapsed$Lambda), 
               gather_array(fit.mongrel$Lambda), .id="Implementation") %>% 
  filter(dim_2 %in% 1:4) %>% 
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
ggsave("output/Real_Data_Comparision.pdf", plot=p, height=7, width=5, units="in")  


fit.mongrel$metadata
fit.stancollapsed$metadata