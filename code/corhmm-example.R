# functions
quickSim <- function(phy, rate_12=0.1, rate_21=0.1, anc = c(0,1)){
  Q <- matrix(c(-rate_21, rate_12, rate_21, -rate_12), 2, 2)
  data <- corHMM:::simMarkov(phy, Q, anc)$TipStates
  data <- data.frame(sp = names(data), d = data)
  return(data)
}

paraBoot <- function(phy, data, nboots = 100){
  corhmm_fit <- corHMM(phy = phy, data = data, rate.cat = 1)
  loglik <- corhmm_fit$loglik
  neg_loglik <- -loglik
  anc <- corhmm_fit$states[1,]
  Q <- corhmm_fit$solution
  diag(Q) <- -rowSums(Q, na.rm = TRUE)
  boot_dat <- lapply(1:nboots, function(x) corHMM:::simMarkov(phy, Q, anc)$TipStates)
  quick_fit <- function(phy, dat){
    dat <- data.frame(sp = names(dat), d = dat)
    corhmm_fit <- corHMM(phy = phy, data = dat, rate.cat = 1, node.states = "none")
    return(corhmm_fit)
  }
  fits <- mclapply(boot_dat, function(x) try(quick_fit(phy, x)), mc.cores = 10)
  fits <- fits[unlist(lapply(fits, class)) == "corhmm"] # remove univariate simulations
  refit_pars <- do.call(rbind, lapply(fits, function(x) x$solution[!is.na(x$solution)]))
  true_ci <- apply(refit_pars, 2, function(x) quantile(x, c(0.025, 0.5, .975)))
  return(true_ci)
}

fn_corHMM <- function(par, phy, data, rate.cat){
  corhmm_fit <- corHMM(phy = phy, data = data, rate.cat = 1, p = par)
  loglik <- corhmm_fit$loglik
  neg_loglik <- -loglik
  return(neg_loglik)
}

singleRun <- function(nsteps, tree, data){
  corhmm_fit <- corHMM(phy = tree, data = data, rate.cat = 1)
  loglik <- corhmm_fit$loglik
  best_neglnL <- -loglik
  best_par <- corhmm_fit$solution[!is.na(corhmm_fit$solution)]
  names(best_par) <- c("rate_12", "rate_21")
  true_ci <- paraBoot(tree, data, nboots = 1000)
  dent_list <- lapply(nsteps, function(x) dent_walk(par=best_par, fn=fn_corHMM, best_neglnL=best_neglnL, nsteps=x, print_freq=1e10, phy = tree, data = data))
  estimtates <- do.call(cbind, lapply(dent_list, function(x) x$all_ranges[1:3,]))
  true_table <- data.frame(method = "parametric-bootstrap", paramater = c("rate_12", "rate_21"), gen_value = c(0.1, 0.01), best = true_ci[2,], lower.CI = true_ci[1,], upper.CI = true_ci[3,], row.names = NULL)
  dentist_table <- data.frame(method = paste0("dentist_",rep(nsteps, each = 2)), paramater = c("rate_12", "rate_21"), gen_value = c(0.1, 0.01), t(estimtates), row.names = NULL)
  out <- rbind(true_table, dentist_table)
  return(out)
}

# run
setwd("~/dentist-paper/")

require(dentist)
require(corHMM)
require(parallel)
require(TreeSim)

# well get our pars from here
pars <- c(0.1, 0.1)
ntaxa <- 50
age <- 50

set.seed(1985)
trees <- sim.bd.taxa.age(ntaxa, 100, 0.1, 0.05, age = age)
datas <- lapply(trees, function(x) quickSim(x, pars[1], pars[2]))

all_data <- list()
for(i in 1:length(datas)){
  all_data[[i]] <- list(tree = trees[[i]], data = datas[[i]])
}

# parametric boostrap approach
nsteps <- c(10, 50, 100, 500, 1000)

# fits <- mclapply(all_data, function(x) try(singleRun(nsteps, x[[1]], x[[2]])), mc.cores = 20)
# save(fits, file = "saves/corhmm-example-fits.rsave")
load(file = "saves/corhmm-example-fits.rsave")
fits <- fits[unlist(lapply(fits, class)) != "try-error"] # remove univariate simulations
many_sims <- do.call(rbind, fits)
# Add simulation identifier
many_sims$simulation <- rep(1:(nrow(many_sims) / 12), each = 12)

library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)

# Calculate absolute differences
many_sims_diff <- many_sims %>%
  group_by(simulation, paramater, gen_value) %>%
  mutate(
    true_lower.CI = lower.CI[method == "parametric-bootstrap"],
    true_upper.CI = upper.CI[method == "parametric-bootstrap"],
    diff_lower = sqrt((log(true_lower.CI+1) - log(lower.CI+1))^2),
    diff_upper = sqrt((log(true_upper.CI+1) - log(upper.CI+1))^2)
  ) %>%
  ungroup()

many_sims_diff <- many_sims_diff[!many_sims_diff$method == "parametric-bootstrap",]

many_sims_diff$nsteps <- as.numeric(gsub("dentist_", "", many_sims_diff$method))

# Reshape data to long format
many_sims_diff_long <- many_sims_diff %>%
  pivot_longer(cols = c(diff_lower, diff_upper),
               names_to = "bound",
               values_to = "abs_diff")

# Change bound names to more readable format
many_sims_diff_long$bound <- recode(many_sims_diff_long$bound, 
                                    abs_diff_lower = "Lower Bound", 
                                    abs_diff_upper = "Upper Bound")

# Plot
ggplot(many_sims_diff_long, aes(x = as.factor(nsteps), y = (abs_diff), fill = bound)) +
  xlab("Number of steps taken") +
  ylab("Mean Squared Logarithmic Error") +
  geom_boxplot(outlier.shape = NA, width = 0.5) +
  theme_classic() +
  scale_fill_brewer(palette = "Set1") +
  coord_cartesian(ylim = c(0, 1.5)) +
  facet_wrap(~paramater)

ggsave("plots/corhmm-bootstrap-vs-dentist.pdf", width = 10, height = 5)

# table 1 results
dent_1000_12 <- many_sims[(many_sims$method == "dentist_1000") &
                            (many_sims$paramater == "rate_12"),]
para_boot_12 <- many_sims[many_sims$method == "parametric-bootstrap" &
                            (many_sims$paramater == "rate_12"),]

dent_1000_21 <- many_sims[many_sims$method == "dentist_1000" &
                            (many_sims$paramater == "rate_21"),]
para_boot_21 <- many_sims[many_sims$method == "parametric-bootstrap" &
                            (many_sims$paramater == "rate_21"),]

quick_test <- function(focal_row){
  res <- (focal_row[4] > focal_row[5]) & (focal_row[4] < focal_row[6])
  return(res)
}

table(apply(para_boot_12, 1, quick_test))
table(apply(para_boot_21, 1, quick_test))

table(apply(dent_1000_12, 1, quick_test))
table(apply(dent_1000_21, 1, quick_test))

