rm(list=ls())

#functions
convert2Lambda <- function(pars){
  if(is.na(pars[1])){
    focal_pars <- sample(which(!is.na(pars)), size = 2, replace = FALSE)
    if(2 %in% focal_pars & 3 %in% focal_pars){
      # mu and div
      lambda <- pars[2] + pars[3]
    }
    if(2 %in% focal_pars & 4 %in% focal_pars){
      # mu and turn
      lambda <- pars[4] - pars[2]
    }
    if(2 %in% focal_pars & 5 %in% focal_pars){
      # mu and ef
      lambda <- pars[2]/pars[5]
    }
    if(3 %in% focal_pars & 4 %in% focal_pars){
      # div and turn
      lambda <- (pars[3]+pars[4])/2
    }
    if(3 %in% focal_pars & 5 %in% focal_pars){
      # div and ef
      lambda <- pars[3]/(1-pars[5])
    }
    if(4 %in% focal_pars & 5 %in% focal_pars){
      # turn and ef
      lambda <- pars[4]/(1+pars[5])
    }
  }else{
    lambda <- pars[1]
  }
  return(lambda)
}

convert2Mu <- function(pars){
  if(is.na(pars[2])){
    focal_pars <- sample(which(!is.na(pars)), size = 2, replace = FALSE)
    if(1 %in% focal_pars & 3 %in% focal_pars){
      # lambda and div
      mu <- pars[1] - pars[3]
    }
    if(1 %in% focal_pars & 4 %in% focal_pars){
      # lambda and turn
      mu <- pars[4] - pars[1]
    }
    if(1 %in% focal_pars & 5 %in% focal_pars){
      # lambda and ef
      mu <- pars[1]*pars[5]
    }
    if(3 %in% focal_pars & 4 %in% focal_pars){
      # div and turn
      mu <- (pars[4] - pars[3])/2
    }
    if(3 %in% focal_pars & 5 %in% focal_pars){
      # div and ef
      mu <- (pars[3]*pars[5])/(1-pars[5])
    }
    if(4 %in% focal_pars & 5 %in% focal_pars){
      # turn and ef
      mu <- (pars[4]*pars[5])/(1+pars[5])
    }
  }else{
    mu <- pars[2]
  }
  return(mu)
  
}

convertBetweenPars <- function(pars){
  # pars <- c("lambda", "mu", "net.div", "turn", "ef")
  if(length(which(!is.na(pars))) >= 3){
    warning("More than 2 paramaters are specified. Randomly choosing 2 for the calculations.")
  }
  if(is.na(pars[1])){
    lambda <- convert2Lambda(pars)
  }else{
    lambda <- pars[1]
  }
  if(is.na(pars[2])){
    mu <- convert2Mu(pars)
  }else{
    mu <- pars[2]
  }
  net.div <- lambda - mu
  turn <- lambda + mu
  ef <- mu/lambda
  out <- c(lambda=lambda, mu=mu, net.div=net.div, turn=turn, ef=ef)
  if(!setequal(round(out[which(!is.na(pars))], 5), round(pars[which(!is.na(pars))], 5))){
    stop("An error occured because the calculated output doesn't match the input. Please check that your input parameters can be combined in a way that is possible.")
  }
  return(out)
}

bd_fn <- function(par, phy){
  lik <- make.bd(phy)
  LnLik <- lik(par)
  return(-LnLik)
}

singleRun <- function(nsteps, phy){
  lik <- make.bd(phy)
  fit <- find.mle(lik, c(1, .5))
  loglik <- fit$lnLik
  best_neglnL <- -loglik
  best_par <- fit$par
  age = max(branching.times(phy))
  taxa = length(phy$tip.label)
  
  many_trees <- sim.bd.taxa.age(n = taxa, numbsim = 1000, lambda = best_par[1], mu = best_par[2], age = age)
  many_fits <- lapply(many_trees, function(x) find.mle(make.bd(x), c(1,.5)))
  many_fits <- do.call(rbind, lapply(many_fits, function(x) x$par))
  true_ci <- apply(many_fits, 2, function(x) quantile(x, c(0.025, 0.5, 0.975)))
  true_table <- data.frame(method = "parametric-bootstrap", paramater = c("lambda", "mu"), gen_value = c(1, 0.5), best = true_ci[2,], lower.CI = true_ci[1,], upper.CI = true_ci[3,], row.names = NULL)
  
  dent_list <- lapply(nsteps, function(x) dent_walk(par = best_par, bd_fn, best_neglnL = best_neglnL, nsteps = x, phy=phy))
  estimtates <- do.call(cbind, lapply(dent_list, function(x) x$all_ranges[1:3,]))
  dentist_table <- data.frame(method = paste0("dentist_",rep(nsteps, each = 2)), paramater = c("lambda", "mu"), gen_value = c(1, 0.5), t(estimtates), row.names = NULL)
  out <- rbind(true_table, dentist_table)
  return(out)
}

# run
setwd("~/dentist-paper/")
set.seed(1)

require(diversitree)
require(dentist)

### #### ##### ### #### ##### ### #### #####
# testing CI convergence
### #### ##### ### #### ##### ### #### #####
trees <- sim.bd.taxa.age(n = 20, numbsim = 100, lambda = 1, mu = 0.5, age = 3)

nsteps <- c(10, 50, 100, 500, 1000)
fits <- lapply(trees, function(x) singleRun(nsteps, x))
# save(fits, file = "saves/bd-example-fits.rsave")
load(file = "saves/bd-example-fits.rsave")

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
    abs_diff_lower = abs(true_lower.CI - lower.CI),
    abs_diff_upper = abs(true_upper.CI - upper.CI)
  ) %>%
  ungroup()

many_sims_diff <- many_sims_diff[!many_sims_diff$method == "parametric-bootstrap",]

many_sims_diff$nsteps <- as.numeric(gsub("dentist_", "", many_sims_diff$method))

# Reshape data to long format
many_sims_diff_long <- many_sims_diff %>%
  pivot_longer(cols = c(abs_diff_lower, abs_diff_upper),
               names_to = "bound",
               values_to = "abs_diff")

# Change bound names to more readable format
many_sims_diff_long$bound <- recode(many_sims_diff_long$bound, 
                                    abs_diff_lower = "Lower Bound", 
                                    abs_diff_upper = "Upper Bound")

# Plot
ggplot(many_sims_diff_long, aes(x = as.factor(nsteps), y = abs_diff, fill = bound)) +
  xlab("Number of steps taken") +
  ylab("Absolute distance to bootstrap CI") +
  geom_boxplot(outlier.shape = NA, width = 0.5) +
  theme_classic() +
  scale_fill_brewer(palette = "Set1") +
  ylim(c(0, 2)) +
  facet_wrap(~paramater)


### #### ##### ### #### ##### ### #### ##### 
## practical identifiability and the like
### #### ##### ### #### ##### ### #### #####
## the maximum likelihood speciation/extinction parameters:
phy <- trees(c(1, .5), "bd", max.taxa=10)[[1]]
lik <- make.bd(phy)

## By default, optimisation gives a lambda close to 0.1 and extremely
## small mu:
fit <- find.mle(lik, c(1, .5))

# dent_res_10 <- dent_walk(par = coef(fit), bd_fn, best_neglnL = -fit$lnLik, nsteps = 1000, phy=phy, sd = c(1, 0.5))
# save(dent_res_10, file="saves/dent_res_10.Rsave")
load("saves/dent_res_10.Rsave")

# save the plot as pdf
pdf("plots/bd-example-1.pdf", width=10, height=10)
dentist:::plot.dentist(dent_res_10)
dev.off()

# add more data
phy <- trees(c(1, .5), "bd", max.taxa=1000)[[1]]
lik <- make.bd(phy)

## By default, optimisation gives a lambda close to 0.1 and extremely
## small mu:
fit <- find.mle(lik, c(1, .5))
# dent_res_1000 <- dent_walk(par = coef(fit), bd_fn, best_neglnL = -fit$lnLik, nsteps = 2000, phy=phy)
# save(dent_res_1000, file="saves/dent_res_1000.Rsave")

load("saves/dent_res_1000.Rsave")
# save the plot as pdf
pdf("plots/bd-example-2.pdf", width=10, height=10)
plot(dent_res_1000)
dev.off()


# set.seed(1)
# phy <- trees(c(.1, .05), "bd", max.taxa=10)[[1]]
# lik <- make.bd(phy)
# brts <- getx(phy)
# 
# dd_fit <- dd_ML(brts = brts, initparsopt = c(.1, .05, 100), idparsopt = c(1:3), ddmodel = 1,
#          cond = 1, tol = c(1E-3,1E-3,1E-4), optimmethod = 'simplex')
# 
# dd_fn <- function(par, brts){
#   LnLik <- dd_loglik(pars1 = par, pars2 = c(100,1,1,1,0,2), brts = brts, missnumspec = 0) 
#   return(-LnLik)
# }
# 
# par <- c(l = dd_fit$lambda, m = dd_fit$mu, k = dd_fit$K)
# -dd_fit$loglik == dd_fn(par, brts) # test it's working
# 
# dent_res <- dent_walk(par = par, fn = dd_fn, best_neglnL = -dd_fit$loglik, brts = brts)
# plot(dent_res, local.only = TRUE)
# 
# 
# set.seed(1)
# phy <- trees(c(.1, .05), "bd", max.taxa=500)[[1]]
# lik <- make.bd(phy)
# brts <- getx(phy)
# 
# dd_fit <- dd_ML(brts = brts, initparsopt = c(.1, .05, 1000), idparsopt = c(1:3), ddmodel = 1,
#                 cond = 1, tol = c(1E-3,1E-3,1E-4), optimmethod = 'simplex')
# 
# dd_fn <- function(par, brts){
#   LnLik <- dd_loglik(pars1 = par, pars2 = c(1000,1,1,1,0,2), brts = brts, missnumspec = 0) 
#   return(-LnLik)
# }
# 
# par <- c(l = dd_fit$lambda, m = dd_fit$mu, k = dd_fit$K)
# -dd_fit$loglik == dd_fn(par, brts) # test it's working
# 
# dent_res <- dent_walk(par = par, fn = dd_fn, best_neglnL = -dd_fit$loglik, brts = brts)
# plot(dent_res)
# 
# 





