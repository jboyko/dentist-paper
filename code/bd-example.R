rm(list=ls())
setwd("~/2022_dentist/")
set.seed(1)

require(diversitree)
require(dentist)

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

## Simulate a tree under a constant rates birth-death model and look at
## the maximum likelihood speciation/extinction parameters:
phy <- trees(c(1, .5), "bd", max.taxa=10)[[1]]
lik <- make.bd(phy)

## By default, optimisation gives a lambda close to 0.1 and extremely
## small mu:
fit <- find.mle(lik, c(1, .5))

# remember the change the output of the function to negative loglik
bd_fn <- function(par, phy){
  lik <- make.bd(phy)
  LnLik <- lik(par)
  return(-LnLik)
}


# dent_res_10 <- dent_walk(par = coef(fit), bd_fn, best_neglnL = -fit$lnLik, nsteps = 1000, phy=phy, sd = c(1, 0.5))
#save(dent_res_10, file="saves/dent_res_10.Rsave")
load("saves/dent_res_10.Rsave")

# a quick parametric bootstrap
library(TreeSim)
age = max(branching.times(phy))
taxa = length(phy$tip.label)
many_trees <- sim.bd.taxa.age(n = taxa, numbsim = 1000, lambda = fit$par[1], mu = fit$par[2], age = age)
many_fits <- lapply(many_trees, function(x) find.mle(make.bd(x), c(1,.5)))
many_fits <- do.call(rbind, lapply(many_fits, function(x) x$par))
save(many_fits, file = "saves/many_fits_10.Rsave")
colMeans(many_fits)
quantile(many_fits[,1], c(0.025, 0.975))
quantile(many_fits[,2], c(0.025, 0.975))
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





