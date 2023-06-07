# functions
quickSim <- function(phy, rate_12=0.1, rate_21=0.01, anc = c(0,1)){
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
  true_ci <- paraBoot(tree, data, nboots = 100)
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
pars <- c(0.1, 0.01)
ntaxa <- 50
age <- 50

set.seed(1)
trees <- sim.bd.taxa.age(ntaxa, 100, 0.1, 0.05, age = age)
datas <- lapply(trees, function(x) quickSim(x, pars[1], pars[2]))

all_data <- list()
for(i in 1:length(datas)){
  all_data[[i]] <- list(tree = trees[[i]], data = datas[[i]])
}

# parametric boostrap approach
nsteps <- c(10, 50, 100, 500, 1000)

fits <- lapply(all_data, function(x) try(singleRun(nsteps, x[[1]], x[[2]])))
save(fits, file = "saves/corhmm-example-fits.rsave")
load(file = "saves/corhmm-example-fits.rsave")
# fits <- fits[unlist(lapply(fits, class)) == "corhmm"] # remove univariate simulations
# many_sims <- do.call(rbind, lapply(all_data, function(x) singleRun(npoints, nsteps)))

refit_pars <- do.call(rbind, lapply(fits, function(x) x$solution[!is.na(x$solution)]))

# using dentist
par <- c(MK_3state$solution[!is.na(MK_3state$solution)])
names(par) <- c("rate_21", "rate_12")
data <- primates[[2]]
data <- data[,c(1,2)]
# corhmm_example_jnt <- dent_walk(par, fn_corHMM, -MK_3state$loglik, phy = phy, data = data, rate.cat = 1, nsteps = 2000)
# save(corhmm_example_jnt, file = "saves/corhmm-example-corhmm_example_jnt.rsave")
load("saves/corhmm-example-corhmm_example_jnt.rsave")

# comparison
corhmm_example_jnt
mean(refit_pars[,1])
median(refit_pars[,1])
quantile(refit_pars[,1], c(.025, .975))
mean(refit_pars[,2])
median(refit_pars[,2])
quantile(refit_pars[,2], c(.025, .975))

# save(corhmm_example, file = "saves/corhmm_example.Rsave")
# 
# load("saves/corhmm_example.Rsave")
all_results <- corhmm_example_jnt$results[-1,]
accepted_results <- corhmm_example_jnt$results[-1,][corhmm_example_jnt$acceptances,]
round(cor(all_results[,-1]), 3)
round(cor(accepted_results[,-1]), 3)


dentist:::summary.dentist(corhmm_example_jnt)

# save the plot as a pdf
pdf(file = "plots/corhmm_example.pdf", width = 10, height = 10)
plot(corhmm_example_jnt)
dev.off()
