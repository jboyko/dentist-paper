setwd("~/dentist-paper/")

require(dentist)
require(corHMM)
require(parallel)
# ?corHMM

data(primates)
phy <- multi2di(primates[[1]])
phy$edge.length <- phy$edge.length + 1e-4 
data <- primates[[2]]
data <- data[,c(1,2)]
MK_3state <- corHMM(phy = phy, data = data, rate.cat = 1)
MK_3state

fn_corHMM <- function(par, phy, data, rate.cat){
  corhmm_fit <- corHMM(phy = phy, data = data, rate.cat = 1, p = par)
  loglik <- corhmm_fit$loglik
  neg_loglik <- -loglik
  return(neg_loglik)
}

quick_fit <- function(phy, dat){
  dat <- data.frame(sp = names(dat), d = dat)
  corhmm_fit <- corHMM(phy = phy, data = dat, rate.cat = 1, node.states = "none")
  return(corhmm_fit)
}

# parametric boostrap approach
anc <- MK_3state$states[1,]
Q <- MK_3state$solution
diag(Q) <- -rowSums(Q, na.rm = TRUE)
data <- lapply(1:1000, function(x) corHMM:::simMarkov(phy, Q, anc)$TipStates)
# fits <- mclapply(data, function(x) try(quick_fit(phy, x)), mc.cores = 10)
# fits <- fits[unlist(lapply(fits, class)) == "corhmm"] # remove univariate simulations
# save(fits, file = "saves/corhmm-example-fits.rsave")
load(file = "saves/corhmm-example-fits.rsave")
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
