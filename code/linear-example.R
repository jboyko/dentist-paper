rm(list=ls())
set.seed(1)
library(diversitree)
library(dentist)

# example 1
neglnL <- function(par) {
  return(abs(10-par[1]-5*par[2]))
}

lnL <- function(par) {
  return(-neglnL(par)) 
}

samples <- mcmc(lik=lnL, x.init=c(3,1), nsteps=10000, w=5, prior=make.prior.exponential(10))
hist(samples$X1)

plot(X2 ~ X1, samples, pch=19, cex=.2, col="#00000055", asp=1)

plot(X2 ~ X1, samples)


col <- c("blue", "red")



dented_results <- dent_walk(par=c(x=3,y=1), fn=neglnL, best_neglnL = neglnL(c(5,1)), sd=1)
plot(dented_results)
profiles.plot(samples[c("X1", "X2")], col.line=col, las=1,
              xlab="Val 1", legend="topright")
density_values <- density(rexp(1e6, 10))
lines(density_values$x, density_values$y, lwd=2)


# example 2
neglnL <- function(par){
  llik <- abs(10-par[1]-par[2]^2)
  return(llik)
}

lnL <- function(par) {
  return(-neglnL(par)) 
}

rate = 10
prior_fn <- make.prior.exponential(rate)

samples <- mcmc(lik=lnL, x.init=c(9,1), nsteps=10000, w=5, prior=prior_fn)

hist(samples$X1)

plot(X2 ~ X1, samples, pch=19, cex=.2, col="#00000055", asp=1)

plot(X2 ~ X1, samples)


col <- c("blue", "red")



dented_results <- dent_walk(par=c(x=9,y=1), fn=neglnL, best_neglnL = neglnL(c(9,1)), lower_bound = -Inf, upper_bound = Inf, sd = 1)
plot(dented_results)
profiles.plot(samples[c("X1", "X2")], col.line=col, las=1,
              xlab="Parameter value", legend="topright")

density_values <- density(rexp(1e6, rate))
lines(density_values$x, density_values$y, lwd=2)
