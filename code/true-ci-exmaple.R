setwd("~/dentist-paper/")

require(reshape2)
require(dentist)
require(ggplot2)

## simulations which show that dentist converges on true CI for a lognormal

# euclidean distance
euclidean <- function(a, b) sqrt(sum((a - b)^2))

# fit a log normal
dlnorm_to_run <- function(par, sims) {
  return(-sum(stats::dlnorm(sims, meanlog=par[1], sdlog=par[2], log=TRUE)))
}

dlnorm_to_run_uni <- function(par, sdlog, sims) {
  return(-sum(stats::dlnorm(sims, meanlog=par, sdlog=sdlog, log=TRUE)))
}


# single simulation which results in a distance calculation for all number of steps
singleRun <- function(npoints, nsteps){
  sims <- stats::rlnorm(npoints, meanlog=1, sdlog=3)
  optimized_results <- stats::optim(c(meanlog=.5, sdlog=1), dlnorm_to_run, sims=sims)
  best_par <- optimized_results$par
  best_neglnL <- optimized_results$value
  # print(best_par)
  sem <- best_par[2]/sqrt(length(sims))
  se <- sem*1.96
  true_ci <- c(lower = best_par[1] - (se), upper = best_par[1] + (se))
  dent_list <- lapply(nsteps, function(x) dent_walk(par=best_par, fn=dlnorm_to_run, best_neglnL=best_neglnL, nsteps=x, print_freq=1e10, sims=sims))
  dent_list <- lapply(nsteps, function(x) dent_walk(par=best_par[1], fn=dlnorm_to_run_uni, best_neglnL=best_neglnL, nsteps=x, print_freq=1e10, sims=sims, sdlog = best_par[2]))
  estimtaes <- do.call(rbind, lapply(dent_list, function(x) c(x$all_ranges[2,1], x$all_ranges[3,1])))
  lower.bound.diff <- estimtaes[,1] - true_ci[1]
  upper.bound.diff <- estimtaes[,2] - true_ci[2]
  # distances <- apply(estimtaes, 1, function(x) euclidean(true_ci, x))
  # names(distances) <- nsteps
  # names(lower.bound.diff) <- names(upper.bound.diff) <- nsteps
  differences <- data.frame(steps = nsteps, lower = lower.bound.diff, upper = upper.bound.diff)
  return(differences)
}

# number of steps for dentist to take
nsteps <- c(10, 50, 100, 500, 1000, 5000)
# number of data points
npoints <- 100

# final data frame
# many_sims <- do.call(rbind, lapply(1:100, function(x) singleRun(npoints, nsteps)))
# save(many_sims, file = "2022_dentist/saves/true-example.Rsave")

setwd("2022_dentist")
load("saves/true-example.Rsave")
plot_df <- melt(many_sims, id.vars = c("steps"))

# plotting the results as a simple boxplot
ggplot(plot_df, aes(x = as.factor(steps), y = abs(value), fill = as.factor(variable))) +
  xlab("Number of steps taken") +
  ylab("Absolute distance to closed form CI") +
  geom_boxplot(outlier.shape = NA, width=0.5) +
  theme_classic()
# save the plot as a pdf
# ggsave("figures/raw/true-example.pdf", width = 6, height = 4)

# tmp <- stats::rlnorm(100000, meanlog=1, sdlog=3)
# sort(tmp)[c(0.025 * 100000, 0.975 * 100000)]

# Optimize the model given the empirical data. We guess at the starting values


# log norm
# sem <- best_par[2]^(1/sqrt(length(sims)))
# se <- best_par[1]/(sem^1.96)
# # normal
# 
# true_ci
# estimtaes
# 
# euclidean
# 
# 
# out <- list()
# for(i in 1:100){
# }
# 
# # boxplot(do.call(rbind, out))
# plot(y = , x = nsteps, xlab = "No. of steps", ylab = "RMSE", main = "Error in SE estimates of the mean")
# abline(h = 0)
# 
# 
# require(corHMM)
# 
# data(primates)
# phy <- multi2di(primates[[1]])
# data <- primates[[2]]
# MK_3state <- corHMM(phy = phy, data = data, rate.cat = 1)
# 
# cor_list <- lapply(nsteps, function(x) ComputeCI(MK_3state, desired.delta = 2, x))
# estimtaes <- do.call(rbind, lapply(cor_list, function(x) c(x$all_ranges[2,], x$all_ranges[3,])))
# 
# # confidence_results <- ComputeCI(MK_3state, desired.delta = 2, 200)
# print(confidence_results)
# corHMM:::plot.dentist(confidence_results)

