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
  # dent_list <- lapply(nsteps, function(x) dent_walk(par=best_par, fn=dlnorm_to_run, best_neglnL=best_neglnL, nsteps=x, print_freq=1e10, sims=sims))
  dent_list <- lapply(nsteps, function(x) dent_walk(par=best_par[1], fn=dlnorm_to_run_uni, best_neglnL=best_neglnL, nsteps=x, print_freq=1e10, sims=sims, sdlog = best_par[2]))
  estimtaes <- do.call(rbind, lapply(dent_list, function(x) c(x$all_ranges[1,1], x$all_ranges[2,1], x$all_ranges[3,1])))
  true_table <- data.frame(method = "t-dist", paramater = "mean", gen_value = 1, best = best_par[1], lower.CI = true_ci[1], upper.CI = true_ci[2], row.names = NULL)
  dentist_table <- data.frame(method = paste0("dentist_", nsteps), paramater = "mean", gen_value = 1, estimtaes, row.names = NULL)
  out <- rbind(true_table, dentist_table)
  # lower.bound.diff <- estimtaes[,1] - true_ci[1]
  # upper.bound.diff <- estimtaes[,2] - true_ci[2]
  # distances <- apply(estimtaes, 1, function(x) euclidean(true_ci, x))
  # names(distances) <- nsteps
  # names(lower.bound.diff) <- names(upper.bound.diff) <- nsteps
  # differences <- data.frame(steps = nsteps, lower = lower.bound.diff, upper = upper.bound.diff)
  return(out)
}

# number of steps for dentist to take
nsteps <- c(10, 50, 100, 500, 1000)
# number of data points
npoints <- 100

# final data frame
# many_sims <- do.call(rbind, lapply(1:100, function(x) singleRun(npoints, nsteps)))
# save(many_sims, file = "saves/true-example.Rsave")

setwd("dentist-paper")
load("saves/true-example.Rsave")
# plot_df <- melt(many_sims, id.vars = c("steps"))

library(dplyr)
library(tidyr)
library(stringr)

# Add simulation identifier
many_sims$simulation <- rep(1:(nrow(many_sims) / 6), each = 6)

# Calculate absolute differences
many_sims_diff <- many_sims %>%
  group_by(simulation, paramater, gen_value) %>%
  mutate(
    true_lower.CI = lower.CI[method == "t-dist"],
    true_upper.CI = upper.CI[method == "t-dist"],
    abs_diff_lower = abs(true_lower.CI - lower.CI),
    abs_diff_upper = abs(true_upper.CI - upper.CI)
  ) %>%
  ungroup()

many_sims_diff <- many_sims_diff[!many_sims_diff$method == "t-dist",]

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
  ylab("Absolute distance to closed form CI") +
  geom_boxplot(outlier.shape = NA, width = 0.5) +
  theme_classic() +
  scale_fill_brewer(palette = "Set1")

# table 1 results

dent_1000 <- many_sims[many_sims$method == "dentist_1000",]
t_dist <- many_sims[many_sims$method == "t-dist",]

quick_test <- function(focal_row){
  res <- (focal_row[4] > focal_row[5]) & (focal_row[4] < focal_row[6])
  return(res)
}

table(apply(dent_1000, 1, quick_test))
table(apply(t_dist, 1, quick_test))


# save the plot as a pdf
# ggsave("plots/true-example.pdf", width = 6, height = 4)

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

