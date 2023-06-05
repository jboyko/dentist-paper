library(dentist)
library(jocre)

n <- 20
data <- rnorm(n, mean=5, sd=2)

## Not run:
# special package
moodsd <- csetMV(dat=data, method="mood", alpha=0.05, scale="sd")
plot(moodsd)

norm_to_run <- function(par, data) {
  return(-sum(stats::dnorm(data, mean=par[1], sd=par[2], log=TRUE)))
}

optimized_results <- stats::optim(c(mean=.5, sd=1), norm_to_run, data=data)
best_par <- optimized_results$par
best_neglnL <- optimized_results$value

# dentist joint
dent_run_jnt <- dent_walk(par=best_par, fn=norm_to_run, best_neglnL=best_neglnL, nsteps=1000, print_freq=1e10, data=data, confidence_level = 0.95)

# dentist univariate
dent_run_uni <- dent_walk(par=best_par, fn=norm_to_run, best_neglnL=best_neglnL, nsteps=1000, print_freq=1e10, data=data, delta = 1.92)

# univariate
sample_mean <- mean(data)
sample_var <- var(data)
se_mean <- sd(data) / sqrt(n)
alpha <- 0.05
ci_mean <- sample_mean + qt(c(alpha/2, 1-alpha/2), df=n-1) * se_mean
ci_sd <- sqrt((n-1) * sample_var / qchisq(c(1-alpha/2, alpha/2), df=n-1))

# Print the confidence intervals
print(paste0("95% confidence interval for the mean: [", ci_mean[1], ", ", ci_mean[2], "]"))
print(paste0("95% confidence interval for the sd: [", ci_sd[1], ", ", ci_sd[2], "]"))

dent_run_jnt
dent_run_uni
summary(moodsd)

