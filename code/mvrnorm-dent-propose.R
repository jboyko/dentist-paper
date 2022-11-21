require(MASS)

#' Propose new values
#' This proposes new values using a normal distribution centered on the original parameter values, with desired standard deviation. If any proposed values are outside the bounds, it will propose again.
#' @param old_params The original parameter values
#' @param lower_bound Minimum parameter values to try. One for all or a vector of the length of par.
#' @param upper_bound Maximum parameter values to try. One for all or a vector of the length of par.
#' @param sd Standard deviation to use for the proposals. One for all or a vector of the length of par.
#' @return A vector of the new parameter values
dent_propose <- function(old_params, lower_bound=-Inf, upper_bound=Inf, sd=1) {
  sd <- abs(sd)
  if(runif(1)<1.1) { #try changing all
    new_params <- stats::rnorm(length(old_params), old_params, sd)
  } else { #try sampling some but not all. Usually just one.
    new_params <- old_params
    focal <- sample.int(length(old_params),min(length(old_params), ceiling(stats::rexp(1, 1/2))))
    new_params[focal] <- stats::rnorm(1, old_params[focal], ifelse(length(sd)==1,sd, sd[focal]))  
  }
  while(any(new_params<lower_bound) | any(new_params>upper_bound)) {
    sd <- sd*0.1
    new_params <- dent_propose(old_params, lower_bound=lower_bound, upper_bound=upper_bound, sd=sd)
  }
  return(new_params)
}

#' Propose new values multivariate normal
#' This proposes new values using a normal distribution centered on the original parameter values, with desired standard deviation. If any proposed values are outside the bounds, it will propose again.
#' @param old_params The original parameter values
#' @param lower_bound Minimum parameter values to try. One for all or a vector of the length of par.
#' @param upper_bound Maximum parameter values to try. One for all or a vector of the length of par.
#' @param Sigma VCV to use for the proposals. NULL sets variance to 1 
#' @return A vector of the new parameter values
dent_propose_mv <- function(old_params, lower_bound=-Inf, upper_bound=Inf, Sigma=NULL) {
  if(is.null(Sigma)){
    Sigma <- matrix(0, length(old_params), length(old_params))
    diag(Sigma) <- 1
  }
  new_params <- MASS::mvrnorm(1, old_params, Sigma)
  while(any(new_params<lower_bound) | any(new_params>upper_bound)) {
    diag(Sigma) <- (sqrt(diag(Sigma)) * 0.1)^2
    new_params <- dent_propose_mv(old_params, lower_bound=lower_bound, upper_bound=upper_bound, Sigma=Sigma)
  }
  return(new_params)
}

old_params <- c(1,1)
standard_df <- do.call(rbind, lapply(1:1000, function(x) dent_propose(old_params)))
mvnorm_df <- do.call(rbind, lapply(1:1000, function(x) dent_propose_mv(old_params, Sigma=matrix(c(10,3,3,2),2,2))))
par(mfrow=c(1,2))
plot(standard_df)
plot(mvnorm_df)

Sigma <- matrix(c(10,3,3,2),2,2)
Sigma
var(mvrnorm(n = 1000, rep(0, 2), Sigma))
var(mvrnorm(n = 1000, rep(0, 2), Sigma, empirical = TRUE))



