require(dentist)
require(pdqr)


fn_xy <- function(par, value){
  llik <- ((par[1] + (par[2]^2)) - value)^2
  return(llik)
}

getCube <- function(results){
  volumes <- vector("numeric", dim(results)[1]-1)
  for(i in 2:dim(results)[1]){
    new_volume <- prod(apply(results[1:i,-1], 2, function(x) diff(range(x))))
    volumes[i-1] <- new_volume
  }
  return(volumes)
}

fn_xy(c(9,1), 10)

par <- setNames(c(9,1), c("x", "y"))

simple_xy <- dent_walk(par, fn_xy, 0, value = 10, nsteps = 5000, lower_bound = -Inf, upper_bound = Inf, sphere_probability = 0, sd = 1)
save(simple_xy, file = "saves/simple_xy.Rsave")
# sphere_0.5 <- dent_walk(par, fn_xy, 0, value = 10, nsteps = 5000, lower_bound = -Inf, upper_bound = Inf, sphere_probability = 0.15, sd = 1)
# sphere_1.0 <- dent_walk(par, fn_xy, 0, value = 10, nsteps = 5000, lower_bound = -Inf, upper_bound = Inf, sphere_probability = 1, sd= 1)


load("saves/simple_xy.Rsave")
dentist:::plot.dentist(simple_xy)

x <- seq(from=8, to=-50, length.out=1000)
y1 <- sqrt(8-x)
y2 <- -sqrt(8-x)
y3 <- sqrt(12-x)
y4 <- -sqrt(12-x)
lines(x, y1, col = "blue")
lines(x, y2, col = "blue")
lines(x, y3, col = "blue")
lines(x, y4, col = "blue")

# r <- new_r(sphere_0.0$results[,2], "continuous")


# dentist:::plot.dentist(sphere_0.5)
# global_sphere_mat <- read.table("global_sphere_mat.txt")
# points(global_sphere_mat[-1,], col = "red", cex = 0.1)
# dentist:::plot.dentist(sphere_1.0)
# global_sphere_mat <- read.table("global_sphere_mat.txt")
# points(global_sphere_mat[-1,], col = "red", cex = 0.1)

# plot(test, local.only = FALSE)


all_results <- test$results[-1,]
accepted_results <- test$results[-1,][test$acceptances,]
cov(all_results[,-1])
cov(accepted_results[,-1])

require(MASS)
mvrnorm(n = 1, mu = c(1,1), Sigma = matrix(c(1.114, -0.270, -0.270, 20.701), 2, 2))
