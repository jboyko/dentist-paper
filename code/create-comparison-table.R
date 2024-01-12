# ugly and hardcoded table.... but at least it's clear

# log normal example
load("saves/true-example.Rsave")

td1 <- many_sims[many_sims[,1] == "t-dist",]
dt1 <- many_sims[many_sims[,1] == "dentist_1000",]

a <- rbind(
  colMeans(td1[,-c(1,2)]),
  colMeans(dt1[,-c(1,2)])
)

b <- rbind(
  td1[1, c(1,2)],
  dt1[1, c(1,2)]
)

lognorm_table <- cbind(model = "lognorm", b, a)

# corhmm example
load(file = "saves/corhmm-example-fits.rsave")
fits <- fits[unlist(lapply(fits, class)) != "try-error"] # remove univariate simulations
many_sims <- do.call(rbind, fits)

pb1 <- many_sims[many_sims[,1] == "parametric-bootstrap" & many_sims[,2] == "rate_12",]
pb2 <- many_sims[many_sims[,1] == "parametric-bootstrap" & many_sims[,2] == "rate_21",]
dt1 <- many_sims[many_sims[,1] == "dentist_1000" & many_sims[,2] == "rate_12",]
dt2 <- many_sims[many_sims[,1] == "dentist_1000" & many_sims[,2] == "rate_21",]

a <- rbind(
  colMeans(pb1[,-c(1,2)]),
  colMeans(pb2[,-c(1,2)]),
  colMeans(dt1[,-c(1,2)]),
  colMeans(dt2[,-c(1,2)])  
)

b <- rbind(
  pb1[1, c(1,2)],
  pb2[1, c(1,2)],
  dt1[1, c(1,2)],
  dt2[1, c(1,2)]
)

corhmm_table <- cbind(model = "corhmm", b, a)

# bd example
load(file = "saves/bd-example-fits.rsave")
fits <- fits[unlist(lapply(fits, class)) != "try-error"] # remove univariate simulations
many_sims <- do.call(rbind, fits)

pb1 <- many_sims[many_sims[,1] == "parametric-bootstrap" & many_sims[,2] == "lambda",]
pb2 <- many_sims[many_sims[,1] == "parametric-bootstrap" & many_sims[,2] == "mu",]
dt1 <- many_sims[many_sims[,1] == "dentist_1000" & many_sims[,2] == "lambda",]
dt2 <- many_sims[many_sims[,1] == "dentist_1000" & many_sims[,2] == "mu",]

a <- rbind(
  colMeans(pb1[,-c(1,2)]),
  colMeans(pb2[,-c(1,2)]),
  colMeans(dt1[,-c(1,2)]),
  colMeans(dt2[,-c(1,2)])  
)

b <- rbind(
  pb1[1, c(1,2)],
  pb2[1, c(1,2)],
  dt1[1, c(1,2)],
  dt2[1, c(1,2)]
)

bd_table <- cbind(model = "bd", b, a)

final_table <- rbind(lognorm_table, corhmm_table, bd_table)

final_table[,-c(1,2,3)] <- round(final_table[,-c(1,2,3)], 5)

write.csv(final_table, file = "compare-table.csv")

