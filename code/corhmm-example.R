setwd("~/dentist-paper/")

require(dentist)
require(corHMM)
# ?corHMM

data(primates)
phy <- multi2di(primates[[1]])
phy$edge.length <- phy$edge.length + 1e-4 
data <- primates[[2]]
data <- data[,c(1,2)]
MK_3state <- corHMM(phy = phy, data = data, rate.cat = 1)
MK_3state

par <- c(1,1,1,1)

fn_corHMM <- function(par, phy, data, rate.cat){
  corhmm_fit <- corHMM(phy = phy, data = data, rate.cat = 1, p = par)
  loglik <- corhmm_fit$loglik
  neg_loglik <- -loglik
  return(neg_loglik)
}

MK_3state$solution

par <- c(0.01689, 0.006224)
names(par) <- c("rate_21", "rate_12")

# par <- 0.02293632
# names(par) <- "rate"

# corhmm_example <- dent_walk(par, fn_corHMM, -MK_3state$loglik, phy = phy, data = data, rate.cat = 1)
# save(corhmm_example, file = "saves/corhmm_example.Rsave")
# 
load("saves/corhmm_example.Rsave")
all_results <- corhmm_example$results[-1,]
accepted_results <- corhmm_example$results[-1,][corhmm_example$acceptances,]
round(cor(all_results[,-1]), 3)
round(cor(accepted_results[,-1]), 3)


dentist:::summary.dentist(corhmm_example)

# save the plot as a pdf
pdf(file = "plots/corhmm_example.pdf", width = 10, height = 10)
plot(corhmm_example)
dev.off()
