rm(list=ls())
setwd("~/2022_dentist/")
set.seed(1)

require(dentist)
require(OUwie)

data(tworegime)

#Plot the tree and the internal nodes to highlight the selective regimes:
select.reg<-character(length(tree$node.label))
select.reg[tree$node.label == 1] <- "black"
select.reg[tree$node.label == 2] <- "red"

run1 <- OUwie(phy = tree,data = trait,model=c("OUM"), get.root.theta = TRUE, algorithm="invert")
pars <- c(t(run1$solution))
pars <- c(pars, run1$theta[,1])
names(pars) <- c("alpha_1", "alpha_2", "sigma.sq_1", "sigma.sq_2", "theta_1", "theta_2", "theta_root")
# some parameters are fixed, so we need to remove them
pars <- pars[c(1,3,5,6,7)]

ouwie_fn <- function(pars, phy, data){
  ouwie_res <- OUwie.fixed(tree,trait,model=c("OUM"), get.root.theta = TRUE, algorithm="three.point", 
                           alpha = pars[c(1,1)], sigma.sq = pars[c(2,2)], theta = pars[3:5])
  lnLik <- ouwie_res$loglik
  neg_lnLik <- -lnLik
  return(neg_lnLik)
}

neg_lnLik <- ouwie_fn(pars, tree, trait)

#ouwie_example <- dent_walk(pars, ouwie_fn, neg_lnLik, phy=tree, data=trait, nsteps = 2000)
#save(ouwie_example, file = "saves/ouwie_example.Rsave")

load("saves/ouwie_example.Rsave")
# save the plot as a pdf
x <- ouwie_example
pdf("plots/ouwie-example.pdf", width=10, height=10)
nparams <- ncol(x$results) - 1
nplots <- nparams + (nparams^2 - nparams)/2
results <- x$results
threshold <- x$best_neglnL + x$delta
results$color <- ifelse(results[, 1] <= threshold, "black", 
    "gray")
results_outside <- subset(results, results$color == "gray")
results_inside <- subset(results, results$color == "black")
graphics::par(mfrow = c(ceiling(nplots/nparams), nparams))
for (i in sequence(nparams)) {
  xlim = range(c(results_inside[, i + 1]), results_outside[,i + 1])
  ylim = c(quantile(c(results_inside[, 1], results_outside[,1]), c(0, 0.95)))
    plot(results[, i + 1], results[, 1], pch = 20, col = results$color, 
        main = colnames(results)[i + 1], xlab = colnames(results)[i + 
            1], ylab = "Negative Log Likelihood", bty = "n", 
        xlim = xlim, ylim = ylim)
    graphics::abline(h = threshold, col = "blue")
    graphics::points(results[which.min(results[, 1]), i + 
        1], results[which.min(results[, 1]), 1], pch = 21, 
        col = "red")
}
for (i in sequence(nparams)) {
    for (j in sequence(nparams)) {
        if (j > i) {
            xlim = quantile(c(results_inside[, i + 1], results_outside[,i + 1]), c(0, 0.99))
            ylim = quantile(c(results_inside[, j + 1], results_outside[,j + 1]), c(0, 0.99))
            plot(results_outside[, i + 1], results_outside[, 
              j + 1], pch = 20, col = results_outside$color, 
              xlab = colnames(results)[i + 1], ylab = colnames(results)[j + 
                1], bty = "n", main = paste0(colnames(results)[j + 
                1], " vs. ", colnames(results)[i + 1]), xlim = xlim, 
              ylim = ylim)
            graphics::points(results_inside[, i + 1], results_inside[, 
              j + 1], pch = 20, col = results_inside$color)
            graphics::points(results[which.min(results[, 
              1]), i + 1], results[which.min(results[, 1]), 
              j + 1], pch = 21, col = "red")
        }
    }
}
dev.off()

# all_results <- dent_res$results[-1,]
# accepted_results <- dent_res$results[-1,][dent_res$acceptances,]
# round(cov(all_results[,-1]), 3)
# round(cov(accepted_results[,-1]), 3)
# 
# 
# plot(tree)
# nodelabels(pch=21, bg=select.reg)
# 
# 
