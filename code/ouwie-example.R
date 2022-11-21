require(dentist)
require(OUwie)

data(tworegime)

#Plot the tree and the internal nodes to highlight the selective regimes:
select.reg<-character(length(tree$node.label))
select.reg[tree$node.label == 1] <- "black"
select.reg[tree$node.label == 2] <- "red"



## Not run: 
#To see the first 5 lines of the data matrix to see what how to
#structure the data:

run1 <- OUwie(phy = tree,data = trait,model=c("OUMV"), get.root.theta = TRUE, algorithm="invert")
pars <- c(t(run1$solution))
pars <- c(pars, run1$theta[,1])
names(pars) <- c("alpha_1", "alpha_2", "sigma.sq_1", "sigma.sq_2", "theta_1", "theta_2", "theta_root")

ouwie_fn <- function(pars.free, pars.fixed, phy, data){
  ouwie_res <- OUwie.fixed(tree,trait,model=c("OUMV"), get.root.theta = TRUE, algorithm="three.point", 
                           alpha = pars.fixed, sigma.sq = pars.free[1:2], theta = pars.free[3:5])
  lnLik <- ouwie_res$loglik
  neg_lnLik <- -lnLik
  return(neg_lnLik)
}

neg_lnLik <- ouwie_fn(pars[3:7], pars[1:2], tree, trait)

pars.fixed <- pars[1:2]
pars.free <- pars[3:7]

dent_res <- dent_walk(pars.free, ouwie_fn, neg_lnLik, pars.fixed=pars.fixed, phy=tree, data=trait, nsteps = 1000)
dentist:::plot.dentist(dent_res, local.only=TRUE)

all_results <- dent_res$results[-1,]
accepted_results <- dent_res$results[-1,][dent_res$acceptances,]
round(cov(all_results[,-1]), 3)
round(cov(accepted_results[,-1]), 3)


plot(tree)
nodelabels(pch=21, bg=select.reg)


