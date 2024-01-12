setwd("~/dentist-paper/")

library(corHMM)
library(phytools)
library(dentist)

# The dataset comes from unpublished empirical work and thus names of species and traits are removed. This dataset is not curated to produce a desired result, but is the most recent example I have of this phenomenon.
# Load data
dat <- read.csv("trait_data.csv")

# Load tree data
tree <- read.tree("tree.tre")

# fit the model
# model_fit <- corHMM(tree, dat, 1)
# save(model_fit, file = "saves/complex-corhmm-model.rsave")
load("saves/complex-corhmm-model.rsave")

fn_corHMM <- function(par, phy, data){
  corhmm_fit <- corHMM(phy = phy, data = data, rate.cat = 1, p = par)
  loglik <- corhmm_fit$loglik
  neg_loglik <- -loglik
  return(neg_loglik)
}

p <- sapply(1:max(model_fit$index.mat, na.rm = TRUE), function(x) 
  na.omit(c(model_fit$solution))[na.omit(c(model_fit$index.mat) == x)][1])

model_fit$index.mat[is.na(model_fit$index.mat)] <- 0
par_names <- c()
for(i in 1:8){
  index <- which(model_fit$index.mat == i, arr.ind = TRUE)
  from <- rownames(model_fit$index.mat)[index[1]]
  to <- colnames(model_fit$index.mat)[index[2]]
  par_names[i] <- paste0(from, "_to_", to)
}
par_names <- gsub("\\(", "", par_names)
par_names <- gsub("\\)", "", par_names)
names(p) <- par_names

# dent_res <- dent_walk(par=p, fn=fn_corHMM, best_neglnL=-model_fit$loglik, nsteps=1000, print_freq=1e10, phy = tree, data = dat)
# save(dent_res, file = "saves/dent_res_complex.rsave")
load("saves/dent_res_complex.rsave")
plot(dent_res)

# dent_res$results[,-1] <- log(dent_res$results[,-1])
# plot(dent_res)
head(dent_res$results)
tmp <- dent_res
tmp$results <- tmp$results[,-c(3,5:9)]
dentist:::plot.dentist(tmp)

# completely randomm stuff i'm doing to examine things
test_mat <- getStateMat4Dat(cor_dat)
focal_model
rate_mat_a <- dropStateMatPars(test_mat$rate.mat, c(2,5))
rate_mat_b <- dropStateMatPars(test_mat$rate.mat, c(6,8))
tmp_a <- corHMM(tree, cor_dat, 1, rate_mat_a)
tmp_b <- corHMM(tree, cor_dat, 1, rate_mat_b)
model_set[[5]] <- tmp_a
model_set[[6]] <- tmp_b
model_table <- corHMM:::getModelTable(model_set)
print(model_table)


plotRECON(tree, tmp_c$states, show.tip.label = FALSE, piecolors = c(1,2,3,4))
tiplabels(pch = 16, col = tmp_c$data.legend$d[match(tree$tip.label, tmp_c$data.legend$sp)], cex = 0.5)

head(tmp_c$data.legend)
head(tree$tip.label)
head()

# our focal candidates
tmp <- corHMM(tree, cor_dat, 1)
# tmp_c <- corHMM(tree, cor_dat, 1, root.p = c(0,0,1,0))
p <- sapply(1:max(tmp$index.mat, na.rm = TRUE), function(x) 
  na.omit(c(tmp_c$solution))[na.omit(c(tmp$index.mat) == x)][1])

library(dentist)
tmp_fixed <- corHMM(tree, cor_dat, 1, tmp_c$index.mat, root.p = c(0,0,1,0), p = tmp_fixed)
ps <- expand.grid(c(1,2,10,15,20,50,75,100), c(1,2,10,15,20,50,75,100))
p_mat <- matrix(rep(p, dim(ps)[1]), nrow = dim(ps)[1], byrow = TRUE)

get_negLnl <- function(par, free_index, p_best){
  p_used <- p_best
  p_used[free_index] <- par
  cor_fixed <- corHMM(tree, cor_dat, 1, tmp_c$index.mat, root.p = c(0,0,1,0), p = p_used)
  return(-cor_fixed$loglik)
}

get_negLnl(c(1,1), c(7,8), p)
par_best <- p[c(6,8)]
names(par_best) <- c("a", "b")

dent_res <- dent_walk(par = par_best, fn = get_negLnl, 
                      best_neglnL = -tmp_c$loglik,
                      free_index = c(6,8), p_best = p)

plot(dent_res)

?dent_walk

