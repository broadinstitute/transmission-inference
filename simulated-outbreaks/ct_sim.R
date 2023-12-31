# MIT License
# 
# Copyright (c) 2023 Ivan Specht
# 
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#   
#   The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

### Script to assess accuracy of outbreak reconstruction with contact tracng data

# Relevant libraries
library(ape)
library(insect)
library(stringr)
library(igraph)
library(parallel)
library(Rfast)
library(reconstructR)
library(outbreaker2)
library(ggplot2)
library(BSDA)
library(EnvStats)
library(cowplot)
library(ggpubr)

setwd("~/Desktop/outbreak-reconstruction/simulated-outbreaks/")
#load("./newsims.RData")

get_adj <- function(i){
  setwd("~/Desktop/outbreak-reconstruction/simulated-outbreaks/")
  dirname <- paste0("newdir_", i)
  setwd(paste0("./", dirname))

  ### Run reconstructR

  results <- run_mcmc(
    N_iters = 2e3
  )

  burnin <- 0.1

  N <- results[[1]]$N
  adj <- matrix(0, nrow = N, ncol = N)
  for (j in 1:length(results)) {
    if (j > burnin * length(results)) {
      inds <- cbind(results[[j]]$anc[2:N], 2:N)
      adj[inds] <- adj[inds] + 1
    }
  }

  adj <- adj / sum(1:length(results) > burnin * length(results))

  return(adj)

}

get_adj_o2 <- function(i){

  setwd("~/Desktop/outbreak-reconstruction/simulated-outbreaks/")
  dirname <- paste0("newdir_", i)
  setwd(paste0("./", dirname))

  ### Run outbreaker2

  dna <- read.FASTA("./input_data/aligned.fasta")
  names(dna) <- NULL
  ctd <- read.csv("./input_data/ct.csv")
  if(nrow(ctd) > 1){ # Broken in outbreaker source code if only one CT link!
    ctd[,1] <- (as.numeric(gsub("sim_", "", ctd[,1])))
    ctd[,2] <- (as.numeric(gsub("sim_", "", ctd[,2])))
    dat <- list(
      dna = dna,
      dates = read.csv("./input_data/date.csv")[,2],
      ctd = ctd,
      w_dens = dgamma(1:20, 3)
    )
    out <- outbreaker(
      data = dat,
      config = list(
        n_iter = 2e3,
        sample_every = 100,
        ctd_directed = F
      )
    )
  }else{
    dat <- list(
      dna = dna,
      dates = read.csv("./input_data/date.csv")[,2],
      w_dens = dgamma(1:20, 3)
    )
    out <- outbreaker(
      data = dat,
      config = list(
        n_iter = 2e3,
        sample_every = 100
      )
    )
  }


  res <- as.matrix(out)
  res <- res[, 9:(9+length(dna)-1)]
  res <- res + 1
  res[is.na(res)] <- 1
  adj_o2 <- matrix(0, nrow = length(dna)+1, ncol = length(dna)+1)
  burnin <- 0.1

  for (j in 1:nrow(res)) {
    if (j > burnin * nrow(res)) {
      inds <- cbind(res[j, ], 2:(length(dna)+1))
      adj_o2[inds] <- adj_o2[inds] + 1
    }
  }

  adj_o2 <- adj_o2 / sum(1:nrow(res) > burnin * nrow(res))

  return(adj_o2)
}

nus <- c(0,0.25,0.5,0.75)
xis <- c(0,0.5)

all_adjs <- list()
all_adjs_o2 <- list()

for (xi in xis) {
  for (nu in nus) {

    for(i in 1:length(newsims)){
      true_pairs <- newsims[[i]][[1]]
      true_pairs <- true_pairs[true_pairs[,1] != 1, ]
      people <- sort(unique(as.vector(true_pairs)))
      all_pairs <- t(combn(people, 2))
      true <- c()
      for (j in 1:nrow(all_pairs)) {
        from <- all_pairs[j, 1]
        to <- all_pairs[j, 2]


        if(any(true_pairs[,1] == from & true_pairs[,2] == to)){
          true <- c(true, j)
        }
      }
      false_pairs <- all_pairs[-true, , drop = F]

      ct_true <- true_pairs[runif(nrow(true_pairs)) < nu, , drop = F]
      ct_false <- false_pairs[runif(nrow(true_pairs)) < xi, , drop = F]

      ct <- rbind(ct_true, ct_false)
      ct <- as.data.frame(ct)
      if(nrow(ct) > 0){
        ct[,1] <- str_pad(as.character(ct[,1] - 1), 3, pad = "0")
        ct[,1] <- paste0("sim_", ct[,1])
        ct[,2] <- str_pad(as.character(ct[,2] - 1), 3, pad = "0")
        ct[,2] <- paste0("sim_", ct[,2])
      }


      setwd("~/Desktop/outbreak-reconstruction/simulated-outbreaks/")
      dirname <- paste0("newdir_", i)
      setwd(paste0("./", dirname))

      write.csv(ct, "./input_data/ct.csv", row.names = F)

      # if(nrow(ct) == 1){
      #   print(i)
      # }
    }

    print(paste("xi =", xi))
    print(paste("nu =", nu))
    adj <- mclapply(1:100, get_adj, mc.cores = 12)
    print("reconstructR complete")
    adj_o2 <- mclapply(1:100, get_adj_o2, mc.cores = 12)
    print("outbreaker2 complete")
    all_adjs <- c(all_adjs, list(adj))
    all_adjs_o2 <- c(all_adjs_o2, list(adj_o2))

  }
}

save(all_means, file = "all_means.RData")
save(all_means_o2, file = "all_means_o2.RData")


get_mean <- function(truth, adj, spec = F){
  unrelated <- which(truth[,1] == 1)
  related <- setdiff(1:nrow(truth), unrelated)

  related_acc <- adj[truth[related, ]]

  # Probability that the ancestor of an index or unrelated case ISNT among our cluster
  unrelated_acc <- colsums(adj[c(1, unrelated+1), unrelated+1, drop = F])

  if(!spec){
    mean(related_acc)
  }else{
    mean(unrelated_acc)
  }

  #mean(c(related_acc, unrelated_acc))
}

get_cov <- function(truth, adj){

  out <- c()
  real <- truth
  unrelated <- real[which(real[,1] == 1), 2]
  m <- adj
  for (j in 2:ncol(m)) {
    col <- m[,j] # extract probabilities of possible ancestors for case j
    cred <- 0 # credibility of subset of ancestors; initialized to empty set
    set <- c() # subset of ancestors

    while (cred < 0.9) {
      to_append <- which.max(col)
      set <- c(set, to_append)
      cred <- sum(m[set, j])
      col[to_append] <- 0
    }

    if(j %in% unrelated){
      if(any(unrelated %in% set)){
        out[j] <- 1
      }else{
        out[j] <- 0
      }
    }else{
      if(real[which(real[,2] == j), 1] %in% set){
        out[j] <- 1
      }else{
        out[j] <- 0
      }
    }
  }

  mean(out[2:ncol(m)])


}

## Get hit rate
get_hit_rate <- function(truth, adj){
  unrelated <- which(truth[,1] == 1)
  related <- setdiff(1:nrow(truth), unrelated)
  
  bin_adj <- matrix(0, nrow = nrow(adj), ncol = ncol(adj))
  bin_adj[cbind(colMaxs(adj)[2:ncol(adj)], 2:ncol(adj))] <- 1
  
  related_hit <- bin_adj[truth[related, ]]
  mean(related_hit)
}

all_means <- list()
all_means_o2 <- list()
all_specs <- list()
all_specs_o2 <- list()
all_covs <- list()
all_covs_o2 <- list()
all_hits <- list()
all_hits_o2 <- list()

truths <- list()
for (i in 1:length(newsims)) {
  truths[[i]] <- newsims[[i]][[1]]
}

for (r in 1:(length(xis) * length(nus))) {
  all_means[[r]] <- mapply(get_mean, truths, all_adjs[[r]])
  all_means_o2[[r]] <- mapply(get_mean, truths, all_adjs_o2[[r]])
  all_specs[[r]] <- mapply(get_mean, truths, all_adjs[[r]], MoreArgs = list(spec = T))
  all_specs_o2[[r]] <- mapply(get_mean, truths, all_adjs_o2[[r]], MoreArgs = list(spec = T))
  all_covs[[r]] <- mapply(get_cov, truths, all_adjs[[r]])
  all_covs_o2[[r]] <- mapply(get_cov, truths, all_adjs_o2[[r]])
  all_hits[[r]] <- mapply(get_hit_rate, truths, all_adjs[[r]])
  all_hits_o2[[r]] <- mapply(get_hit_rate, truths, all_adjs_o2[[r]])
}

params <- c()
for (xi in xis) {
  for(nu in nus){
    params <- c(params, paste0("nu=", nu, ", xi=", xi))
  }
}

addline_format <- function(x,...){
  gsub('\\,\\s','\n',x)
}

setwd("~/Desktop/outbreak-reconstruction/simulated-outbreaks/")

acc_plot <- data.frame(
  Accuracy = c(unlist(all_means), unlist(all_means_o2)),
  Parameters = rep(rep(params, each = length(newsims)), 2),
  Model = c(rep("This Paper", length(newsims) * length(xis) * length(nus)), rep("outbreaker2 (2018)", length(newsims) * length(xis) * length(nus)))
)
acc_plot$Parameters <- factor(acc_plot$Parameters, levels = params)
acc_plot$Model <- factor(acc_plot$Model, levels = c("This Paper", "outbreaker2 (2018)"))

ct1 <- ggplot(acc_plot, aes(x=Parameters, y=Accuracy, fill=Model, color = Model)) +
  geom_boxplot() +
  scale_fill_manual(values=c( "#EECCAA", "#AACCEE")) +
  scale_color_manual(values=c("#BB9977", "#7799BB")) +
  scale_x_discrete(labels = addline_format(params)) +
  theme_minimal() +
  theme(legend.position = "none")

print(ct1)

ggsave("./figs/ct_acc.png", width = 12, height = 6)
ggsave("./figs/ct_acc.pdf", width = 12, height = 6)

cov_plot <- data.frame(
  Coverage = c(unlist(all_covs), unlist(all_covs_o2)),
  Parameters = rep(rep(params, each = length(newsims)), 2),
  Model = c(rep("This Paper", length(newsims) * length(xis) * length(nus)), rep("outbreaker2 (2018)", length(newsims) * length(xis) * length(nus)))
)
cov_plot$Parameters <- factor(cov_plot$Parameters, levels = params)
cov_plot$Model <- factor(cov_plot$Model, levels = c("This Paper", "outbreaker2 (2018)"))


ct2 <- ggplot(cov_plot, aes(x=Parameters, y=Coverage, fill=Model, color = Model)) +
  geom_boxplot() +
  scale_fill_manual(values=c( "#EECCAA", "#AACCEE")) +
  scale_color_manual(values=c("#BB9977", "#7799BB")) +
  scale_x_discrete(labels = addline_format(params)) +
  theme_minimal() +
  theme(legend.position = "none")


print(ct2)

ggsave("./figs/ct_cov.png", width = 12, height = 6)
ggsave("./figs/ct_cov.pdf", width = 12, height = 6)


spec_plot <- data.frame(
  Specificity = c(unlist(all_specs), unlist(all_specs_o2)),
  Parameters = rep(rep(params, each = length(newsims)), 2),
  Model = c(rep("This Paper", length(newsims) * length(xis) * length(nus)), rep("outbreaker2 (2018)", length(newsims) * length(xis) * length(nus)))
)
spec_plot$Parameters <- factor(spec_plot$Parameters, levels = params)
spec_plot$Model <- factor(spec_plot$Model, levels = c("This Paper", "outbreaker2 (2018)"))


ct3 <- ggplot(spec_plot, aes(x=Parameters, y=Specificity, fill=Model, color = Model)) +
  geom_boxplot() +
  scale_fill_manual(values=c( "#EECCAA", "#AACCEE")) +
  scale_color_manual(values=c("#BB9977", "#7799BB")) +
  scale_x_discrete(labels = addline_format(params)) +
  theme_minimal() +
  theme(legend.position = "none")

print(ct3)

ggsave("./figs/ct_spec.png", width = 12, height = 6)
ggsave("./figs/ct_spec.pdf", width = 12, height = 6)

hit_plot <- data.frame(
  Hit = c(unlist(all_hits), unlist(all_hits_o2)),
  Parameters = rep(rep(params, each = length(newsims)), 2),
  Model = c(rep("This Paper", length(newsims) * length(xis) * length(nus)), rep("outbreaker2 (2018)", length(newsims) * length(xis) * length(nus)))
)
hit_plot$Parameters <- factor(hit_plot$Parameters, levels = params)
hit_plot$Model <- factor(hit_plot$Model, levels = c("This Paper", "outbreaker2 (2018)"))


ct4 <- ggplot(hit_plot, aes(x=Parameters, y=Hit, fill=Model, color = Model)) +
  geom_boxplot() +
  scale_fill_manual(values=c( "#EECCAA", "#AACCEE")) +
  scale_color_manual(values=c("#BB9977", "#7799BB")) +
  scale_x_discrete(labels = addline_format(params)) +
  ylab("Hit Rate") +
  theme_minimal() + 
  theme(legend.position = "none")

print(ct4)

ggsave("./figs/ct_hit_rate.png", width = 12, height = 6)
ggsave("./figs/ct_hit_rate.pdf", width = 12, height = 6)




