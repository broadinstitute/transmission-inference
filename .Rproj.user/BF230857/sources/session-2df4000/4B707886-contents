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

### Simulate outbreak

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

set.seed(211)

### Generate transmission network

min_n_cluster = 10
max_n_cluster = 15
min_n_unrelated = 0
max_n_unrelated = 5
window <- 10 # max number of days between in-group and out-group case

# Epi params
tau_T <- 2
tau_E <- 2
tau_I <- 6

var_T <- 2
var_E <- 2
var_I <- 6

# Genetics params
N_bases <- 29903 # length of viral genome
p <- 3.38e-6 # probability of mutation per site per RNA replication
w <- 1000 # virions produced per cycle
mu <- 2.74e-6 # probability of substitution per site per unit time
k <- 1/sqrt(p) # number of virions at end of exponential growth phase
lambda <- (mu/p) - mu

# Probability of a mutation in growth phase
p_growth_mut <- 1 - (1-p)^k

# Number of generations to simulate
n_gen <- 4

# Generate frequency of iSNV per site, from 100k dataset
context <- read.csv(system.file("extdata", "context.csv", package = "reconstructR"))
tot_isnvs <- data.frame(MUT = context$aa_change_full, POS = context$POS, COUNT = Rfast::rowsums(as.matrix(context[,4:ncol(context)]), na.rm = T))
# Probability per position
N_context <- 172519
pos_probs <- c()
hgeom_probs <- c()
for (i in 1:N_bases) {
  sub <- tot_isnvs[tot_isnvs$POS == i, ]
  pos_probs[i] <- sum(sub$COUNT)
  hgeom_probs[i] <- 1 - phyper(1, pos_probs[i], N_context - pos_probs[i], 1000)
}
pos_probs <- pos_probs / sum(pos_probs)



run_sim <- function(lol){
  
  setwd("~/Desktop/sim-outbreak/")
  
  count <- 1
  subtree_sizes = 0
  
  while (count < max_n_cluster + max_n_unrelated | all(subtree_sizes < min_n_cluster | subtree_sizes > max_n_cluster)) {
    # Ancestor
    anc <- NA
    
    # Bottleneck size = 1
    A = 1
    
    # Set exposure time to 0 for first case
    t_E <- 0
    
    # Time to reach k virions
    t_G <- t_E + (1/(lambda * log(w)))*log(k)
    
    # Sample infectious time
    t_I <- t_G + rgamma(1, tau_E^2/var_E, tau_E/var_E)
    
    # Sample test time
    t_test <- t_I + rgamma(1, tau_T^2/var_T, tau_T/var_T)
    
    # Ancestral heritage
    heritage <- list()
    heritage[[1]] <- integer(0)
    
    trans <- data.frame(id = 1, anc = NA, A = A, t_E = t_E, t_G = t_G, t_I = t_I, t_test = t_test)
    count <- 1 # what row are we on?
    gen <- 1 # what row ends the generation?
    
    for(g in 1:n_gen){
      if(count <= gen){
        for (i in count:gen) {
          n_offspring <- rnbinom(1, 0.1, 0.1/2.6)
          if(n_offspring > 0){
            A <- 1 + rpois(n_offspring, 0.5)
            t_E <- trans$t_I[i] + rexp(n_offspring, 1/tau_I)
            t_G <- t_E + (1/(lambda * log(w)))*log(k/A)
            
            # Sample infectious time
            t_I <- t_G + rgamma(n_offspring, tau_E^2/var_E, tau_E/var_E)
            
            # Sample test time
            t_test <- t_I + rgamma(n_offspring, tau_T^2/var_T, tau_T/var_T)
            
            for (j in (max(trans[,1]) + 1):(max(trans[,1]) + n_offspring)) {
              heritage[[j]] <- c(heritage[[i]], i)
            }
            
            trans <- rbind(trans, data.frame(
              id = (max(trans[,1]) + 1):(max(trans[,1]) + n_offspring),
              anc = trans$id[i],
              A = A,
              t_E = t_E,
              t_G = t_G,
              t_I = t_I,
              t_test = t_test
            ))
            
          }
          
          count <- count + 1
          gen <- gen + n_offspring
        }
      }
    }
    
    #print(trans)
    
    subtree_sizes <- c()
    for (i in trans$id) {
      subtree_sizes[i] <- sum(unlist(heritage) == i) + 1 # adding one for root
    }
    
    #print(subtree_sizes)
  }
  
  
  # Proportions of mutated particles at end of growth phase
  x <- list()
  
  # Which sites have a mutation in the exponential growth phase?
  n_growth_mut <- rbinom(1, N_bases, p_growth_mut)
  which_growth_mut <- sample(1:N_bases, n_growth_mut, prob = pos_probs)
  
  # If so, what fraction of particles have said mutation at end of expo growth phase?
  x[[1]] <- rep(0, N_bases)
  x[[1]][which_growth_mut] <- rbeta(n_growth_mut, 1, sample(1:k, n_growth_mut, replace = T))
  
  # Depth
  depths <- list()
  depths[[1]] <- sample(100:20000, N_bases, replace = T)
  
  # Reads of the alternate allele (C)
  reads <- list()
  reads[[1]] <- rbinom(N_bases, depths[[1]], (x[[1]] * (1-mu*(trans$t_test[1] - trans$t_G[1])) +  (1 - x[[1]])*mu*(trans$t_test[1] - trans$t_G[1])))
  
  # Consensus genomes
  cons <- list()
  cons[[1]] <- reads[[1]] > (depths[[1]] - reads[[1]])
  cons[[1]] <- c("A", "C")[cons[[1]] + 1]
  
  for (i in 2:nrow(trans)) {
    
    # Bottleneck
    B <- sapply(1:trans$A[i], function(r){runif(N_bases) < (x[[trans$anc[i]]] * (1 - mu*(trans$t_E[i] - trans$t_G[trans$anc[i]])) + (1 - x[[trans$anc[i]]]) * mu*(trans$t_E[i] - trans$t_G[trans$anc[i]]))})
    
    # Which sites have a mutation in the exponential growth phase?
    n_growth_mut <- rbinom(1, N_bases, p_growth_mut)
    which_growth_mut <- sample(1:N_bases, n_growth_mut, prob = pos_probs)
    
    x[[i]] <- rep(0, N_bases)
    x[[i]][which_growth_mut] <- rbeta(n_growth_mut, 1, sample(1:k, n_growth_mut, replace = T))
    
    x[[i]][rowSums(B) == trans$A[i]] <- 1 - x[[i]][rowSums(B) == trans$A[i]] # If all bottleneck particles were mutated, flip the fraction of mutated particles
    
    # Finally, we need a correction if the bottleneck itself was split...
    split <- (rowSums(B) > 0) & (rowSums(B) < trans$A[i])
    #print(sum(split))
    x[[i]][split] <- rbeta(sum(split), rowSums(B)[split], trans$A[i] - rowSums(B)[split])
    
    depths[[i]] <- sample(100:20000, N_bases, replace = T)
    reads[[i]] <- rbinom(N_bases, depths[[i]], (x[[i]] * (1-mu*(trans$t_test[i] - trans$t_G[i])) +  (1 - x[[i]])*mu*(trans$t_test[i] - trans$t_G[i])))
    
    cons[[i]] <- reads[[i]] > (depths[[i]] - reads[[i]])
    cons[[i]] <- c("A", "C")[cons[[i]] + 1]
    
  }
  
  
  
  ### Now, extract the cluster we want to work with
  viable_roots <- which(subtree_sizes >= min_n_cluster & subtree_sizes <= max_n_cluster)
  if(length(viable_roots) == 1){
    root <- viable_roots
  }else{
    root <- sample(viable_roots, 1)
  }
  who <- root
  for (i in trans$id) {
    if(root %in% heritage[[i]]){
      who <- c(who, i)
    }
  }
  cluster <- trans[who, ]
  
  unrelated <- which(trans$t_test > min(cluster$t_test) - 10 & trans$t_test < max(cluster$t_test) + 10)
  unrelated <- setdiff(unrelated, c(who, trans$anc[root]))
  
  # If unrelated has more than max_n_unrelated people, sample who it is. Otherwise keep unrelated as it is
  if(length(unrelated) > max_n_unrelated){
    unrelated <- sample(unrelated, sample(min_n_unrelated:max_n_unrelated, 1))
  }
  
  cluster <- rbind(cluster, trans[unrelated, ])
  who <- c(who, unrelated)
  
  cluster$id <- match(cluster$id, who)
  cluster$anc <- match(cluster$anc, who)
  rownames(cluster) <- NULL
  
  # Filter genetic info
  x <- x[who]
  depths <- depths[who]
  reads <- reads[who]
  cons <- cons[who]
  
  # Number of people
  N <- length(who)
  
  names <- paste(1:N)
  names <- str_pad(names, 3, pad = "0")
  names <- paste0("sim_", names)
  
  # For writing files, move into new directory
  dirname <- paste0("newdir_", lol)
  dir.create(paste0("./", dirname), showWarnings = F)
  setwd(paste0("./", dirname))
  dir.create("./input_data", showWarnings = F)
  dir.create("./input_data/vcf", showWarnings = F)
  
  # Write fastas
  names(cons) <- names
  write.dna(cons, file = "./input_data/aligned.fasta", format = "fasta")
  
  # Write ref genome
  ref <- list(rep("A", N_bases))
  names(ref) <- "ref"
  write.dna(ref, file = "./input_data/ref.fasta", format = "fasta")
  
  # Write read depth table
  depth_table <- matrix(unlist(depths), ncol = N)
  colnames(depth_table) <- names(cons)
  position <- 1:N_bases
  depth_table <- cbind(position, depth_table)
  write.csv(depth_table, "./input_data/depth.csv", row.names = F)
  
  # Write test date table
  dates <- cbind(names(cons), round(cluster$t_test - min(cluster$t_test)))
  write.csv(dates, "./input_data/date.csv", row.names = F)
  
  # Write VCFs
  
  # First, remove existing vcfs, just in case
  do.call(file.remove, list(list.files("./input_data/vcf/", full.names = TRUE)))
  
  for (i in 1:N) {
    # Which positions have iSNVs in read data?
    pos <- which(reads[[i]] > 0 & reads[[i]] < depths[[i]])
    info <- paste0(
      "DP=",
      depths[[i]][pos],
      ";AF=",
      format(round(reads[[i]][pos] / depths[[i]][pos], 6), scientific = F),
      ";SB=0;DP4=0,", # strand bias wouldn't actually be 0 here, but assuming no strand bias in synthetic data
      depths[[i]][pos] - reads[[i]][pos],
      ",0,",
      reads[[i]][pos]
    )
    vcf <- data.frame(
      CHROM = "ref",
      POS = pos,
      ID = ".",
      REF = "A",
      ALT = "C",
      QUAL = 5000,
      FILTER = "PASS",
      INFO = info
    )
    
    write.table(vcf, file = paste0("./input_data/vcf/", names[i], ".vcf"), col.names = F, row.names = F)
  }
  
  truth <- cbind(cluster$anc, 1:N) + 1
  truth[is.na(truth)] <- 1
  colnames(truth) <- NULL
  
  dir.create("./badtrip", showWarnings = F)
  
  # First, remove existing files, just in case
  do.call(file.remove, list(list.files("./badtrip/", full.names = TRUE)))
  
  
  gen_input <- matrix(ncol = N, nrow = N_bases)
  diverse <- c()
  which_ref <- sample(1:4, N_bases, replace = T)
  which_alt <- ((which_ref - 1 + sample(1:3, N_bases, replace = T)) %% 4) + 1
  for (i in 1:N) {
    ref <- depths[[i]] - reads[[i]]
    alt <- reads[[i]]
    alt[(alt / depths[[i]]) < 0.03] <- 0
    ref[(ref / depths[[i]]) < 0.03] <- 0
    diverse <- union(diverse, which(alt > 0))
    
    out <- matrix(0, nrow = N_bases, ncol = 4)
    
    out[cbind(1:N_bases, which_ref)] <- ref
    out[cbind(1:N_bases, which_alt)] <- alt
    gen_input[,i] <- apply(out, 1, function(x){paste(x, collapse = "-")})
  }
  #gen_input <- gen_input[diverse, ]
  colnames(gen_input) <- paste0("S", 1:N)
  
  write.table(gen_input, file = "./badtrip/inputAlignment.txt", sep = "   ", quote = F, row.names = F)
  
  # Badtrip sample input
  sample_input <- data.frame(paste0("H", 1:N), paste0("S", 1:N), round(cluster$t_test))
  write.table(sample_input, file = "./badtrip/inputSample.txt", sep = "   ", quote = F, col.names = F, row.names = F)
  
  # Badtrip
  epi_input <- data.frame(paste0("H", 1:N), round(cluster$t_test) - 10, round(cluster$t_test) + 10)
  write.table(epi_input, file = "./badtrip/inputEpi.txt", sep = "   ", quote = F, col.names = F, row.names = F)
  
  ### Now, let's reconstruct it!
  
  setwd("~/Desktop/sim-outbreak/")
  setwd(paste0("./", dirname))
  
  results <- run_mcmc(
    N_iters = 1e4
  )
  burnin <- 0.1
  
  N <- results[[1]]$N
  adj <- matrix(0, nrow = N, ncol = N)
  for (i in 1:length(results)) {
    if (i > burnin * length(results)) {
      inds <- cbind(results[[i]]$anc[2:N], 2:N)
      adj[inds] <- adj[inds] + 1
    }
  }
  
  adj <- adj / sum(1:length(results) > burnin * length(results))
  
  return(list(truth, adj))
}

#newsims[failed] <- mclapply(failed, run_sim, mc.silent = T, mc.cores = 12)
#save(newsims, file = "~/Desktop/sim-outbreak/newsims2.RData")




#####-------Analysis
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

## Get hit rate
get_hit_rate <- function(truth, adj){
  unrelated <- which(truth[,1] == 1)
  related <- setdiff(1:nrow(truth), unrelated)
  
  bin_adj <- matrix(0, nrow = nrow(adj), ncol = ncol(adj))
  bin_adj[cbind(colMaxs(adj)[2:ncol(adj)], 2:ncol(adj))] <- 1
  
  related_hit <- bin_adj[truth[related, ]]
  mean(related_hit)
}


means <- c()
truths <- list()
adjs <- list()
specs <- c()
hit_rates <- c()
for (i in 1:100) {
  truth <- newsims[[i]][[1]]
  adj <- newsims[[i]][[2]]
  adjs[[i]] <- adj
  truths[[i]] <- truth
  
  means[i] <- get_mean(truth, adj)
  specs[i] <- get_mean(truth, adj, spec = T)
  hit_rates[i] <- get_hit_rate(truth, adj)
}


## Badtrip comparison
run_badtrip <- function(lol){
  
  setwd("~/Desktop/sim-outbreak/")
  dirname <- paste0("newdir_", lol)
  setwd(paste0("./", dirname))
  setwd("./badtrip")
  
  file.remove("BADTRIP_setup.log")
  file.remove("BADTRIP_setup.trees")
  file.remove("BADTRIP_setup.xml.state")
  file.remove("summary_network.txt")
  file.remove("BADTRIP_setup.xml")
  
  
  system(
    "python3 ~/Desktop/sim-outbreak/create_BADTRIP_xml.py -a inputAlignment.txt -e inputEpi.txt -s inputSample.txt -m 100000 -o BADTRIP_setup"
  )
  
  # Run badtrip
  
  system(
    "/Applications/BEAST_2.7.5/bin/beast -threads 8 BADTRIP_setup.xml"
  )
  
  # Process tree into summary file
  system(
    "python3 ~/Desktop/sim-outbreak/Make_transmission_tree_alternative.py -i BADTRIP_setup.trees -b 10 -o summary"
  )
  
}

#
# failed <- c()
# for(lol in 1:100){
#   setwd("~/Desktop/sim-outbreak/")
#   dirname <- paste0("newdir_", lol)
#   setwd(paste0("./", dirname))
#   setwd("./badtrip")
#
#   if(!("summary_network.txt" %in% list.files())){
#     failed <- c(failed, lol)
#   }
# }
#
#
## TO RUN BADTRIP
# for (lol in failed) {
#   print(lol)
#   run_badtrip(lol)
# }
#
#mclapply(failed, run_badtrip, mc.cores = 12)




### Outbreaker comparison
compare_outbreaker <- function(sims, lol){
  
  setwd(paste0("~/Desktop/sim-outbreak/newdir_", lol))
  
  dna <- read.FASTA("./input_data/aligned.fasta")
  names(dna) <- NULL
  dat <- list(
    dna = dna,
    dates = read.csv("./input_data/date.csv")[,2],
    w_dens = dgamma(1:20, 3)
  )
  out <- outbreaker(
    data = dat,
    config = list(n_iter = 1e4, sample_every = 100)
  )
  
  res <- as.matrix(out)
  res <- res[, 9:(9+length(dna)-1)]
  res <- res + 1
  res[is.na(res)] <- 1
  adj <- matrix(0, nrow = length(dna)+1, ncol = length(dna)+1)
  burnin <- 0.1
  
  for (i in 1:nrow(res)) {
    if (i > burnin * nrow(res)) {
      inds <- cbind(res[i, ], 2:(length(dna)+1))
      adj[inds] <- adj[inds] + 1
    }
  }
  
  adj <- adj / sum(1:nrow(res) > burnin * nrow(res))
  
  return(
    list(
      get_mean(sims[[lol]][[1]], adj),
      adj
    )
  )
}

## our method with consensus genomes
compare_cons <- function(sims, lol){
  
  dirname <- paste0("newdir_", lol)
  
  setwd("~/Desktop/sim-outbreak/")
  setwd(paste0("./", dirname))
  
  results <- run_mcmc(
    N_iters = 1e4,
    filters = list(
      omit = character(0), # any samples we should omit?
      coverage = 0.75, # minimum coverage
      sb = Inf, # filter for strand bias (Inf if not using)
      dp4 = 0, # p-value cutoff for fishers exact test on ref-fwd, ref-rev, alt-fwd, alt-rev
      af = 0.5, # filter for minor allele frequency
      call = 10, # minimum reads of the alternate allele to call iSNV
      depth = 100, # minimum total depth
      problem_qtile = 0.05
    )
  )
  burnin <- 0.1
  
  N <- results[[1]]$N
  adj <- matrix(0, nrow = N, ncol = N)
  for (i in 1:length(results)) {
    if (i > burnin * length(results)) {
      inds <- cbind(results[[i]]$anc[2:N], 2:N)
      adj[inds] <- adj[inds] + 1
    }
  }
  
  adj <- adj / sum(1:length(results) > burnin * length(results))
  
  return(
    list(
      get_mean(sims[[lol]][[1]], adj),
      adj
    )
  )
}

## Compare to badtrip
compare_badtrip <- function(sims, lol){
  setwd(paste0("~/Desktop/sim-outbreak/newdir_", lol, "/badtrip"))
  
  net <- readLines("summary_network.txt")
  # Get number of hosts
  N <- length(unlist(strsplit(net[1], ","))) - 1
  start <- which(net == "Probabilities of direct transmittor to each sampled host: ")
  
  adj <- matrix(0, nrow = N+1, ncol = N+1)
  for (r in (start+1):length(net)) {
    if(grepl("To host", net[r])){
      to <- sub("To host HH", "", net[r])
      to <- sub(" from : ", "", to)
      to <- as.numeric(to)
      
      from <- unlist(strsplit(net[r+1], ", "))
      probs <- sub(".* ", "", from)
      probs <- as.numeric(probs)
      
      from <- sub(" .*", "", from)
      from <- sub("HH", "", from)
      from[from == "Unsampled"] <- "0" # then adding 1
      from <- as.numeric(from)
      
      adj[from + 1, to + 1] <- probs
      
    }
  }
  
  return(
    list(
      get_mean(sims[[lol]][[1]], adj),
      adj
    )
  )
  
}


control_means <- mclapply(1:100, compare_outbreaker, sims = newsims, mc.cores = 12)
badtrip_means <- lapply(1:100, compare_badtrip, sims = newsims)

cons_means <- mclapply(1:100, compare_cons, sims = newsims, mc.cores = 12)


outbreaker_adj <- list()
badtrip_adj <- list()
cons_adj <- list()
for (i in 1:100) {
  outbreaker_adj[[i]] <- control_means[[i]][[2]]
  control_means[[i]] <- control_means[[i]][[1]]
  
  badtrip_adj[[i]] <- badtrip_means[[i]][[2]]
  badtrip_means[[i]] <- badtrip_means[[i]][[1]]
  
  cons_adj[[i]] <- cons_means[[i]][[2]]
  cons_means[[i]] <- cons_means[[i]][[1]]
}

## get accuracy and specificity
control_spec <- c()
badtrip_spec <- c()
cons_spec <- c()
for (i in 1:100) {
  control_means[i] <- get_mean(newsims[[i]][[1]], outbreaker_adj[[i]])
  badtrip_means[i] <- get_mean(newsims[[i]][[1]], badtrip_adj[[i]])
  cons_means[i] <- get_mean(newsims[[i]][[1]], cons_adj[[i]])
  
  control_spec[i] <- get_mean(newsims[[i]][[1]], outbreaker_adj[[i]], spec = T)
  badtrip_spec[i] <- get_mean(newsims[[i]][[1]], badtrip_adj[[i]], spec = T)
  cons_spec[i] <- get_mean(newsims[[i]][[1]], cons_adj[[i]], spec = T)
  
}



control_means <- unlist(control_means)
badtrip_means <- unlist(badtrip_means)
cons_means <- unlist(cons_means)


mean(means)
sd(means)

mean(control_means)
sd(control_means)

mean(badtrip_means)
sd(badtrip_means)

mean(cons_means)
sd(cons_means)

t.test(means, control_means)


mean(specs)
sd(specs)
mean(control_spec)
sd(control_spec)
mean(badtrip_spec)
sd(badtrip_spec)
mean(cons_spec)
sd(cons_spec)


### Compute the coverage probability for the new method
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

covs <- mapply(get_cov, truths, adjs)
outbreaker_covs <- mapply(get_cov, truths, outbreaker_adj)
cons_covs <- mapply(get_cov, truths, cons_adj)
badtrip_covs <- mapply(get_cov, truths, badtrip_adj)


mean(covs)
sd(covs)

mean(outbreaker_covs)
sd(outbreaker_covs)

mean(cons_covs)
sd(cons_covs)

mean(badtrip_covs)
sd(badtrip_covs)

wilcox.test(covs, cons_covs)

## Hit rate analysis
outbreaker_hit_rates <- c()
cons_hit_rates <- c()
badtrip_hit_rates <- c()

for (i in 1:100) {
  outbreaker_hit_rates[i] <- get_hit_rate(truths[[i]], outbreaker_adj[[i]])
  cons_hit_rates[i] <- get_hit_rate(truths[[i]], cons_adj[[i]])
  badtrip_hit_rates[i] <- get_hit_rate(truths[[i]], badtrip_adj[[i]])
}

mean(hit_rates)
sd(hit_rates)

mean(badtrip_hit_rates)
sd(badtrip_hit_rates)

mean(outbreaker_hit_rates)
sd(outbreaker_hit_rates)

mean(cons_hit_rates)
sd(cons_hit_rates)

### Box plots of accuracy and coverage

setwd("~/Desktop/outbreak-reconstruction/simulated-outbreaks//")

Model <-  c(
  rep("This Paper", 100),
  rep("This Paper, Consensus", 100),
  rep("BadTrIP (2018)", 100),
  rep("outbreaker2 (2018)", 100)
)
Model <- factor(Model, levels = c("This Paper", "BadTrIP (2018)", "outbreaker2 (2018)", "This Paper, Consensus"))

box1 <- data.frame(
  Accuracy = c(means, cons_means, badtrip_means, control_means),
  Model = Model
)

b1 <- ggplot(box1, aes(x = Model, y = Accuracy, fill = Model, color = Model)) +
  geom_boxplot(outlier.size = -1) +
  geom_jitter(shape=16, position=position_jitter(0.2)) +
  scale_fill_manual(values=c('#EECCAA', "#AAEECC", "#AACCEE", "#EEAACC")) +
  scale_color_manual(values=c('#BB9977', "#77BB99", "#7799BB", "#BB7799")) +
  xlab("Model") +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  ) +
  theme(legend.position = "none")

ggsave("./figs/acc.pdf", width = 6, height = 6)
ggsave("./figs/acc.png", width = 6, height = 6)

box2 <- data.frame(
  Coverage = c(covs, cons_covs, badtrip_covs, outbreaker_covs),
  Model = Model
)

b2 <- ggplot(box2, aes(x = Model, y = Coverage, fill = Model, color = Model)) +
  geom_boxplot(outlier.size = -1) +
  geom_jitter(shape=16, position=position_jitter(0.2)) +
  scale_fill_manual(values=c('#EECCAA', "#AAEECC", "#AACCEE", "#EEAACC")) +
  scale_color_manual(values=c('#BB9977', "#77BB99", "#7799BB", "#BB7799")) +
  xlab("Model") +
  ylab("Coverage Probability") +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  ) +
  theme(legend.position = "none")

ggsave("./figs/cov.pdf", width = 6, height = 6)
ggsave("./figs/cov.png", width = 6, height = 6)

## Box plot of specificity

box3 <- data.frame(
  Specificity = c(specs, cons_spec, badtrip_spec, control_spec),
  Model = Model
)

b3 <- ggplot(box3, aes(x = Model, y = Specificity, fill = Model, color = Model)) +
  geom_boxplot(outlier.size = -1) +
  geom_jitter(shape=16, position=position_jitter(0.2)) +
  scale_fill_manual(values=c('#EECCAA', "#AAEECC", "#AACCEE", "#EEAACC")) +
  scale_color_manual(values=c('#BB9977', "#77BB99", "#7799BB", "#BB7799")) +
  xlab("Model") +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  ) +
  theme(legend.position = "none")

ggsave("./figs/spec.pdf", width = 6, height = 6)
ggsave("./figs/spec.png", width = 6, height = 6)

box4 <- data.frame(
  Hit = c(hit_rates, cons_hit_rates, badtrip_hit_rates, outbreaker_hit_rates),
  Model = Model
)

b4 <- ggplot(box4, aes(x = Model, y = Hit, fill = Model, color = Model)) +
  geom_boxplot(outlier.size = -1) +
  geom_jitter(shape=16, position=position_jitter(0.2)) +
  scale_fill_manual(values=c('#EECCAA', "#AAEECC", "#AACCEE", "#EEAACC")) +
  scale_color_manual(values=c('#BB9977', "#77BB99", "#7799BB", "#BB7799")) +
  xlab("Model") +
  ylab("Hit Rate") +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  ) +
  theme(legend.position = "top")

legend <- as_ggplot(ggpubr::get_legend(b4))

b4 <- b4 + theme(legend.position = "none")

ggsave("~/Desktop/outbreak-reconstruction/simulated-outbreaks/figs/hit_rate.pdf", width = 6, height = 6)
ggsave("~/Desktop/outbreak-reconstruction/simulated-outbreaks/figs/hit_rate.png", width = 6, height = 6)


plot_grid(
  plot_grid(legend, labels = "", ncol = 1),
  plot_grid(b1, b2, b3, b4, labels = c(LETTERS[1:4], ""), ncol = 4),
  plot_grid(ct1, ct2, ct3, ct4, labels = c(LETTERS[5:8], ""), ncol = 2),
  ncol = 1,
  rel_heights = c(0.2, 1, 1)
)

ggsave("~/Desktop/outbreak-reconstruction/simulated-outbreaks/figs/fig3.pdf", width = 12, height = 8)
ggsave("~/Desktop/outbreak-reconstruction/simulated-outbreaks/figs/fig3.png", width = 12, height = 8)




