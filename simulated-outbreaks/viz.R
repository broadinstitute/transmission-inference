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

## Visualize the epidemiology of a single simulated outbreak

# Run a simulated outbreak

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
library(cowplot)

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
dirname <- "demo"
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

### Compile relevant epidemiological information

dir.create("./figs", showWarnings = F)

N <- N-1

padding <- 0.06
linewidth <- 0.1

e1 <- ggplot() +
  geom_rect(aes(xmin = cluster$t_E - min(cluster$t_E), xmax = cluster$t_I - min(cluster$t_E), ymin = 1:N - 1/2 + padding/2, ymax = 1:N + 1/2 - padding/2, fill = "Exposed Stage")) +
  geom_rect(aes(xmin = cluster$t_I - min(cluster$t_E), xmax = cluster$t_I + 13 - min(cluster$t_E), ymin = 1:N - 1/2 + padding/2, ymax = 1:N + 1/2 - padding/2, fill = "Infectious Stage")) +
  geom_rect(aes(xmin = cluster$t_E[2:N] - min(cluster$t_E), xmax = cluster$t_E[2:N] + linewidth - min(cluster$t_E), ymin = pmin(2:N, cluster$anc[2:N]), ymax = pmax(2:N, cluster$anc[2:N]) + 1/2 - padding/2), fill = 'black') +
  geom_point(aes(x = cluster$t_E[2:N] + linewidth / 2 - min(cluster$t_E), y = pmin(2:N, cluster$anc[2:N]), color = 'Transmission'), size = 10*linewidth) +
  scale_y_continuous(breaks = 1:N, minor_breaks = NULL) +
  scale_fill_manual(
    name = "",
    values =c("Exposed Stage"= '#DDBB99',"Infectious Stage"='#DD9999')
  ) +
  scale_color_manual(
    name = "",
    values = 'black'
  ) +
  xlab("Time (days)") +
  ylab("Case ID") +
  labs(color=NULL) +
  theme_minimal() +
  theme(legend.position = "top")

print(e1)


ggsave("./figs/epi.pdf", width = 8, height = 6)
ggsave("./figs/epi.png", width = 8, height = 6)

### Genomics of outbreak
# Get all sites of interest (iSNV/SNV in at least one case)
read_mtx <- matrix(nrow = N_bases, ncol = 0)
for (i in 1:N) {
  newcol <- reads[[i]]/depths[[i]]
  newcol[newcol < 0.03] <- 0
  newcol[newcol > 0.97] <- 1
  read_mtx <- cbind(read_mtx, newcol)
}
sites <- which(
  rowsums(read_mtx > 0) >= 2
)
snvs <- paste0(sample(c("A", "C", "G", "T"), length(sites), replace = T), sites, sample(c("A", "C", "G", "T"), length(sites), replace = T))

gen_data <- data.frame()
for (i in 1:N) {
  newdata <- data.frame(site = snvs, Frequency = reads[[i]][sites]/depths[[i]][sites], id = i)
  newdata$Frequency[reads[[i]][sites] < 10 | newdata$Frequency < 0.03] <- 0
  newdata$Frequency[depths[[i]][sites] - reads[[i]][sites] < 10 | newdata$Frequency > 0.97] <- 1
  gen_data <- rbind(
    gen_data,
    newdata
  )
}

gen_data$id <- as.factor(gen_data$id)
gen_data$site <- factor(gen_data$site, levels = snvs)

gen_data <- gen_data[gen_data$Frequency > 0, ]


e2 <- ggplot(gen_data) +
  geom_tile(aes(x = site, y = id, alpha = Frequency), fill = '#BB5522', color = 'white', linewidth = 1) +
  scale_alpha_continuous(range = c(0.25,1)) +
  scale_y_discrete(limits = paste(1:N)) +
  xlab("Substitution") +
  ylab("Case ID") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), legend.position = "top")

print(e2)

ggsave("./figs/genomics.pdf", width = 4, height = 6)
ggsave("./figs/genomics.png", width = 4, height = 6)


unrelated <- truth[which(truth[,1] == 1), 2]
### Next figures: reconstructed transmission network with correct links shown, outbreaker2 vs. new method
N <- N+1
inferred <- data.frame()
for (i in 1:N) {
  for (j in 1:N) {
    newrow <- data.frame(from = i, to = j, prob = adj[i,j])
    if(sum(truth[,1] == i & truth[,2] == j) > 0){
      newrow$truth = 1
    }else{
      newrow$truth = 0
    }

    if(newrow$truth == 1 | newrow$prob > 0){
      inferred <- rbind(inferred, newrow)
    }
  }
}

inferred$prob[inferred$to %in% unrelated & inferred$from %in% unrelated] <- 0
inferred <- inferred[inferred$from != 1, ]
inferred$from <- paste(inferred$from)
inferred$to <- paste(inferred$to)


g <- graph_from_edgelist(as.matrix(inferred[, 1:2]))


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
adj_o2 <- matrix(0, nrow = length(dna)+1, ncol = length(dna)+1)
burnin <- 0.1

for (i in 1:nrow(res)) {
  if (i > burnin * nrow(res)) {
    inds <- cbind(res[i, ], 2:(length(dna)+1))
    adj_o2[inds] <- adj_o2[inds] + 1
  }
}

adj_o2 <- adj_o2 / sum(1:nrow(res) > burnin * nrow(res))


inferred_o2 <- data.frame()
for (i in 1:N) {
  for (j in 1:N) {
    newrow <- data.frame(from = i, to = j, prob = adj_o2[i,j])
    if(sum(truth[,1] == i & truth[,2] == j) > 0){
      newrow$truth = 1
    }else{
      newrow$truth = 0
    }

    if(newrow$truth == 1 | newrow$prob > 0){
      inferred_o2 <- rbind(inferred_o2, newrow)
    }
  }
}

inferred_o2$prob[inferred_o2$to %in% unrelated & inferred_o2$from %in% unrelated] <- 0
inferred_o2 <- inferred_o2[inferred_o2$from != 1, ]
inferred_o2$from <- paste(inferred_o2$from)
inferred_o2$to <- paste(inferred_o2$to)

inferred$truth <- factor(inferred$truth)
inferred_o2$truth <- factor(inferred_o2$truth)


coords <- layout_nicely(graph_from_data_frame(inferred))

layout <- create_layout(graph_from_data_frame(inferred), coords)
e3 <- ggraph(layout) + 
  geom_edge_link(aes(color = truth, alpha = prob), width = 1, arrow = arrow(type = "closed", length = unit(3, 'mm')), 
                 start_cap = circle(5, 'mm'),
                 end_cap = circle(5, 'mm')) + 
  geom_node_point(size = 6, color = "#555555") +
  geom_node_text(aes(label = as.numeric(name) - 1), color = "white", size = 4) + 
  scale_edge_color_manual(values= c("#BB5522", '#55BB22')) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_blank(), legend.position = "none")
print(e3)

new_layout <- create_layout(graph_from_data_frame(inferred_o2), coords)
new_layout <- create_layout(graph_from_data_frame(inferred_o2), coords[match(new_layout$name, layout$name), ])

e4 <- ggraph(new_layout) + 
  geom_edge_link(aes(color = truth, alpha = prob), width = 1, arrow = arrow(type = "closed", length = unit(3, 'mm')), 
                 start_cap = circle(5, 'mm'),
                 end_cap = circle(5, 'mm')) + 
  geom_node_point(size = 6, color = "#555555") +
  geom_node_text(aes(label = as.numeric(name) - 1), color = "white", size = 4) + 
  scale_edge_color_manual(values= c("#BB5522", '#55BB22')) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_blank(), legend.position = "none")
print(e4)


figs1 <- plot_grid(e1,e2,e3,e4, labels = LETTERS[1:4], ncol = 2)
figs1
ggsave("./figs/figs1.pdf", width = 10, height = 8)
ggsave("./figs/figs1.png", width = 10, height = 8)



