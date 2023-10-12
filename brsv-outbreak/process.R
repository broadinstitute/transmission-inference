### This script runs the four outbreak reconstruction algorithms on the BRSV outbreak.

# Relevant libraries
library(ape)
library(reconstructR)
library(outbreaker2)
library(igraph)
library(cowplot)
library(ggraph)

# Set the working directory
setwd("./brsv-outbreak/")

# Read in the consensus genomes
cons <- read.FASTA("./raw/expt_G_all_consensus.fasta")

# Raw read-level data files
filepaths <- list.files("./raw/vcf/")

# Number of cows
N <- length(filepaths)

# Length of BRSV genome
N_bases <- 14562

# Initialize read depth table
depth_table <- 1:N_bases

# Initialize genetic input table for BadTrIP
gen_input <- matrix(ncol = N, nrow = N_bases)

# Initialize vector of sites with any within-host divesrity
relevant <- c()

# Vector of nucleotides
ntides <- c("A", "C", "G", "T")

# Loop through read-level data files and process them into the correct format for outbreak reconstruction
for (i in 1:length(filepaths)) {
  
  # Read the data
  vcf <- read.csv(paste0("./raw/vcf/", filepaths[i]), sep = "\t")
  
  # Read depth at every site
  all_dp <- vcf$NonRefN + vcf$RefN

  # Update the read depth table
  depth_table <- cbind(depth_table, all_dp)

  # For Badtrip input data, apply 3% filter to all read counts
  relevant <- c(relevant, which(vcf$NonRefN / all_dp > 0.03))
  A_mod <- vcf$A
  A_mod[vcf$A / all_dp < 0.03] <- 0
  C_mod <- vcf$C
  C_mod[vcf$C / all_dp < 0.03] <- 0
  G_mod <- vcf$G
  G_mod[vcf$G / all_dp < 0.03] <- 0
  T_mod <- vcf$T
  T_mod[vcf$T / all_dp < 0.03] <- 0

  # Update BadTrIP input file
  gen_input[,i] <- paste(A_mod, C_mod, G_mod, T_mod, sep = "-")

  # For speed, eliminate rows with no iSNV
  vcf <- vcf[vcf$NonRefN > 0, ]

  # Calculate number of reads of alternate nucleotide and strand bias
  ref <- vcf$Ref
  pos <- vcf$Pos
  alt <- c()
  sb <- c()

  ref_fwd <- c()
  ref_rev <- c()
  alt_fwd <- c()
  alt_rev <- c()

  for (j in 1:nrow(vcf)) {
    counts <- c(vcf$A[j], vcf$C[j], vcf$G[j], vcf$T[j])
    which_ref <- match(vcf$Ref[j], ntides)
    counts[which_ref] <- 0
    which_alt <- which.max(counts) # Whichever nucleotide has the second-highest read count is the "alternate" allele. We don't account for multiallelic sites.
    alt[j] <- ntides[which_alt]

    # Find column for ref_fwd, ref_rev, alt_fwd, alt_rev
    ref_fwd_col <- paste0(ntides[which_ref], "fwd")
    ref_rev_col <- paste0(ntides[which_ref], "rev")
    alt_fwd_col <- paste0(ntides[which_alt], "fwd")
    alt_rev_col <- paste0(ntides[which_alt], "rev")

    ref_fwd[j] <- vcf[j, ref_fwd_col]
    ref_rev[j] <- vcf[j, ref_rev_col]
    alt_fwd[j] <- vcf[j, alt_fwd_col]
    alt_rev[j] <- vcf[j, alt_rev_col]

    # Phred scoring for SB
    sb[j] <- round(-10 * log10(fisher.test(matrix(c(ref_fwd[j], ref_rev[j], alt_fwd[j], alt_rev[j]), ncol = 2))$p))

  }

  # Depth at these sites (after masking multiallelic sites, in case it changes)
  dp <- ref_fwd + ref_rev + alt_fwd + alt_rev

  # Create "info" column in standard VCF
  info <- paste0(
    "DP=",
    dp,
    ";AF=",
    format(round((alt_fwd + alt_rev) / dp, 6), scientific = F),
    ";SB=",
    sb,
    ";DP4=",
    paste(
      ref_fwd, ref_rev, alt_fwd, alt_rev, sep = ","
    )
  )
  
  # Create vcf table
  vcf <- data.frame(
    CHROM = "ref",
    POS = pos,
    ID = ".",
    REF = ref,
    ALT = alt,
    QUAL = 5000,
    FILTER = "PASS",
    INFO = info
  )

  # Write VCF file
  write.table(vcf, file = paste0("./input_data/vcf/", names(cons)[i], ".vcf"), col.names = F, row.names = F)

}

# Rename columns of depth_table to match case IDs, and write
colnames(depth_table) <- c("position", names(cons))
write.csv(depth_table, "./input_data/depth.csv", row.names = F)

# Time at which each cow tested positive (extracted from paper)
t_test <- c(
  7,
  21,
  21,
  21,
  7,
  14,
  14,
  21,
  0
)

# Write table of test dates
dates <- cbind(names(cons), t_test)
write.csv(dates, "./input_data/date.csv", row.names = F)

#Badtrip input
colnames(gen_input) <- paste0("S", 1:N)
relevant <- sort(unique(relevant))
gen_input <- gen_input[relevant, ]

write.table(gen_input, file = "./badtrip/inputAlignment.txt", sep = "   ", quote = F, row.names = F)

# Badtrip sample input
sample_input <- data.frame(paste0("H", 1:N), paste0("S", 1:N), t_test)
write.table(sample_input, file = "./badtrip/inputSample.txt", sep = "   ", quote = F, col.names = F, row.names = F)

# Badtrip
epi_input <- data.frame(paste0("H", 1:N), t_test - 10, t_test + 10)
write.table(epi_input, file = "./badtrip/inputEpi.txt", sep = "   ", quote = F, col.names = F, row.names = F)


### Run reconstructR

results <- run_mcmc(
  1e4,
  filters = list(
    omit = character(0), # any samples we should omit?
    coverage = 0.75, # minimum coverage
    sb = Inf, # filter for strand bias (Inf if not using)
    dp4 = 0.05, # p-value cutoff for fishers exact test on ref-fwd, ref-rev, alt-fwd, alt-rev
    af = 0.03, # filter for minor allele frequency
    call = 10, # minimum reads of the alternate allele to call iSNV
    depth = 100, # minimum total depth
    problem_qtile = 0 # if an iSNV is one of the problem_qtile-most prevalent in the Broad dataset, throw it out
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

### Run outbreaker2

dna <- read.FASTA("./input_data/aligned.fasta")
names(dna) <- NULL

dat <- list(
  dna = dna,
  dates = read.csv("./input_data/date.csv")[,2],
  w_dens = dgamma(1:20, 3)
)
out <- outbreaker(
  data = dat,
  config = list(
    n_iter = 1e4,
    sample_every = 100
  )
)

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

### Run reconstructR with consensus genomes only

results_cont <- run_mcmc(
  1e4,
  filters = list(
    omit = character(0), # any samples we should omit?
    coverage = 0.75, # minimum coverage
    sb = Inf, # filter for strand bias (Inf if not using)
    dp4 = 0.05, # p-value cutoff for fishers exact test on ref-fwd, ref-rev, alt-fwd, alt-rev
    af = 0.5, # filter for minor allele frequency
    call = 10, # minimum reads of the alternate allele to call iSNV
    depth = 100, # minimum total depth
    problem_qtile = 0 # if an iSNV is one of the problem_qtile-most prevalent in the Broad dataset, throw it out
  )
)

burnin <- 0.1

N <- results_cont[[1]]$N
adj_cont <- matrix(0, nrow = N, ncol = N)
for (i in 1:length(results_cont)) {
  if (i > burnin * length(results_cont)) {
    inds <- cbind(results_cont[[i]]$anc[2:N], 2:N)
    adj_cont[inds] <- adj_cont[inds] + 1
  }
}

adj_cont <- adj_cont / sum(1:length(results_cont) > burnin * length(results_cont))


### Run Badtrip

setwd("./badtrip")
# 
# file.remove("BADTRIP_setup.log")
# file.remove("BADTRIP_setup.trees")
# file.remove("BADTRIP_setup.xml.state")
# file.remove("summary_network.txt")
# file.remove("BADTRIP_setup.xml")
# 
# 
# system(
#   "python3 ~/Desktop/sim-outbreak/create_BADTRIP_xml.py -a inputAlignment.txt -e inputEpi.txt -s inputSample.txt -m 100000 -o BADTRIP_setup"
# )
# 
# # Run badtrip
# 
# system(
#   "/Applications/BEAST_2.7.5/bin/beast -threads 8 BADTRIP_setup.xml"
# )
# 
# # Process tree into summary file
# system(
#   "python3 ~/Desktop/sim-outbreak/Make_transmission_tree_alternative.py -i BADTRIP_setup.trees -b 10 -o summary"
# )

net <- readLines("summary_network.txt")
# Get number of hosts
N <- length(unlist(strsplit(net[1], ","))) - 1
start <- which(net == "Probabilities of direct transmittor to each sampled host: ")

adj_badtrip <- matrix(0, nrow = N+1, ncol = N+1)
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

    adj_badtrip[from + 1, to + 1] <- probs

  }
}

#### End Badtrip----------

### What is the true transmission network?
truth <- cbind(
  c(9,7,7,7,9,5,1,7,0),
  1:length(cons)
)

truth_names <- cbind(
  c(names(cons)[c(9,7,7,7,9,5,1,7,0)], NA),
  names(cons)
)

# Accuracy
mean(adj[truth + 1])
mean(adj_badtrip[truth + 1])
mean(adj_o2[truth + 1])
mean(adj_cont[truth + 1])


### Compute the coverage probability
get_cov <- function(truth, adj){
  
  out <- c()
  real <- truth
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
    
    if(real[which(real[,2] == j), 1] %in% set){
      out[j] <- 1
    }else{
      out[j] <- 0
    }
  
  }
  
  mean(out[2:ncol(m)])
}

## Get hit rate
get_hit_rate <- function(truth, adj){
  
  bin_adj <- matrix(0, nrow = nrow(adj), ncol = ncol(adj))
  bin_adj[cbind(colMaxs(adj)[2:ncol(adj)], 2:ncol(adj))] <- 1
  
  related_hit <- bin_adj[truth]
  mean(related_hit)
}

get_cov(truth + 1, adj)
get_cov(truth + 1, adj_badtrip)
get_cov(truth + 1, adj_o2)
get_cov(truth + 1, adj_cont)

get_hit_rate(truth + 1, adj)
get_hit_rate(truth + 1, adj_badtrip)
get_hit_rate(truth + 1, adj_o2)
get_hit_rate(truth + 1, adj_cont)


adj_true <- matrix(0, nrow = length(cons) + 1, ncol = length(cons) + 1)
adj_true[truth + 1] <- 1

g_true <- graph_from_adjacency_matrix(adj_true[2:(N+1), 2:(N+1)])
coords <- layout_as_tree(g_true)



g1 <- ggraph(create_layout(g_true, coords)) + 
  geom_edge_link(color = '#55BB22', width = 1, arrow = arrow(type = "closed", length = unit(3, 'mm')), 
                 start_cap = circle(5, 'mm'),
                 end_cap = circle(5, 'mm')) + 
  geom_node_point(size = 6, color = "#555555") +
  geom_node_text(aes(label = V(g_true)), color = "white", size = 4) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
  panel.background = element_blank(), axis.line = element_blank())

print(g1)

graph_to_df <- function(g){
  es <- as_edgelist(g)
  df <- as.data.frame(es)
  
  for (i in 1:nrow(df)) {
    sub <- truth[truth[,1] == es[i,1] & truth[,2] == es[i,2], , drop = F]
    if(nrow(sub) > 0){
      df$true[i] <- T
    }else{
      df$true[i] <- F
    }
  }
  
  df <- cbind(df, E(g)$weight)
  
  colnames(df) <- c("from", "to", "true", "probability")
  
  layout <- create_layout(graph_from_data_frame(df), coords)
  layout <- create_layout(graph_from_data_frame(df), coords[as.numeric(layout$name), ])
  
  g2 <- ggraph(layout) + 
    geom_edge_link(aes(color = true, alpha = probability), width = 1, arrow = arrow(type = "closed", length = unit(3, 'mm')), 
                   start_cap = circle(5, 'mm'),
                   end_cap = circle(5, 'mm')) + 
    geom_node_point(size = 6, color = "#555555") +
    geom_node_text(aes(label = name), color = "white", size = 4) + 
    scale_edge_color_manual(values= c("black", '#55BB22')) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_blank(), legend.position = "none")
  g2
}

# Our tool
g <- graph_from_adjacency_matrix(adj[2:(N+1), 2:(N+1)], weighted = T)
g2 <- graph_to_df(g)
print(g2)

# Badtrip
g_badtrip <- graph_from_adjacency_matrix(adj_badtrip[2:(N+1), 2:(N+1)], weighted = T)
g3 <- graph_to_df(g_badtrip)
print(g3)

### Outbreaker

g_o2 <- graph_from_adjacency_matrix(adj_o2[2:(N+1), 2:(N+1)], weighted = T)
g4 <- graph_to_df(g_o2)
print(g4)

### Consensus alone

g_cont <- graph_from_adjacency_matrix(adj_cont[2:(N+1), 2:(N+1)], weighted = T)
g5 <- graph_to_df(g_cont)
print(g5)

row1 <- plot_grid(g1, g2, g3, g4, g5, labels = LETTERS[1:5], ncol = 5)
row1
