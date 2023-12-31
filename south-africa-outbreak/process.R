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

library(ape)
library(readxl)
library(reconstructR)
library(igraph)
library(ggplot2)
library(ggraph)
library(outbreaker2)
library(cowplot)

setwd("~/Desktop/outbreak-reconstruction/south-africa-outbreak/")

N_bases <- 29903

cons <- read.FASTA("./raw/nextclade.aligned.fasta")
cons_accession <- sub(".*?\\|", "", names(cons))
cons_accession <- sub("\\|.*", "", cons_accession)
cons_date <- sub(".*\\|", "", names(cons))
cons_date <- as.Date(cons_date)

leger <- read_xlsx("./raw/leger.xlsx", skip = 1)

matching <- match(cons_accession, leger$`Gisaid Accession`)
names(cons) <- leger$Name[matching]

# Reorder cons and cons_date to match leger
cons_date <- cons_date[match(names(cons), leger$Name)]
cons <- cons[match(names(cons), leger$Name)]

# Split up cons into the two outbreaks
ch1 <- which(leger$Outbreak == "CH1")
ch3 <- which(leger$Outbreak == "CH3")

# Write fastas by outbreak
write.FASTA(cons[ch1], "./ch1/input_data/aligned.fasta")

big_vcf <- as.data.frame(read_xlsx("./raw/vcf.xlsx", skip = 1))

# Rename substituting - for _
big_vcf$Sample <- gsub("\\-", "\\_", big_vcf$Sample)

depth <- 1:N_bases

# BadTrIP genetic input
gen_input <- matrix(ncol = 0, nrow = N_bases)

# For badtrip: convert NC_045512.2 to nucleotide (1,2,3,4) at each position
ref <- read.FASTA("./ch1/input_data/ref.fasta")
raw_to_base <- function(r){
  if(r == as.raw(136)){
    "A"
  }else if(r == as.raw(40)){
    "C"
  }else if(r == as.raw(72)){
    "G"
  }else if(r == as.raw(24)){
    "T"
  }else{
    "N"
  }
}
ref_nt <- sapply(ref[[1]], raw_to_base)
ref_nt <- as.numeric(sapply(ref_nt, function(x){match(x, c("A","C","G","T"))}))
relevant <- c()

# Make one vcf for each sample
for (i in 1:length(cons)) {
  sub <- big_vcf[big_vcf$Sample == names(cons)[i], ]

  pos <- sub$POS
  info <- character(0)
  if(nrow(sub) > 0){
    info <- paste0(
      "DP=",
      sub$Depth,
      ";AF=",
      sub$AF,
      ";SB=",
      sub$SB,
      ";DP4=",
      paste(
        sub$`R1+`,
        sub$`R1-`,
        sub$`R2+`,
        sub$`R2-`,
        sep = ","
      )
    )
  }

  vcf <- data.frame(
    CHROM = sub$CHROM,
    POS = pos,
    ID = rep(".", nrow(sub)),
    REF = sub$REF,
    ALT = sub$ALT,
    QUAL = rep(5000, nrow(sub)),
    FILTER = rep("PASS", nrow(sub)),
    INFO = info
  )

  if(i %in% ch1){
    write.table(vcf, file = paste0("./ch1/input_data/vcf/", names(cons)[i], ".vcf"), col.names = F, row.names = F)
    relevant <- c(relevant, pos)
  }

  ## Read depth
  dp <- rep(1000, N_bases)
  missing <- which(!(cons[[i]] %in% as.raw(c(136,40,72,24))))
  dp[missing] <- 0
  dp[pos] <- sub$Depth

  depth <- cbind(depth, dp)

  ### Prepare genetic input for BadTrIP
  newcol <- matrix(0, ncol = 4, nrow = N_bases)
  newcol[cbind(1:N_bases, ref_nt)] <- 1000
  newcol[cbind(1:N_bases, ref_nt)][missing] <- 0
  newcol[cbind(1:N_bases, ref_nt)][pos] <- sub$R1

  which_alt <- as.numeric(sapply(sub$ALT, function(x){match(x, c("A","C","G","T"))}))

  newcol[cbind(pos, which_alt)] <- sub$R2

  gen_input <- cbind(gen_input, paste(newcol[,1], newcol[,2], newcol[,3], newcol[,4], sep = "-"))

}

colnames(depth) <- c("position", names(cons))
depth1 <- depth[, c(1, ch1 + 1)]

write.csv(depth1, "./ch1/input_data/depth.csv", row.names = F)

ids <- leger$`Outbreak ID`[match(names(cons)[ch1], leger$Name)]
facility <- leger$Facility[match(names(cons)[ch1], leger$Name)]

t_test <- as.Date(c(
  "03-27-2020",
  "04-02-2020",
  "04-02-2020",
  "03-26-2020",
  "03-21-2020",
  "04-04-2020",
  "03-29-2020",
  "04-03-2020",
  #"04-03-2020",
  "04-06-2020",
  "04-04-2020",
  "04-04-2020",
  "04-04-2020",
  "04-04-2020",
  "04-04-2020",
  "04-05-2020",
  "04-03-2020",
  "03-23-2020",
  "04-02-2020",
  "03-09-2020",
  #"03-27-2020",
  "03-24-2020",
  "04-01-2020",
  #"04-03-2020",
  "04-02-2020",
  "04-03-2020",
  #"04-01-2020",
  "04-03-2020"#,
  #"04-03-2020"
  #"04-03-2020"
), format = '%m-%d-%Y')

date1 <- cbind(names(cons)[ch1], as.numeric(difftime(t_test, min(t_test), units = "days")))
write.csv(date1,"./ch1/input_data/date.csv", row.names = F)

setwd("./ch1/")

### Badtrip input

#Badtrip input
relevant <- sort(unique(relevant))
gen_input <- gen_input[relevant, ch1]
colnames(gen_input) <- paste0("S", 1:length(ch1))
write.table(gen_input, file = "./badtrip/inputAlignment.txt", sep = "   ", quote = F, row.names = F)

# Badtrip sample input
sample_input <- data.frame(paste0("H", 1:length(ch1)), paste0("S", 1:length(ch1)), as.numeric(difftime(t_test, min(t_test), units = "days")))
write.table(sample_input, file = "./badtrip/inputSample.txt", sep = "   ", quote = F, col.names = F, row.names = F)

# Badtrip
epi_input <- data.frame(paste0("H", 1:length(ch1)), as.numeric(difftime(t_test, min(t_test), units = "days")) - 10, as.numeric(difftime(t_test, min(t_test), units = "days")) + 10)
write.table(epi_input, file = "./badtrip/inputEpi.txt", sep = "   ", quote = F, col.names = F, row.names = F)


### Run reconstructR
res <- run_mcmc(1e4)

burnin <- 0.1
N <- res[[1]]$N
adj_sa <- matrix(0, nrow = N, ncol = N)
for (i in 1:length(res)) {
  if (i > burnin * length(res)) {
    inds <- cbind(res[[i]]$anc[2:N], 2:N)
    adj_sa[inds] <- adj_sa[inds] + 1
  }
}
adj_sa <- adj_sa / sum(1:length(res) > burnin * length(res))

##### Run outbreaker2
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

res_o2 <- as.matrix(out)
res_o2 <- res_o2[, 9:(9+length(dna)-1)]
res_o2 <- res_o2 + 1
res_o2[is.na(res_o2)] <- 1
adj_sa_o2 <- matrix(0, nrow = length(dna)+1, ncol = length(dna)+1)
burnin <- 0.1

for (j in 1:nrow(res_o2)) {
  if (j > burnin * nrow(res_o2)) {
    inds <- cbind(res_o2[j, ], 2:(length(dna)+1))
    adj_sa_o2[inds] <- adj_sa_o2[inds] + 1
  }
}

adj_sa_o2 <- adj_sa_o2 / sum(1:nrow(res_o2) > burnin * nrow(res_o2))

## Run reconstructR, control settings
res_cont <- run_mcmc(
  N_iters = 1e4,
  filters = list(
    omit = character(0), # any samples we should omit?
    coverage = 0.75, # minimum coverage
    sb = Inf, # filter for strand bias (Inf if not using)
    dp4 = 0, # p-value cutoff for fishers exact test on ref-fwd, ref-rev, alt-fwd, alt-rev
    af = 0.5, # filter for minor allele frequency
    call = 10, # minimum reads of the alternate allele to call iSNV
    depth = 100, # minimum total depth
    problem_qtile = 0 # if an iSNV is one of the problem_qtile-most prevalent in the Broad dataset, throw it out
  )
)


burnin <- 0.1

N <- res_cont[[1]]$N
adj_sa_cont <- matrix(0, nrow = N, ncol = N)
for (i in 1:length(res_cont)) {
  if (i > burnin * length(res_cont)) {
    inds <- cbind(res_cont[[i]]$anc[2:N], 2:N)
    adj_sa_cont[inds] <- adj_sa_cont[inds] + 1
  }
}

adj_sa_cont <- adj_sa_cont / sum(1:length(res_cont) > burnin * length(res_cont))

### Finally, run BadTrIP

setwd("./badtrip")

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

adj_sa_badtrip <- matrix(0, nrow = N+1, ncol = N+1)
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

    adj_sa_badtrip[from + 1, to + 1] <- probs

  }
}

#######

dec <- decipher(res)

setwd("~/Desktop/outbreak-reconstruction/south-africa-outbreak/")
truth <- as.data.frame(read_xlsx("./raw/links.xlsx"))
leger_id <- gsub("\\_.*", "", leger$`Outbreak ID`)

truth <- cbind(leger$Name[match(truth[,1], leger_id)], leger$Name[match(truth[,2], leger_id)])

truth_numbers <- cbind(dec$Alias[match(truth[,1], dec$Name)], dec$Alias[match(truth[,2], dec$Name)])

N <- res[[1]]$N
reads <- res[[1]]$reads
all_pos <- res[[1]]$all_pos

### measure accuracy on putative pairs and visualize
get_vcol <- function(s){
  if(s == "MICU"){
    '#BB9977'
  }else if(s == "MW1"){
    '#BB7799'
  }else if(s == "NH"){
    '#77BB99'
  }else if(s == "CICU"){
    '#7799BB'
  }else if(s == "SICU"){
    '#9977BB'
  }else if(s == "NW"){
    '#BBBB77'
  }else{
    '#BBBBBB'
  }
}

get_coords <- function(adj){
  N <- ncol(adj)
  g <- igraph::graph_from_adjacency_matrix(adj[2:N, 2:N], weighted = T)
  layout_nicely(g)
}

coords <- get_coords(adj_sa)

### Compute the coverage probability
get_cov <- function(truth, adj){
  
  out <- c()
  real <- truth
  m <- adj
  for (j in truth[,2]) {
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
      out <- c(out, 1)
    }else{
      out <- c(out, 0)
    }
    
  }
  
  mean(out)
}

## Get hit rate
get_hit_rate <- function(truth, adj){
  
  bin_adj <- matrix(0, nrow = nrow(adj), ncol = ncol(adj))
  bin_adj[cbind(colMaxs(adj)[2:ncol(adj)], 2:ncol(adj))] <- 1
  
  related_hit <- bin_adj[truth]
  mean(related_hit)
}

analyze <- function(adj){
  N <- ncol(adj)

  g <- igraph::graph_from_adjacency_matrix(adj, weighted = T)
  es <- as_edgelist(g)
  traced <- c()
  for (i in 1:nrow(es)) {
    sub <- truth_numbers[truth_numbers[,1] == es[i,1] & truth_numbers[,2] == es[i,2], , drop = F]
    if(nrow(sub) > 0){
      traced[i] <- 1
    }else{
      traced[i] <- 0
    }
  }

  linked <- c()
  for (i in 1:nrow(es)) {
    from <- es[i, 1]
    to <- es[i, 2]
    if(from == 1 | to == 1){
      linked[i] <- 0
    }else{
      if(grepl(facility[from-1], facility[to-1]) | grepl(facility[to-1], facility[from-1])){
        linked[i] <- 1
      }else{
        linked[i] <- 0
      }
    }
  }

  print(mean(adj[truth_numbers]))
  print(get_cov(truth_numbers, adj))
  print(get_hit_rate(truth_numbers, adj))
  
  
  print(sum(E(g)$weight[traced == 1]) / N)
  print(sum(E(g)$weight[linked == 1 | traced == 1]) / N)
  
  # Convert to data frame suitable for ggraph
  df <- as.data.frame(es)
  df$e_type <- 0
  df$e_type[traced == 1] <- 1
  df$e_type[linked == 1 & traced == 0] <- 2
  df$e_type <- factor(df$e_type)
  df$probability <- E(g)$weight
  colnames(df)[1:2] <- c("from", "to")
  
  new_g <- graph_from_data_frame(df)
  new_g <- delete_vertices(new_g, 1)
  
  layout <- create_layout(new_g, coords)
  layout <- create_layout(new_g, coords[as.numeric(layout$name) - 1, ])
  
  facility[grepl(",", facility) | facility == "ED"] <- "Multiple"
  
  g2 <- ggraph(layout) + 
    geom_edge_link(aes(color = e_type, alpha = probability), width = 0.5, arrow = arrow(type = "closed", length = unit(1.5, 'mm')), 
                   start_cap = circle(3, 'mm'),
                   end_cap = circle(3, 'mm')) + 
    geom_node_point(aes(color = facility[as.numeric(layout$name) - 1]), size = 6) +
    geom_node_text(aes(label = gsub("B", "", leger_id[ch1][as.numeric(layout$name) - 1])), color = "white", size = 2) + 
    scale_edge_color_manual(values= c("black", '#55BB22', '#2255BB')) +
    scale_color_manual(values = c('#7799BB', '#BB9977', '#BBBBBB', '#BB7799', '#77BB99','#BBBB77', '#9977BB')) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_blank(), legend.position = "none")
  
  g2

}


adj_sa_true <- matrix(0, nrow = N, ncol = N)
adj_sa_true[truth_numbers] <- 1

h <- igraph::graph_from_adjacency_matrix(adj_sa_true[2:ncol(adj_sa_true), 2:ncol(adj_sa_true)], weighted = T)
layout <- create_layout(h, coords)
facility[grepl(",", facility) | facility == "ED"] <- "Multiple"
h1 <- ggraph(layout) + 
  geom_edge_link(color = '#55BB22', width = 0.5, arrow = arrow(type = "closed", length = unit(1.5, 'mm')), 
                 start_cap = circle(3, 'mm'),
                 end_cap = circle(3, 'mm')) + 
  geom_node_point(aes(color = facility), size = 6) +
  geom_node_text(aes(label = gsub("B", "", leger_id[ch1])), color = "white", size = 2) + 
  scale_color_manual(values = c('#7799BB', '#BB9977', '#BBBBBB', '#BB7799', '#77BB99','#BBBB77', '#9977BB')) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_blank(), legend.position = "none", plot.title = element_text(hjust = 0.5))

h2 <- analyze(adj_sa)

h3 <- analyze(adj_sa_badtrip)

h4 <- analyze(adj_sa_o2)

h5 <- analyze(adj_sa_cont)

fig4 <- plot_grid(g1,g2,g3,g4,g5,h1,h2,h3,h4,h5, labels = LETTERS[1:10], ncol = 5)
ggsave("~/Desktop/outbreak-reconstruction/south-africa-outbreak/figs/fig4.pdf", width = 12, height = 6)
ggsave("~/Desktop/outbreak-reconstruction/south-africa-outbreak/figs/fig4.png", width = 12, height = 6)

#






