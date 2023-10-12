### Investigate all tables of reconstruction results
library(reconstructR)
library(ggplot2)
library(igraph)
library(parallel)
library(Rfast)

setwd("./cluster-analyses/")


amalgam <- data.frame()
dirs <- list.files("./all_reconstructions/")

max_cluster_size <- 200

isnv_info <- c() # What proportion of transmission links have an iSNV signature?
isnv_info_nonzero <- c() # What proportion of transmission links with SNV distance > 0 have an iSNV signature?
cluster_size <- c()
n_facilities <- c()
facility_type <- c() # if only one facility; otherwise NA
has_minor_to_fixed <- c() # does any link with a minor allele -> fixed have >90% posterior probability?
has_fixed_to_minor <- c() # does any link with a fixed allele -> minor have >90% posterior probability?
has_split <- c() # does any link with a split bottleneck have >90% posterior probability?


# Read in metadata to get info on clusters
meta <- read.csv("./metadata.txt", sep = "\t")

for (i in 1:length(dirs)) {
  
  d <- dirs[i]
  
  if("results.csv" %in% list.files(paste0("./all_reconstructions/", d))){
    tab <- read.csv(paste0("./all_reconstructions/", d, "/results.csv"))
    tab[,1] <- d
    
    
    
    if(length(unique(tab$To)) <= max_cluster_size){
      isnv_info[i] <- (sum(tab$Probability[tab$Relationship == "Split Bottleneck"]) +
                         sum(tab$Probability[tab$Relationship == "Minor Allele to Fixation"]) +
                         sum(tab$Probability[tab$Relationship == "Fixation to Minor Allele"])) /
        sum(tab$Probability[tab$Relationship != "Index Case"])
      isnv_info_nonzero[i] <- sum(tab$Probability[tab$Relationship %in% c("Split Bottleneck", "Minor Allele to Fixation", "Fixation to Minor Allele") & tab$Distance > 0]) /
        sum(tab$Probability[tab$Relationship != "Index Case" & tab$Distance > 0])
      cluster_size[i] <- length(unique(tab$To))
      
      # Add a new column stating how many posterior samples the link is present in
      tab$N_Post_Samples <- cluster_size[i] * 10 * 0.9 * tab$Probability
      
      amalgam <- rbind(amalgam, tab)
      
      names <- unique(tab$To)
      facilities <- meta$Practice.Name[match(names, meta$Sequence.Name)]
      types <- meta$Program[match(names, meta$Sequence.Name)]
      n_facilities[i] <- length(unique(facilities))
      
      if(n_facilities[i] == 1){
        facility_type[i] <- unique(types)
      }else{
        facility_type[i] <- NA
      }
      
      
      if(any(tab$Relationship == "Minor Allele to Fixation" & tab$Probability > 0.9)){
        has_minor_to_fixed[i] <- T
      }else{
        has_minor_to_fixed[i] <- F
      }
      
      if(any(tab$Relationship == "Fixation to Minor Allele" & tab$Probability > 0.9)){
        has_fixed_to_minor[i] <- T
      }else{
        has_fixed_to_minor[i] <- F
      }
      
      if(any(tab$Relationship == "Split Bottleneck" & tab$Probability > 0.9)){
        has_split[i] <- T
      }else{
        has_split[i] <- F
      }
      
    }else{
      isnv_info[i] <- NA
      isnv_info_nonzero[i] <- NA
      cluster_size[i] <- NA
      n_facilities[i] <- NA
      has_minor_to_fixed[i] <- NA
      has_fixed_to_minor[i] <- NA
      has_split[i] <- NA
      facility_type[i] <- NA
    }
  }
}


# How many clusters are there?
n_clusters <- sum(!is.na(n_facilities))
print(n_clusters)

# Mean/median size and IQR
mean(cluster_size, na.rm = T)
median(cluster_size, na.rm = T)
quantile(cluster_size, 0.25, na.rm = T)
quantile(cluster_size, 0.75, na.rm = T)

# SD
sd(cluster_size, na.rm = T)


# Calculate the proportion of links with iSNV signature per cluster
mean_signature <- mean(isnv_info, na.rm = T)
print(mean_signature)

# 95% CI by bootstrapping
n_boot <- 1000

get_bootstrap <- function(data){
  n_obs <- length(data)
  thetas <- c()
  for (i in 1:n_boot) {
    samp <- sample(data, n_obs, replace = T)
    thetas[i] <- mean(samp)
  }
  c(quantile(thetas, 0.025), quantile(thetas, 0.975))
}

print(get_bootstrap(isnv_info[!is.na(isnv_info)]))

# Proportion of cases in dataset with 0 SNV distance
prop_snv_zero <- sum(amalgam$Probability[amalgam$Distance==0 & amalgam$Relationship != "Index Case"]) / sum(amalgam$Probability[amalgam$Relationship != "Index Case"])
print(prop_snv_zero)

# Proportion of cases in dataset with nonzero SNV distance
prop_snv_nonzero <- sum(amalgam$Probability[amalgam$Distance>0 & amalgam$Relationship != "Index Case"]) / sum(amalgam$Probability[amalgam$Relationship != "Index Case"])
print(prop_snv_nonzero)

# What's the average posterior probability for transmission links with SNV distance = 0 and a split bottleneck?
mean(amalgam$Probability[amalgam$Relationship == "Split Bottleneck" & amalgam$Distance == 0 & amalgam$Probability > 0.01])

# 95% CI via bootstrap
get_bootstrap(amalgam$Probability[amalgam$Relationship == "Split Bottleneck" & amalgam$Distance == 0 & amalgam$Probability > 0.01])

# What's the average posterior probability for transmission links with SNV distance = 0 and NO split bottleneck?
mean(amalgam$Probability[amalgam$Relationship == "Same Consensus" & amalgam$Probability > 0.01])

# 95% CI via bootstrap
get_bootstrap(amalgam$Probability[amalgam$Relationship == "Same Consensus" & amalgam$Probability > 0.01])

### What proporiton links with a SNV distance >0 have an iSNV signature?
p_signature_nonzero <- mean(isnv_info_nonzero[!is.na(isnv_info_nonzero)])
print(p_signature_nonzero)

# 95% CI
print(get_bootstrap(isnv_info_nonzero[!is.na(isnv_info_nonzero)]))





