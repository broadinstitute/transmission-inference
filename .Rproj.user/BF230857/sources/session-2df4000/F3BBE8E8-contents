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

### Fit de novo iSNV frequency model to 100k genomes data

library(data.table)
library(ggplot2)
library(FuzzyNumbers)
library(gridExtra)
library(cowplot)

setwd("./cluster-analyses/")

# Minimum read depth
dp_filter <- 100

# Maximum p-value for deleting observation due to strand bias
sb_filter <- 0.05

# Number of alternate reads needed
alt_read_filter <- 10

# Minimum alternate allele frequency
min_af_filter <- 0.03

# Maximum alternate allele frequency
max_af_filter <- 0.5

# Quantile for problematic sites
problem_qtile <- 0.05

#### Lines 55-106 rely on a large file not stored in this repo, available upon request. 
#### Their purpose is to generate a smaller file, maf.RData, which is stored in this repo.
#### This script can be run from line 106 onwards and produce the same results.

# Merged vcf file
vcf <- fread("./allvariants.tsv")

# Extract strand bias
sb <- vcf$INFO
sb <- gsub(".*;SB=", "", sb)
sb <- sub(";DP4.*", "", sb)
sb <- as.numeric(sb)

# Extract DP4 and process
dp4 <- vcf$INFO
dp4 <- gsub(".*;DP4=", "", dp4)
dp4 <- strsplit(dp4, ",")

dp4 <- sapply(dp4, function(v){
  v <- as.numeric(v)
  c(v[1] + v[2], v[3] + v[4])
})

# Extract position
pos <- vcf$POS

## Obtain a list of problematic sites on the genome, based on those that consistently demonstrate iSNVs
context <- read.csv(system.file("extdata", "context.csv", package = "reconstructR"))
tot_isnvs <- data.frame(MUT = context$aa_change_full, POS = context$POS, COUNT = Rfast::rowsums(as.matrix(context[,4:ncol(context)]), na.rm = T))
N_isnv_cases <- sum(tot_isnvs$COUNT)

# What sites have an unusually-high count?
problem_sites <- tot_isnvs$POS[tot_isnvs$COUNT >= quantile(tot_isnvs$COUNT, 1 - problem_qtile)]
problem_sites <- sort(unique(problem_sites))

ref <- dp4[1, ]
alt <- dp4[2, ]

# Subset by strand bias, read depth, # of ref reads, # of alt reads, problematic sites, and minor allele frequency
dp <- ref + alt
af <- alt / dp
maf <- pmin(af, 1-af)

pass <- which(
  sb < -10*log10(sb_filter) &
    ref >= alt_read_filter &
    alt >= alt_read_filter &
    dp >= dp_filter &
    maf <= max_af_filter &
    maf > min_af_filter &
    !(pos %in% problem_sites)
)

maf <- maf[pass]

save(maf, file = "maf.RData")
load("./maf.RData")

### Histogram and fitted density function

# PDF, CDF, and quantile funtions, using approximation in Methods
approx_pdf <- function(x){
  max_af_filter * min_af_filter / ((max_af_filter - min_af_filter) * x^2)
}

approx_cdf <- function(x){
  (max_af_filter * (min_af_filter - x))/((min_af_filter - max_af_filter) * x)
}

approx_qtile <- function(p){
  min_af_filter * max_af_filter / (max_af_filter + p*(min_af_filter - max_af_filter))
}

df <- data.frame(log(maf))

p1 <- ggplot(df, aes(x = (maf))) +
  geom_histogram(aes(y=after_stat(density)), binwidth = 0.01, boundary = 0.01, color = "white", fill = "#CCCCCC") +
  geom_vline(xintercept=min_af_filter + 0.00001, linetype=3, color = "#BB5522") +
  geom_function(fun = approx_pdf, color = "#2255BB", linewidth = 1.5, linetype = "dashed") +
  xlab("Minor Allele Frequency") +
  ylab("Probability Density") +
  theme_minimal()

print(p1)

ggsave("./figs/hist.png", width = 6, height = 6)
ggsave("./figs/hist.pdf", width = 6, height = 6)

### Histogram with log-transform

round_label <- function(x){
  round(x, 2)
}

# Empirical log density
bw <- 0.01
rights <- seq(min_af_filter + bw, max_af_filter, bw)
lefts <- rights - bw
counts <- c()
for (i in 1:length(rights)) {
  counts[i] <- sum(maf < rights[i] & maf >= lefts[i])
}
counts <- counts/sum(counts)/bw

p2 <- ggplot(data.frame(x = lefts), aes(x=x)) +
  geom_line(aes(y = counts), stat = 'identity', linewidth = 1.5, color = "#CCCCCC") +
  geom_vline(xintercept=min_af_filter + 0.00001, linetype=3, color = "#BB5522") +
  geom_function(fun = approx_pdf, color = "#2255BB", linewidth = 1.5, linetype = "dashed") +
  xlab("Log Minor Allele Frequency") +
  scale_x_continuous(trans = 'log', labels = round_label) +
  scale_y_continuous(trans = 'log', labels = round_label) +
  ylab("Log Probability Density") +
  theme_minimal()

print(p2)

ggsave("./figs/loglog.png", width = 6, height = 6)
ggsave("./figs/loglog.pdf", width = 6, height = 6)



### Q-Q Plot

breaks <- approx_cdf(
  seq(min_af_filter, max_af_filter, 0.01)
)

# Rewrite CDF to correctly compute values below min_af_filter and above max_af_filter

full_cdf <- function(x){
  below <- which(x <= min_af_filter)
  above <- which(x  > max_af_filter)
  other <- which(x > min_af_filter & x <= max_af_filter)

  x[below] <- 0
  x[above] <- 1
  x[other] <- approx_cdf(x[other])

  x
}


ks.test(maf + rnorm(length(maf), 0, 0.000000001), full_cdf)


### CDF comparison
p3 <- ggplot(data.frame(x=maf), aes(x=x)) +
  stat_ecdf(geom = "line", color = "#BB5522", linewidth = 1.5, pad = F) +
  geom_function(fun = full_cdf, color = "#2255BB", linewidth = 1.5, linetype = "dashed", xlim = c(min_af_filter, max_af_filter)) +
  xlab("Minor Allele Frequency") +
  ylab("Cumulative Probability Density") +
  # xlim(c(0,0.01)) +
  # ylim(c(0,0.01)) +
  # scale_x_continuous(trans='log') +
  # scale_y_continuous(trans='log') +
  theme_minimal()

print(p3)

ggsave("./figs/cdf.png", width = 6, height = 6)
ggsave("./figs/cdf.pdf", width = 6, height = 6)


ks.test(maf + rnorm(length(maf), 0, 0.00000001), approx_cdf)


fig2 <- plot_grid(p1, p2, p3, labels = "AUTO", ncol = 3)
ggsave("./figs/fig2.pdf", width = 12, height = 4)
ggsave("./figs/fig2.png", width = 12, height = 4)

