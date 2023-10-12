### Show that we don't need to deal with mutations going from template to new virions
library(ggplot2)
library(reshape2)

w <- 1000 # burst size
N_bursts <- 100 # number of bursts to simulate
p <- 3e-6 # probability of mutation per site per cycle
N_start <- 1 # number of inoculating virions

N_tot <- N_start + w*N_bursts # number of virions at end of process

N_white <- 0 # Number of white (mutated) virions in urn

## New approach: what if we literally just compute the PMF, recursively?
# row i, column j is the probability that after i burst events, there are j-1 white particles (may be 0)

pmf <- matrix(0, nrow = N_bursts, ncol = N_tot)
pmf[1,1] <- 1 - p
pmf[1, w + 1] <- p
for (i in 2:N_bursts) {
  for (j in seq(1, N_tot, w)) {
    # Probability of getting j-1 whites after i-1 events
    # pick white ball and it does mutate, or pick black ball and it doesn't mutate
    pmf[i,j] <- pmf[i-1,j] * (((j-1) / (w*(i-1) + 1)) * p + ((1 - (j-1) / (w*(i-1) + 1))) * (1-p))
    if(j - w >= 1){
      pmf[i,j] <- pmf[i,j] +
        # pick white ball and it doesn't mutate, or pick black ball and it does mutate
        pmf[i-1,j-w] * (((j-w-1) / (w*(i-1) + 1)) * (1-p) + ((1 - (j-w-1) / (w*(i-1) + 1))) * (p))
    }
  }
  print(i)
}



binoms <- ((1-p)*dbinom(0:w, w, p) + p*dbinom(w:0, w, p)) # Probability for adding color DIFFERENT than ancestor
binoms_same <- ((1-p)*dbinom(w:0, w, p) + p*dbinom(0:w, w, p)) # Probability for adding color SAME as ancestor

pmf_current <- rep(0, N_tot)
pmf_current[1:(w+1)] <- binoms
for (i in 1:(N_bursts - 1)) {
  pmf_new <- rep(0, N_tot)
  for (j in 1:(1 + i*w)) {
    pmf_new[j:(j+w)] = pmf_new[j:(j+w)] +
      pmf_current[j] * binoms * (1 - (j-1) / (w*i + 1)) +
      pmf_current[j] * binoms_same * ((j-1) / (w*i + 1))
  }
  pmf_current <- pmf_new
  print(i)
}


plot(log(pmf_current))
plot(log(pmf_current[seq(1, N_tot, w)]))
plot(log(pmf[N_bursts, seq(1, N_tot, w)]))

df <- data.frame(Frequency = (seq(1, N_tot, w) - 1)/N_tot, Bidirectional = log(pmf_current[seq(1, N_tot, w)]), Unidirectional = log(pmf[N_bursts, seq(1, N_tot, w)]))
df <- melt(df, id = c("Frequency"))
colnames(df)[2] <- "Model"

ggplot(df, aes(x = Frequency, y = value, color = Model)) +
  geom_line(linewidth = 1) +
  xlim(c(0,0.5)) +
  xlab("Proportion of Mutated Particles") +
  ylab("Log Probability") +
  scale_color_manual(values = c("#BB5522", "#2255BB")) +
  theme_minimal()

ggsave("./figs/supp.png", width = 6, height = 6)

