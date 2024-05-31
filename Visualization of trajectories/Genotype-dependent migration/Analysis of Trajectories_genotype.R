Gillespie_fct_geno <- function(Nit, fA, gA, NA10, NA20, fB, gB, NB10, NB20, K1, K2, mu, mWT, mM, theta1, theta2, tlist) {
  
  NA1list <- matrix(NA, nrow = length(tlist), ncol = Nit) # Vector the population size of A individuals in deme 1
  NB1list <- matrix(NA, nrow = length(tlist), ncol = Nit) # Vector the population size of B individuals in deme 1
  NA2list <- matrix(NA, nrow = length(tlist), ncol = Nit) # Vector the population size of A individuals in deme 2
  NB2list <- matrix(NA, nrow = length(tlist), ncol = Nit) # Vector the population size of B individuals in deme 2
  
  for (i in 1:Nit) { # Loop over the number of replicates
    #  Initialization:  
    NA1 <- NA10  # Initialization of the number of A individuals in deme 1
    NB1 <- NB10  # Initialization of the number of B individuals in deme 1
    NA2 <- NA20  # Initialization of the number of A individuals in deme 2
    NB2 <- NB20  # Initialization of the number of B individuals in deme 2
    t <- 0       # Initialization of time
    q <- 1       # Index for the vector of time points
    
    NA1list[q, i] <- NA1  #Vector of the population size of A individuals in deme 1
    NB1list[q, i] <- NB1  #Vector of the population size of B individuals in deme 1
    NA2list[q, i] <- NA2  #Vector of the population size of A individuals in deme 2
    NB2list[q, i] <- NB2  #Vector of the population size of B individuals in deme 2
    q <- q + 1
    cumul <- rep(0, 14)   # initialize sampling tower vector
    
    while ((NA1 + NA2 + NB1 + NB2 != 0) && t < max(tlist)) { # running untill tmax reached or no individuals left
      
      # Compute the transition rates
      # deme 1
      if (t < theta1) {
        TrepA1 <- fA * NA1 * (1 - mu) # new number of A individuals in deme 1 (not yet deteriorated)
        TmutA1 <- fA * NA1 * mu       # new number of B individuals in deme 1
      } else {
        TrepA1 <- 0                   # once deme 1 deteriorated, no new A individuals (& no new mutations)
        TmutA1 <- 0
      }
      
      TrepB1 <- fB * NB1              # new number of B individuals in deme 1 (from existing B)
      
      TdeathA1 <- gA * NA1 * (NA1 + NB1) / K1 # number of dead A, according to total pop. size and carrying capacity of deme 1
      TdeathB1 <- gB * NB1 * (NA1 + NB1) / K1 # number of dead B, according to total pop. size and carrying capacity of deme 1
      
      
      # deme 2
      if (t < theta2) {
        TrepA2 <- fA * NA2 * (1 - mu) # new number of A individuals in deme 2 (not yet deteriorated)
        TmutA2 <- fA * NA2 * mu       # new number of B individuals in deme 2 (not yet deteriorated) ????
      } else {
        TrepA2 <- 0                   # once deme 2 deteriorated, no new A individuals (& no new mutations)
        TmutA2 <- 0
      }
      
      TrepB2 <- fB * NB2              # new number of B individuals in deme 2 (from existing B)
      
      TdeathA2 <- gA * NA2 * (NA2 + NB2) / K2 # number of dead A, according to total pop. size and carrying capacity of deme 2
      TdeathB2 <- gB * NB2 * (NA2 + NB2) / K2 # number of dead B, according to total pop. size and carrying capacity of deme 2
      
      # Migration dependent on genotype
      TmigA12 <- NA1 * mWT            # number of A migrants from deme 1 to deme 2
      TmigB12 <- NB1 * mM            # number of B migrants from deme 1 to deme 2
      TmigA21 <- NA2 * mWT            # number of A migrants from deme 2 to deme 1
      TmigB21 <- NB2 * mM            # number of B migrants from deme 2 to deme 1
      
      # used later for probabilities for different actions??
      T <- TrepA1 + TmutA1 + TrepB1 + TdeathA1 + TdeathB1 + TrepA2 + TmutA2 + TrepB2 + TdeathA2 + TdeathB2 + TmigA12 + TmigB12 + TmigA21 + TmigB21
      
      # Update time  
      r1 <- runif(1)                 # generate one uniform random number; 0 < r1 > 1
      tau <- 1 / T * log(1 / r1)     # finding interval length (when will the next action happen) 
      t <- t + tau                   # passed time plus new time interval
      
      # update list   
      while ((t > tlist[q]) && (q < length(tlist))) {
        NA1list[q, i] <- NA1         # note number of A individuals in deme 1 at position q in tlist vector
        NB1list[q, i] <- NB1         # note number of B individuals in deme 1 at position q in tlist vector  
        NA2list[q, i] <- NA2         # note number of A individuals in deme 2 at position q in tlist vector
        NB2list[q, i] <- NB2         # note number of B individuals in deme 2 at position q in tlist vector
        q <- q + 1                   # update position in list
      }
      
      # Build a sampling tower  
      ir2 <- 1
      r2 <- runif(1)                # generate one uniform random number; 0 < r1 > 1 
      cumul[1]  <- TrepA1
      cumul[2]  <- TrepA1 + TmutA1
      cumul[3]  <- TrepA1 + TmutA1 + TrepB1
      cumul[4]  <- TrepA1 + TmutA1 + TrepB1 + TdeathA1
      cumul[5]  <- TrepA1 + TmutA1 + TrepB1 + TdeathA1 + TdeathB1
      cumul[6]  <- TrepA1 + TmutA1 + TrepB1 + TdeathA1 + TdeathB1 + TrepA2
      cumul[7]  <- TrepA1 + TmutA1 + TrepB1 + TdeathA1 + TdeathB1 + TrepA2 + TmutA2
      cumul[8]  <- TrepA1 + TmutA1 + TrepB1 + TdeathA1 + TdeathB1 + TrepA2 + TmutA2 + TrepB2
      cumul[9]  <- TrepA1 + TmutA1 + TrepB1 + TdeathA1 + TdeathB1 + TrepA2 + TmutA2 + TrepB2 + TdeathA2
      cumul[10] <- TrepA1 + TmutA1 + TrepB1 + TdeathA1 + TdeathB1 + TrepA2 + TmutA2 + TrepB2 + TdeathA2 + TdeathB2
      cumul[11] <- TrepA1 + TmutA1 + TrepB1 + TdeathA1 + TdeathB1 + TrepA2 + TmutA2 + TrepB2 + TdeathA2 + TdeathB2 + TmigA12
      cumul[12] <- TrepA1 + TmutA1 + TrepB1 + TdeathA1 + TdeathB1 + TrepA2 + TmutA2 + TrepB2 + TdeathA2 + TdeathB2 + TmigA12 + TmigB12
      cumul[13] <- TrepA1 + TmutA1 + TrepB1 + TdeathA1 + TdeathB1 + TrepA2 + TmutA2 + TrepB2 + TdeathA2 + TdeathB2 + TmigA12 + TmigB12 + TmigA21
      cumul[14] <- TrepA1 + TmutA1 + TrepB1 + TdeathA1 + TdeathB1 + TrepA2 + TmutA2 + TrepB2 + TdeathA2 + TdeathB2 + TmigA12 + TmigB12 + TmigA21 + TmigB21  # corresponds to T
      
      # Determine which reaction occurs  
      while (cumul[ir2] < r2 * T) {
        ir2 <- ir2 + 1 
      }
      
      # update the population size according to chosen reaction
      if (ir2 == 1) {
        NA1 <- NA1 + 1         # add one A individual to deme 1 (birth of NA1 individual)
      } else if (ir2 == 2) {
        NB1 <- NB1 + 1         # add one B individual to deme 1 (mutation from A to B)? why not NA1-1?
      } else if (ir2 == 3) {
        NB1 <- NB1 + 1         # add one A individual to deme 1 (birth of NB1 individual)
      } else if (ir2 == 4) {
        NA1 <- NA1 - 1         # subtract one A individual from deme 1 (death of NA1 individual)
      } else if (ir2 == 5) {
        NB1 <- NB1 - 1         # subtract one B individual from deme 1 (death of NB1 individual)
      } else if (ir2 == 6) {
        NA2 <- NA2 + 1         # add one A individual to deme 2 (birth of NA2 individual)
      } else if (ir2 == 7) {
        NB2 <- NB2 + 1         # add one B individual to deme 2 (mutation from A to B)
      } else if (ir2 == 8) {
        NB2 <- NB2 + 1         # add one B individual to deme 2 (birth of NB1 individual)
      } else if (ir2 == 9) {
        NA2 <- NA2 - 1         # subtract one A individual from deme 2 (death of NA2 individual)
      } else if (ir2 == 10) {
        NB2 <- NB2 - 1         # subtract one A individual from deme 2 (death of NB2 individual)
      } else if (ir2 == 11) {
        NA1 <- NA1 - 1         # transfer one A individual from deme 1 to deme 2 (migration 1->2)
        NA2 <- NA2 + 1
      } else if (ir2 == 12) {
        NB1 <- NB1 - 1         # transfer one B individual from deme 1 to deme 2 (migration 1->2)
        NB2 <- NB2 + 1
      } else if (ir2 == 13) {
        NA2 <- NA2 - 1         # transfer one A individual from deme 2 to deme 1 (migration 2->1)
        NA1 <- NA1 + 1
      } else if (ir2 == 14) {
        NB2 <- NB2 - 1         # transfer one B individual from deme 2 to deme 1 (migration 2->1)
        NB1 <- NB1 + 1
      }
    }
    
    #store the population size of each time interval
    while (q <= length(tlist)) {
      NA1list[q, i] <- NA1
      NB1list[q, i] <- NB1
      NA2list[q, i] <- NA2
      NB2list[q, i] <- NB2
      q <- q + 1
    }
  }
  return(list(NA1list, NB1list, NA2list, NB2list))
}




# Run it:
# Parameters

fA <- 0.9     # fA: wild-type birth rate
gA <- 1       # gA: wild-type death rate
NA10 <- 90    # NA10: initial number of wild-type individuals in deme 1
NA20 <- 90    # NA20: initial number of wild-type individuals in deme 2
fB <- 0.9     # fB: mutant birth rate
gB <- 1       # gB: mutant death rate
NB10 <- 0     # NB10: initial number of mutants in deme 1
NB20 <- 0     # NB20: initial number of mutants in deme 2
K1 <- 100     # K1: carrying capacity of deme 1
K2 <- 100     # K2: carrying capacity of deme 2
mu <- 0.0001       # mu: mutation probability upon reproduction
mWT <- 0.01     # mWT: migration rate from deme 1 to deme 2
mM <- 0.05      # mM: migration rate from deme 2 to deme 1
theta1 <- 500 # theta1: time at which deme 1 deteriorates
theta2 <- 1000 # theta2: time at which deme 2 deteriorates
time_step <- 1 # time_step: time interval between each time point
tmax <- 1500  # tmax: maximum time
tlist <- seq(0, tmax, by = time_step) # tlist: vector of time points at which I want to save the population sizes
Nit <- 1      # Nit: number of stochastic replicates




library(ggplot2) # for nice plots
library(dplyr)   # for sorting data

# Define how many plots you want to create
Iter <- 10  

# Loop through iterations
for (i in 1:Iter) {
  # Call function Gillespie_fct_geno with appropriate arguments
  result <- Gillespie_fct_geno(Nit, fA, gA, NA10, NA20, fB, gB, NB10, NB20, K1, K2, mu, mWT, mM, theta1, theta2, tlist)
  
  # Process output data
  df_all <- do.call(rbind, lapply(seq_along(result), function(i) {
    data.frame(
      tlist = seq(0, length(result[[i]]) - 1),
      Population = rep(paste("Group", i), each = length(result[[i]])),
      value = as.vector(result[[i]])
    )
  }))
  
  # Add effective time and replicate columns
  df_all <- df_all %>%
    mutate(effective_time = ((seq_along(tlist) - 1) %% nrow(result[[1]])) * (time_step)) %>%
    mutate(replicate = cumsum(c(1, diff(effective_time) < 0)))
  
  # Generate plot
  plot <- ggplot(df_all, aes(x = effective_time, y = value, color = Population, group = replicate)) +
    labs(x = "Time", y = "# Individuals") +
    scale_color_manual(name = "Group",
                       values = c("darkblue", "red", "lightblue", "pink"),
                       labels = c("NA1", "NB1", "NA2", "NB2")) + 
    theme_minimal() +
    theme(plot.title = element_text(size = 24, hjust = 0.5),
          axis.title = element_text(size = 16),
          axis.text = element_text(size = 14),
          legend.text = element_text(size = 14),
          legend.title = element_text(size = 16),
          legend.key.size = unit(1.5, "lines")) +
    annotate("rect", xmin = 500, xmax = tmax, ymin = 0, ymax = 150, fill = "lightgrey", alpha = 0.5) +
    annotate("rect", xmin = 1000, xmax = tmax, ymin = 0, ymax = 150, fill = "grey", alpha = 0.4) +
    geom_line(size = 0.5) +
    ggtitle("Population trajectory genotype dependent, alpha 5, Tmig 2e0") +
    annotate("text", x = 750, y = 140, label = "deme 1 deteriorated", color = "black", linewidth = 4.5) +
    annotate("text", x = 1250, y = 140, label = "deme 1 & 2 deteriorated", color = "black", linewidth = 4.5)
  
  # Save plot as image
  ggsave(paste0("D:/FS24/BSc/Work in Progress/Plotting trajectories/Analysis of Trajectories/Genotype-dependent/alpha 5 _Tmig 2e0/plot_", i, ".jpg"), plot, width = 10, height = 6)
}
