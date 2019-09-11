# This script simply defines a series of common parameter values used across
#    different simulations

DeterministicK <- 100
R0 <- 2.1
# Calculate alpha according to the carrying capacity and growth rate
alpha <- log(R0) / DeterministicK
z <- 0.5

TauSeq <- seq(from = 2, to = 24, by = 2)
