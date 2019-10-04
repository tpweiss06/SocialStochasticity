# This script will evaluate the effect of stochastic group formation on expected
#    population growth. Since it deals with *expected* growth, demographic stochasticity
#    is irrelevant. Similarly, it will start with a fixed population of males and
#    females so stochasticity in the sex ratio is irrelevant.

library(here)
setwd(here())
source("ModelFunctions.R")
source("CommonParams.R")

# Set the starting population size and the number of simulations
Nt <- 30
Nfem <- round(z*Nt)
Nmale <- Nt - Nfem
Nsims <- 10000

# Set up result objects to hold the expected population sizes
DetermNhat <- array(data = NA, dim = c(length(TauSeq), 3))
StochNhat <- array(data = NA, dim = c(length(TauSeq), 3, Nsims))

# Loop through everything to populate the objects
TotalIterations <- length(TauSeq) * 3 * Nsims
Counter <- 1
for(j in 1:length(TauSeq)){
     for(ms in 1:3){
          if(ms == 1){
               Fg <- z * TauSeq[j]
               R <- Fg * R0
          } else{
               R <- R0
          }
          ########## First calculate the expected population size for a deterministic sex ratio
          # Form groups and calculate Nhat
          Groups <- FormGroups(Tau = TauSeq[j], Females = Nfem, Males = Nmale, 
                               Stochastic = FALSE, z = z)
          Ng <- dim(Groups)[2]
          Nf <- Groups[1,]
          Nm <- Groups[2,]
          AllPop <- sum(Nf) + sum(Nm)
          # Determine the number of mating females in the population
          MatedFems <- MatingFemales(ms = ms, Ng = Ng, GroupMales = Nm, GroupFemales = Nf)
          # Now calculate and store Nhat
          DetermNhat[j,ms] <- MatedFems * (R/z) * exp(-1 * alpha * AllPop)
          
          ########## Now simulate stochastic sex ratios and store those values
          for(s in 1:Nsims){
               # Simulate a random number of females and males
               Stoch_Nfem <- sum(rbinom(n = Nt, size = 1, prob = z))
               Stoch_Nmale <- Nt - Stoch_Nfem
               # Now deterministically form groups and calculate Nhat
               Groups <- FormGroups(Tau = TauSeq[j], Females = Stoch_Nfem, Males = Stoch_Nmale, 
                                    Stochastic = FALSE, z = z)
               Ng <- dim(Groups)[2]
               Nf <- Groups[1,]
               Nm <- Groups[2,]
               AllPop <- sum(Nf) + sum(Nm)
               MatedFems <- MatingFemales(ms = ms, Ng = Ng, GroupMales = Nm, GroupFemales = Nf)
               StochNhat[j,ms,s] <- MatedFems * (R/z) * exp(-1 * alpha * AllPop)
               
               # Print out a progress indicator
               Progress <- round((Counter / TotalIterations) * 100, digits = 1)
               print(paste(Progress, "% complete", sep = ""))
               
               Counter <- Counter + 1
          }
     }
}


# Now calculate the difference between the deterministic and stochastic results
DiffResults <- array(data = NA, c(length(TauSeq), 3, Nsims))
DiffSummary <- array(data = NA, c(length(TauSeq), 3, 3))
DiffStdDev <- array(data = NA, c(length(TauSeq), 3))

for(j in 1:length(TauSeq)){
     for(ms in 1:3){
          for(n in 1:Nsims){
               DiffResults[j,ms,n] <- StochNhat[j,ms,n] - DetermNhat[j,ms]
          }
          DiffSummary[j,ms,1] <- mean(DiffResults[j,ms,])
          DiffSummary[j,ms,2:3] <- quantile(DiffResults[j,ms,], probs = c(0.25, 0.75))
          DiffStdDev[j,ms] <- sd(DiffResults[j,ms,])
     }
}

# Finally, plot the results together
# Location
InnerMar <- c(0.5, 0.5, 0.5, 0.5)
OuterMar <- c(3, 5, 0.5, 0.5)
LetterAdj <- 0.025
LetterLine <- -1.25
xLabLine <- 2
yLabLine1 <- 3.5
yLabLine2 <- 2
# Size
LineWidth <- 1.25
AxisSize <- 1
LabelSize <- 1
# Labels
LetterVec <- c("a", "b", "c")
# Range
xRange <- c(1, 25)
yRange <- c(-20, 10)

FigName <- "Figures/EffectOfSexRatio.pdf"
FigWidth <- 6
FigHeight <- 2
pdf(file = FigName, width = FigWidth, height = FigHeight, onefile = FALSE, paper = "special")
     par(mar = InnerMar, oma = OuterMar, mfrow = c(1,3))
     # Plot the graphs in order of mating system
     for(ms in 1:3){
          plot(x = TauSeq, y = DiffSummary[,ms,1], xlim = xRange, ylim = yRange, 
               main = "", xlab = "", ylab = "", yaxt = "n", xaxt = "n")
          segments(x0 = TauSeq, y0 = DiffSummary[,ms,1] + DiffStdDev[,ms], x1 = TauSeq,
                   y1 = DiffSummary[,ms,1] - DiffStdDev[,ms])
          # Add the axes
          axis(side = 1, at = seq(from = 2, to = 24, by = 4), cex.axis = AxisSize)
          axis(side = 1, at = TauSeq, labels = FALSE, tcl = -0.25)
          if(ms == 1){
               axis(side = 2, at = seq(from = -20, to = 10, by = 5), las = 1, cex.axis = AxisSize)
               axis(side = 2, at = seq(from = -20, to = 10, by = 1), tcl = -0.25, labels = FALSE)
          } else{
               axis(side = 2, at = seq(from = -20, to = 10, by = 5), las = 1, labels = FALSE)
               axis(side = 2, at = seq(from = -20, to = 10, by = 1), tcl = -0.25, labels = FALSE)
          }
          # Add a horizontal line at 0
          abline(h = 0, col = "grey", lwd = LineWidth)
          # Add the graph letter
          mtext(LetterVec[ms], side = 3, line = LetterLine, adj = LetterAdj,
                cex = LabelSize)
     }
     # Add the axis labels
     mtext("Group size", side = 1, outer = TRUE, line = xLabLine, cex = LabelSize)
     mtext("Difference in", side = 2, outer = TRUE, line = yLabLine1, cex = LabelSize)
     mtext(expression(paste("expected ", N[t+1])), side = 2, outer = TRUE, line = yLabLine2, cex = LabelSize)
dev.off()



