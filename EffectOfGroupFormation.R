# This script will evaluate the effect of stochastic group formation on expected
#    population growth. Since it deals with *expected* growth, demographic stochasticity
#    is irrelevant. Similarly, it will start with a fixed population of males and
#    females so stochasticity in the sex ratio is irrelevant.

library(here)
setwd(here())
source("ModelFunctions.R")
source("CommonParams.R")
NumSims <- 10000

# Set the starting population composisiton for three different populations:
#    female biased, even, and male biased
InitPop <- matrix(data = c(20, 10, 15, 15, 10, 20), nrow = 3, ncol = 2, byrow = TRUE)
colnames(InitPop) <- c("nFem", "nMale")

# Set up result objects to hold the expected population sizes
DetermNhat <- array(data = NA, dim = c(nrow(InitPop), length(TauSeq), 3))
StochNhat <- array(data = NA, dim = c(nrow(InitPop), length(TauSeq), 3, NumSims))

# Loop through everything to populate the objects
TotalIterations <- nrow(InitPop) * length(TauSeq) * 3 * NumSims
Counter <- 1
for(i in 1:nrow(InitPop)){
     for(j in 1:length(TauSeq)){
          for(ms in 1:3){
               if(ms == 1){
                    Fg <- z * TauSeq[j]
                    R <- Fg * R0
               } else{
                    R <- R0
               }
               ########## First calculate the expected population size for deterministic group formation
               DetermGroups <- FormGroups(Tau = TauSeq[j], Females = InitPop[i,1],
                                          Males = InitPop[i,2], Stochastic = FALSE,
                                          z = z, RemainderMortality = FALSE)
               DetermNg <- dim(DetermGroups)[2]
               DetermNf <- DetermGroups[1,]
               DetermNm <- DetermGroups[2,]
               AllPop <- sum(DetermNf) + sum(DetermNm)
               # Determine the number of mating females in the population
               MatedFems <- MatingFemales(ms = ms, Ng = DetermNg, GroupMales = DetermNm, 
                                          GroupFemales = DetermNf)
               # Now calculate and store Nhat
               GroupSizes <- DetermNf + DetermNm
               GroupAdjust <- GroupSizes / TauSeq[j]
               DetermNhat[i,j,ms] <- (sum((MatedFems*GroupAdjust*R)) / z) * exp(-1 * alpha * AllPop)
               
               ########## Now simulate stochastic group formation and store those values
               for(s in 1:NumSims){
                    StochGroups <- FormGroups(Tau = TauSeq[j], Females = InitPop[i,1],
                                              Males = InitPop[i,2], Stochastic = TRUE,
                                              RemainderMortality = FALSE)
                    StochNg <- dim(StochGroups)[2]
                    StochNf <- StochGroups[1,]
                    StochNm <- StochGroups[2,]
                    AllPop <- sum(StochNf) + sum(StochNm)
                    # Determine the number of mating females in the population
                    MatedFems <- MatingFemales(ms = ms, Ng = StochNg, GroupMales = StochNm, 
                                               GroupFemales = StochNf)
                    # Now calculate and store Nhat
                    GroupSizes <- StochNf + StochNm
                    GroupAdjust <- GroupSizes / TauSeq[j]
                    StochNhat[i,j,ms,s] <- (sum((MatedFems*GroupAdjust*R)) / z) * exp(-1 * alpha * AllPop)
                    
                    # Print out a progress indicator
                    Progress <- round((Counter / TotalIterations) * 100, digits = 1)
                    print(paste(Progress, "% complete", sep = ""))
                    
                    Counter <- Counter + 1
               }
          }
     }
}

# Now calculate the difference between the deterministic and stochastic results
DiffResults <- array(data = NA, c(nrow(InitPop), length(TauSeq), 3, NumSims))
DiffSummary <- array(data = NA, c(nrow(InitPop), length(TauSeq), 3, 3))
DiffStdDev <- array(data = NA, c(nrow(InitPop), length(TauSeq), 3))
for(i in 1:nrow(InitPop)){
     for(j in 1:length(TauSeq)){
          for(ms in 1:3){
               for(n in 1:NumSims){
                    DiffResults[i,j,ms,n] <- StochNhat[i,j,ms,n] - DetermNhat[i,j,ms]
               }
               DiffSummary[i,j,ms,1] <- mean(DiffResults[i,j,ms,])
               DiffSummary[i,j,ms,2:3] <- quantile(DiffResults[i,j,ms,], probs = c(0.25, 0.75))
               DiffStdDev[i,j,ms] <- sd(DiffResults[i,j,ms,])
          }
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
LegendSize <- 1
FigWidth <- 6
FigHeight <- 2
# Range
xRange <- c(1, 25)
yRange <- c(-35, 35)
# Color
CompColSeq <- c(rgb(red = 0.9, green = 0.6, blue = 0), rgb(red = 0, green = 0, blue = 0),
                rgb(red = 0, green = 0.45, blue = 0.7))
# Labels
LetterVec <- c("a", "b", "c")
# Offset x coordinates to reduce overlap
Adj <- 0.25
xSeqs <- matrix(NA, nrow = nrow(InitPop), ncol = length(TauSeq))
xSeqs[1,] <- TauSeq - Adj 
xSeqs[2,] <- TauSeq
xSeqs[3,] <- TauSeq + Adj

FigName <- "Figures/EffectOfGroupFormation.pdf"
pdf(file = FigName, width = FigWidth, height = FigHeight, onefile = FALSE, paper = "special")
     par(mar = InnerMar, oma = OuterMar, mfrow = c(1,3))
     # Plot the graphs in order of mating system
     for(ms in 1:3){
          plot(x = NA, y = NA, xlim = xRange, ylim = yRange, main = "", xlab = "",
               ylab = "", yaxt = "n", xaxt = "n")
          for(i in 1:nrow(InitPop)){
               points(x = xSeqs[i,], y = DiffSummary[i,,ms,1], col = CompColSeq[i])
               #segments(x0 = xSeqs[i,], y0 = DiffSummary[i,,ms,2], x1 = xSeqs[i,],
               #        y1 = DiffSummary[i,,ms,3], col = CompColSeq[i])
               segments(x0 = xSeqs[i,], y0 = DiffSummary[i,,ms,1] + DiffStdDev[i,,ms], x1 = xSeqs[i,],
                       y1 = DiffSummary[i,,ms,1] - DiffStdDev[i,,ms], col = CompColSeq[i])
          }
          # Add the axes
          axis(side = 1, at = seq(from = 2, to = 24, by = 4), cex.axis = AxisSize)
          axis(side = 1, at = TauSeq, labels = FALSE, tcl = -0.25)
          if(ms == 1){
               axis(side = 2, at = seq(from = -35, to = 35, by = 15), las = 1, cex.axis = AxisSize)
               axis(side = 2, at = seq(from = -35, to = 35, by = 5), tcl = -0.25, labels = FALSE)
          } else{
               axis(side = 2, at = seq(from = -35, to = 35, by = 15), las = 1, labels = FALSE)
               axis(side = 2, at = seq(from = -35, to = 35, by = 5), tcl = -0.25, labels = FALSE)
          }
          # Add a horizontal line at 0
          abline(h = 0, col = "grey", lwd = LineWidth)
          # Add the graph letter
          mtext(LetterVec[ms], side = 3, line = LetterLine, adj = LetterAdj,
                cex = LabelSize)
          # Add the legend to the middle plot
          if(ms == 2){
               legend(x = "topright", bty = "n", legend = c("Female biased", "Even", "Male biased"),
                      col = CompColSeq, cex = LegendSize, pch = 1)
          }
     }
     # Add the axis labels
     mtext("Group size", side = 1, outer = TRUE, line = xLabLine, cex = LabelSize)
     mtext("Difference in", side = 2, outer = TRUE, line = yLabLine1, cex = LabelSize)
     mtext(expression(paste("expected ", N[t+1])), side = 2, outer = TRUE, line = yLabLine2, cex = LabelSize)
dev.off()









