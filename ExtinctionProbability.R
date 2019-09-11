# This code will perform 1,000 simulations of each model and mating system at
#    varying group sizes to determine the extinction probability (calculated 
#    over 100 generations) resulting from the different stochastic processes.

# Source the function scripts
library(here)
setwd(here())
source("ModelFunctions.R")
source("CommonParams.R")

# Set the parameter values for the upcoming simulations
NumSims <- 1000
NumGens <- 1000
F_init <- 12
M_init <- 12
MateSysVec <- c(1,2,3)

# Make an array to hold the proportion of simulations that go extinct for each
#    mating system (1st dimension), stochastic model (2nd dimension), and
#    group size (3rd dimension)
ExtProb <- array(NA, dim = c(length(MateSysVec), 4, length(TauSeq)))
# Make another arry with a similar layout, but to hold the time of extinction
#    for simulations that do go extinct
ExtTime <- array(NA, dim = c(length(MateSysVec), 4, length(TauSeq), NumSims))

# Now create a loop counting system to track progress during the simulations
NumScenarios <- length(TauSeq) * length(MateSysVec)
LoopCounter <- 0

# Now loop through each mating system, group size, and z value to perform the simulations
for(i in 1:length(MateSysVec)){
     for(j in 1:length(TauSeq)){
          # Set the group size and determine the R value based on mating system
          Tau <- TauSeq[j]
          ms <- MateSysVec[i]
          if(ms == 1){
               # If the mating system is alpha pair monogamy, calculate the
               #    adjusted R value so the deterministic model dynamics would
               #    be equivalent to the other mating systems. First, calculate
               #    the expected number of females per group and then multiply
               #    that by R0
               Fg <- z * Tau
               R <- Fg * R0
          } else{
               # In all other mating systems, R will simply be R0
               R <- R0
          }
          
          # Create a temporary object to hold all the simulations for each model
          #    and fill in the extinction times as we go
          TempResults <- array(NA, dim = c(4, NumSims, NumGens))
          for(l in 1:NumSims){
               # Demographic stochasticity model
               TempResults[1,l,] <- Demographic(ms = ms, R = R, alpha = alpha, z = z, Tau = Tau, 
                                                NumGens = NumGens, F_init = F_init, M_init = M_init, 
                                                ReturnN = TRUE)
               if(TempResults[1,l,NumGens] == 0){
                    ExtTime[i,1,j,l] <- sum(TempResults[1,l,] != 0) + 1
               }
               # Sex ratio stochasticity model
               TempResults[2,l,] <- SexRatio(ms = ms, R = R, alpha = alpha, z = z, Tau = Tau, 
                                             NumGens = NumGens, F_init = F_init, M_init = M_init, 
                                             ReturnN = TRUE)
               if(TempResults[2,l,NumGens] == 0){
                    ExtTime[i,2,j,l] <- sum(TempResults[2,l,] != 0) + 1
               }
               # Group formation stochasticity model
               TempResults[3,l,] <- GroupFormation(ms = ms, R = R, alpha = alpha, z = z, Tau = Tau, 
                                                   NumGens = NumGens, F_init = F_init, M_init = M_init, 
                                                   ReturnN = TRUE)
               if(TempResults[3,l,NumGens] == 0){
                    ExtTime[i,3,j,l] <- sum(TempResults[3,l,] != 0) + 1
               }
               # Full model
               TempResults[4,l,] <- Full(ms = ms, R = R, alpha = alpha, z = z, Tau = Tau, 
                                         NumGens = NumGens, F_init = F_init, M_init = M_init, 
                                         ReturnN = TRUE)
               if(TempResults[4,l,NumGens] == 0){
                    ExtTime[i,4,j,l] <- sum(TempResults[4,l,] != 0) + 1
               }
          }
          # Now calculate the proportion of simulations that went extinct for
          #    each model and put that in the corresponding entry of the Results
          #    array. 
          for(l in 1:4){
               FinalPopSizes <- TempResults[l,,NumGens]
               ExtProb[i, l, j] <- sum(FinalPopSizes == 0)/NumSims
          }
               
          # Advance the loop counter and print out the progress
          LoopCounter <- LoopCounter + 1
          PercentComplete <- (LoopCounter/NumScenarios) * 100
          Message <- paste(round(PercentComplete, digits = 1), "% complete", sep = "")
          print(Message)
          
     }
}
save(ExtProb, ExtTime, file = "SimData/ExtDynamics.rdata")

# Now make the graphs
load("SimData/ExtDynamics.rdata")
# Create a graph to display the extinction probability under each model and mating system
library(RColorBrewer)
FigName <- "Figures/ExtinctionProbability.pdf"
xRange <- range(TauSeq)
yRange <- c(0,1)
ModelColorSeq <- brewer.pal(n = 4, name = "Dark2")
ModelLineSeq <- 2:5
MatingWordVec <- c("Alpha Pair Monogamy", "Monogamy", "Polygynandry")
MateWordAdj <- 0.95
MateWordLine <- -2.1
AxisSize <- 1.5
LabelSize <- 1.5
xLabLine <- 3
yLabLine <- 4
LegendSize <- 1.5
pdf(file = FigName, width = 8, height = 8, onefile = FALSE, paper = "special")
     par(mfcol = c(3,1), oma = c(6,6,2,2), mar = c(0,0,0,0))
     for(i in 1:length(MateSysVec)){
          plot(x = NA, y = NA, xlim = xRange, ylim = yRange, main = "", 
               ylab = "", xlab = "", axes = FALSE)
          for(j in 1:4){
               lines(x = TauSeq, y = ExtProb[i,j,], col = ModelColorSeq[j], 
                     lty = ModelLineSeq[j])
          }
          # Add an indicator for the mating system corresponding to each graph
          mtext(MatingWordVec[i], side = 3, adj = MateWordAdj, line = MateWordLine)
          # Add the axes to the figure
          box()
          if(i == 3){
               axis(1, cex.axis = AxisSize)
          }else{
               axis(1, labels = FALSE)
          }
          if((i == 1) | (i == 3)){
               axis(2, las = 1, at = seq(0, 1, by = 0.2), cex.axis = AxisSize)
               axis(2, at = seq(0, 1, by = 0.05), labels = FALSE, tcl = -0.25)
          }else{
               axis(2, las = 1, at = seq(0.2, 0.8, by = 0.2), cex.axis = AxisSize)
               axis(2, las = 1, at = seq(0, 1, by = 0.2), labels = FALSE)
               axis(2, at = seq(0, 1, by = 0.05), labels = FALSE, tcl = -0.25)
          }
     }
     # Add y and x axis labels
     mtext("Group size", side = 1, outer = TRUE, line = xLabLine, cex = LabelSize)
     mtext("Extinction probability", side = 2, outer = TRUE, line = yLabLine, cex = LabelSize)
     # Add a legend to the last plot to define the different lines
     legend("right", legend = c("Demographic", "SexRatio", "GroupFormation", "Full"),
            col = ModelColorSeq, lty = ModelLineSeq, bty = "n", cex = LegendSize)
dev.off()

# Create a graph to display the time to extinction for the different models
FigName <- "Figures/ExtinctionTimes.pdf"
xRange <- c(1.5, 24.5)
yRange <- c(0, max(ExtTime, na.rm = TRUE))
Xseq <- matrix(NA, nrow = 4, ncol = length(TauSeq))
Xseq[1,] <- seq(from = 1.7, to = 23.7, by = 2)
Xseq[2,] <- seq(from = 1.9, to = 23.9, by = 2)
Xseq[3,] <- seq(from = 2.1, to = 24.1, by = 2)
Xseq[4,] <- seq(from = 2.3, to = 24.3, by = 2)
ModelPointSeq <- c(17, 15, 18, 16)
pdf(file = FigName, width = 8, height = 8, onefile = FALSE, paper = "special")
     par(mfcol = c(3,1), oma = c(6,6,2,2), mar = c(0,0,0,0))
     for(i in 1:length(MateSysVec)){
          plot(x = NA, y = NA, xlim = xRange, ylim = yRange, main = "", 
               ylab = "", xlab = "", axes = FALSE)
          for(j in 1:4){
               points(x = Xseq[j,], y = rowMeans(ExtTime[i,j,,], na.rm = TRUE),
                      col = ModelColorSeq[j], pch = ModelPointSeq[j])
               # Now add the interquartile ranges
               for(g in 1:length(TauSeq)){
                    IQR <- quantile(ExtTime[i,j,g,], na.rm = TRUE, probs = c(0.25, 0.75))
                    segments(x0 = Xseq[j,g], y0 = IQR[1], x1 = Xseq[j,g], y1 = IQR[2],
                             col = ModelColorSeq[j])
               }
          }
          # Add an indicator for the mating system corresponding to each graph
          mtext(MatingWordVec[i], side = 3, adj = MateWordAdj, line = MateWordLine)
          # Add the axes to the figure
          box()
          if(i == 3){
               axis(1, cex.axis = AxisSize)
          }else{
               axis(1, labels = FALSE)
          }
          if((i == 1) | (i == 3)){
               axis(2, las = 1, at = seq(0, 1, by = 0.2), cex.axis = AxisSize)
               axis(2, at = seq(0, 1, by = 0.05), labels = FALSE, tcl = -0.25)
          }else{
               axis(2, las = 1, at = seq(0.2, 0.8, by = 0.2), cex.axis = AxisSize)
               axis(2, las = 1, at = seq(0, 1, by = 0.2), labels = FALSE)
               axis(2, at = seq(0, 1, by = 0.05), labels = FALSE, tcl = -0.25)
          }
     }
     # Add y and x axis labels
     mtext("Group size", side = 1, outer = TRUE, line = xLabLine, cex = LabelSize)
     mtext("Time to Extinction", side = 2, outer = TRUE, line = yLabLine, cex = LabelSize)
     # Add a legend to the last plot to define the different lines
     legend("right", legend = c("Demographic", "SexRatio", "GroupFormation", "Full"),
            col = ModelColorSeq, pch = ModelPointSeq, bty = "n", cex = LegendSize)
dev.off()
