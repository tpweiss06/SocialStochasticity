# This script will create a graph with temporal population dynamics for a given
#    group size (10) and mating system (monogamous) to use as an example figure.

# Get that function from Lauren to set the working directory. I think that will
#    be the best option.

# Source the function scripts
library(here)
setwd(here())
source("ModelFunctions.R")
source("CommonParams.R")

# Set some initial parameter values
NumGens <- 100
F_init <- 12
M_init <- 12
Tau <- TauSeq[5]
ms <- 3

# Create objects to hold the results for each model
PopSize <- matrix(NA, nrow = 5, ncol = NumGens)
rownames(PopSize) <- c("Deterministic", "Demographic", "SexRatio", "Group Formation", "Full")

# First perform the model runs
PopSize[1,] <- Deterministic(ms = ms, R = R0, alpha = alpha, z = z, Tau = Tau, 
                             NumGens = NumGens, F_init = F_init, M_init = M_init, 
                             ReturnN = TRUE, RemainderMortality = FALSE)
PopSize[2,] <- Demographic(ms = ms, R = R0, alpha = alpha, z = z, Tau = Tau, 
                           NumGens = NumGens, F_init = F_init, M_init = M_init, 
                           ReturnN = TRUE, RemainderMortality = FALSE)
PopSize[3,] <- SexRatio(ms = ms, R = R0, alpha = alpha, z = z, Tau = Tau, 
                        NumGens = NumGens, F_init = F_init, M_init = M_init, 
                        ReturnN = TRUE, RemainderMortality = FALSE)
PopSize[4,] <- GroupFormation(ms = ms, R = R0, alpha = alpha, z = z, Tau = Tau, 
                              NumGens = NumGens, F_init = F_init, M_init = M_init, 
                              ReturnN = TRUE, RemainderMortality = FALSE)
PopSize[5,] <- Full(ms = ms, R = R0, alpha = alpha, z = z, Tau = Tau, 
                    NumGens = NumGens, F_init = F_init, M_init = M_init, 
                    ReturnN = TRUE, RemainderMortality = FALSE)

# Now create the figure
# Size
DetermLineWidth <- 1.5
StochLineWidth <- 1
AxisSize <- 1
LabSize <- 1.25
LegSize <- 1
# Range
xRange <- c(1, NumGens)
yRange <- c(-5, 130)
# Color and line type
library(RColorBrewer)
ModelColorSeq <- brewer.pal(n = 4, name = "Dark2")
ModelLineSeq <- 2:5
pdf(file = "Figures/ExampleRuns.pdf", width = 6, height = 4, onefile = FALSE, paper = "special")
     par(mar = c(4, 4, 1, 1))
     plot(x = 1:NumGens, y = PopSize[1,], type = "l", lwd = DetermLineWidth, col = "black",
          las = 1, xlab = "Time", ylab = "Abundance", main = "", cex.axis = AxisSize,
          cex.lab = LabSize, xlim = xRange, ylim = yRange)
     # Now add the stochastic population sizes
     for(i in 1:4){
          lines(x = 1:NumGens, y = PopSize[1+i,], lty = ModelLineSeq[i], col = ModelColorSeq[i],
                lwd = StochLineWidth)
     }
     # Finally add a legend
     legend(x = 15, y = 30, legend = row.names(PopSize), col = c("black", ModelColorSeq),
            lty = c(1, ModelLineSeq), lwd = 1.25, cex = LegSize, bty = "n", inset = 0,
            ncol = 2)
dev.off()




