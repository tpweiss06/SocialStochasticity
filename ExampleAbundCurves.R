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
ms <- 2

# Create objects to hold the results for each model
PopSize <- matrix(NA, nrow = 5, ncol = NumGens)
rownames(PopSize) <- c("Deterministic", "Demographic", "SexRatio", "GroupFormation", "Full")

# First perform the model runs
PopSize[1,] <- Deterministic(ms = ms, R = R0, alpha = alpha, z = z, Tau = Tau, 
                             NumGens = NumGens, F_init = F_init, M_init = M_init, 
                             ReturnN = TRUE)
PopSize[2,] <- Demographic(ms = ms, R = R0, alpha = alpha, z = z, Tau = Tau, 
                           NumGens = NumGens, F_init = F_init, M_init = M_init, 
                           ReturnN = TRUE)
PopSize[3,] <- SexRatio(ms = ms, R = R0, alpha = alpha, z = z, Tau = Tau, 
                        NumGens = NumGens, F_init = F_init, M_init = M_init, 
                        ReturnN = TRUE)
PopSize[4,] <- GroupFormation(ms = ms, R = R0, alpha = alpha, z = z, Tau = Tau, 
                              NumGens = NumGens, F_init = F_init, M_init = M_init, 
                              ReturnN = TRUE)
PopSize[5,] <- Full(ms = ms, R = R0, alpha = alpha, z = z, Tau = Tau, 
                    NumGens = NumGens, F_init = F_init, M_init = M_init, 
                    ReturnN = TRUE)

# Now create the figure
library(RColorBrewer)
ModelColorSeq <- brewer.pal(n = 4, name = "Dark2")
ModelLineSeq <- 2:5
xRange <- c(1, NumGens)
yRange <- c(0, max(PopSize) + 5)
DetermLineWidth <- 1.5
StochLineWidth <- 1
AxisSize <- 1.25
LabSize <- 1.25
LegSize <- 1.25
pdf(file = "Figures/ExampleRuns.pdf", width = 10, height = 6, onefile = FALSE, paper = "special")
     plot(x = 1:NumGens, y = PopSize[1,], type = "l", lwd = DetermLineWidth, col = "black",
          las = 1, xlab = "Time", ylab = "Abundance", main = "", cex.axis = AxisSize,
          cex.lab = LabSize, xlim = xRange, ylim = yRange1)
     # Now add the stochastic population sizes
     for(i in 1:4){
          lines(x = 1:NumGens, y = PopSize[1+i,], lty = ModelLineSeq[i], col = ModelColorSeq[i],
                lwd = StochLineWidth)
     }
     
     # Finally add a legend
     legend(x = "bottom", legend = row.names(PopSize), col = c("black", ModelColorSeq),
            lty = c(1, ModelLineSeq), lwd = 1.25, cex = LegSize, bty = "n", inset = 0)
dev.off()




