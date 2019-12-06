# This script will run extinction simulations on MSI for the project with Damon

# Set the number of processors to use
nProc <- 24*2

# Load the necessary libraries and set the working directory
library(parallel)
library(Rmpi)
setwd("~/SocialStochasticity/")

# Load the necessary scripts
source("CommonParams.R")

# Set the parameter values for the upcoming simulations
NumSims <- 10000
NumGens <- 1000
F_init <- 12
M_init <- 12
MateSysVec <- c(1,2,3)

# Create vectors to pass to the nodes with simulation information
SimTauVec <- rep(TauSeq, each = length(MateSysVec))
SimMateSys <- rep(MateSysVec, length(TauSeq))
SimVec <- 1:length(SimTauVec)

SimFunction <- function(i){
     # Create empy objects to hold abundance vectors, a boolean extinction indicator,
     #    and the time to extinction for each simulation
     TempResults <- matrix(NA, nrow = 4, ncol = NumGens)
     ExtResults <- matrix(0, nrow = 4, ncol = NumSims)
     ExtTime <- matrix(NA, nrow = 4, ncol = NumSims)
     # Set the group size and determine the R value based on mating system
     Tau <- SimTauVec[i]
     ms <- SimMateSys[i]
     # Calculate R from R0 based on the mating system
     if(ms == 1){
          Fg <- z * Tau
          R <- Fg * R0
     } else{
          R <- R0
     }
     for(j in 1:NumSims){
          # Demographic stochasticity model
          TempResults[1,] <- Demographic(ms = ms, R = R, alpha = alpha, z = z, Tau = Tau, 
                                           NumGens = NumGens, F_init = F_init, M_init = M_init, 
                                           ReturnN = TRUE, RemainderMortality = FALSE)
          # Sex ratio stochasticity model
          TempResults[2,] <- SexRatio(ms = ms, R = R, alpha = alpha, z = z, Tau = Tau, 
                                        NumGens = NumGens, F_init = F_init, M_init = M_init, 
                                        ReturnN = TRUE, RemainderMortality = FALSE)
          # Group formation stochasticity model
          TempResults[3,] <- GroupFormation(ms = ms, R = R, alpha = alpha, z = z, Tau = Tau, 
                                              NumGens = NumGens, F_init = F_init, M_init = M_init, 
                                              ReturnN = TRUE, RemainderMortality = FALSE)
          # Full model
          TempResults[4,] <- Full(ms = ms, R = R, alpha = alpha, z = z, Tau = Tau, 
                                    NumGens = NumGens, F_init = F_init, M_init = M_init, 
                                    ReturnN = TRUE, RemainderMortality = FALSE)
          # Now store the extinction results and times
          for(m in 1:4){
               if(TempResults[m,NumGens] == 0){
                    ExtResults[m,j] <- 1
                    ExtTime[m,j] <- sum(TempResults[m,] != 0) + 1
               }
          }
     }
     # Finally calculate the proportion of simulations to go extinct 
     ExtProp <- rowSums(ExtResults) / NumSims
     # Put both objects together and return them as a list
     ExtResults <- list(ExtProp = ExtProp, ExtTime = ExtTime)
     return(ExtResults)
}

###### Now create the cluster and run the simulations
cl <- makeCluster(nProc - 1, type = "MPI")

# Export the necessary objects to each node
clusterExport(cl, c("NumGens", "NumSims", "SimTauVec", "SimMateSys", "R0", "z", "alpha", "F_init", "M_init") )

# Change the working directory of the worker nodes
source("CommonParams.R")
source("ModelFunctions.R")
temp1 <- clusterEvalQ(cl, source("~/SocialStochasticity/CommonParams.R") )
temp2 <- clusterEvalQ(cl, source("~/SocialStochasticity/ModelFunctions.R"))

# Run the simulations
SimData <- clusterApply(cl, x = SimVec, fun = SimFunction)

# Now extract the data and save it in a format for later graphing
ExtProb <- array(NA, dim = c(length(MateSysVec), 4, length(TauSeq)))
ExtTime <- array(NA, dim = c(length(MateSysVec), 4, length(TauSeq), NumSims))

for(i in 1:length(SimData)){
     CurMateSys <- SimMateSys[i]
     CurTau <- SimTauVec[i]
     TauIndex <- which(TauSeq == CurTau)
     ExtProb[CurMateSys, , TauIndex] <- SimData[[i]]$ExtProp
     ExtTime[CurMateSys, , TauIndex,] <- SimData[[i]]$ExtTime
}

# Save the output
save(ExtProb, ExtTime, file = "~/SocialStochasticity/ExtDynamics.rdata")



