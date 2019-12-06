# Make a new graphing script without insets and with ln(Tm) or ln(1/Tm).
#    Pros of ln(Tm): It might be easier since it's explicitly what the reviewer asked for and won't require new explanations
#    Pros of ln(1/Tm): It's not quite what they asked for, but intuitively it makes way more sense and might clarify the writing (e.g. no more "extinction risk was higher (i.e. Tm was lower)

library(here)
library(RColorBrewer)
setwd(here())
source("ModelFunctions.R")
source("CommonParams.R")
NumGens <- 1000
load("SimData/ExtDynamics.rdata")

#----------------- Calculate and plot the intrinsic mean time to extinction (Tm)
# Grimm, V. and C. Wissel. 2004. The intrinsic mean time to extinction: a unifying
#    approach to analysing persistence and viability of populations. Oikos 105:501-511

# Create a function to take in a vector of extinction times (and NA values) and
#    use the methodology of Grimm and Wissel to calculate Tm
CalcTm <- function(ExtTimes, BreakPoints = seq(from = 0, to = 1000, by = 25),
                   NumSims = 10000, LinearityThreshold = 0.8){
     # First, use the histogram function to bin the extinction times and then
     #    divide by the number of simulations to calculate the probability of 
     #    extinction through time
     ExtCounts <- hist(ExtTimes, breaks = BreakPoints, plot = FALSE)
     # Only use points for which we have extinctions (i.e. counts > 0)
     TimeVals <- ExtCounts$mids[ExtCounts$counts > 0]
     ExtByTime <- ExtCounts$counts[ExtCounts$counts > 0] / NumSims
     # Calculate the cumulative probability of extinction at those time periods
     #    and use it to calculate -ln(1 - P0(t)) from Grimm and Wissel (2004)
     CumulativeExt <- cumsum(ExtByTime)
     LogVal <- -1 * log(1 - CumulativeExt)
     # Remove infinite values that can occur if ther eis total extinction
     TimeVals <- TimeVals[!is.infinite(LogVal)]
     LogVal <- LogVal[!is.infinite(LogVal)]
     # Finally check for sufficient linearity and if it exists, calculate Tm
     if(length(LogVal) > 1){
          if(cor(TimeVals, LogVal) >= LinearityThreshold){
               LinMod <- lm(LogVal ~ TimeVals)
               Tm <- 1 / coef(LinMod)[2]
          } else{
               Tm <- NA
          }
     } else{
          Tm <- NA
     }
     return(Tm)
}

# Create an object to store all the intrinsic mean times to extinction
Tm <- array(NA, dim = c(3, length(TauSeq), 4))
NumSims <- 10000
for(i in 1:3){
     for(j in 1:length(TauSeq)){
          for(k in 1:4){
               CurExtTimes <- ExtTime[i,k,j,]
               if(sum(is.na(CurExtTimes)) < NumSims){
                    Tm[i,j,k] <- CalcTm(ExtTimes = CurExtTimes, LinearityThreshold = 0.5)
               }
          }
     }
}

# Now create some useful objects to plot the Tm values and create the graphic
FigName <- "Figures/Inverse_Tm.pdf"
# Sizing objects
LabelSize <- 1.15
AxisSize <- 1.25
LegendSize <- 1.15
# Location objects
yLabLine <- 2.25
xLabLine <- 2
LetterLine <- -1.25
LetterAdj <- 0.025
# Ranges
xRange <- range(TauSeq)
yRange <-c(-16, 0)
# Labels
yLabel <- expression(paste("ln(1/", T[m], ")", sep = ""))
LetterVec <- c("a", "b", "c")
# Colors
ModelColorSeq <- brewer.pal(n = 4, name = "Dark2")
ModelLineSeq <- 2:5

# Make the figure
pdf(file = FigName, width = 3, height = 5, onefile = FALSE, paper = "special")
     par(mfcol = c(3,1), oma = c(4,4.5,0.5,0.5), mar = c(0.5,0.5,0.5,0.5))
     for(i in 1:3){
          plot(x = NA, y = NA, xlim = xRange, ylim = yRange, main = "", 
               ylab = "", xlab = "", axes = FALSE)
          for(j in 1:4){
               lines(x = TauSeq, y = log(1/Tm[i,,j]), col = ModelColorSeq[j], 
                     lty = ModelLineSeq[j])
          }
          # Add a letter for each graph
          mtext(LetterVec[i], side = 3, adj = LetterAdj, line = LetterLine)
          # Add the axes to the figure
          box()
          if(i == 3){
               axis(1, cex.axis = AxisSize)
          }else{
               axis(1, labels = FALSE)
          }
          # Now add a zoomed-in inset or the legend to each graph and the y axis
          if(i == 1){
               legend("topright", legend = c("Demographic", "SexRatio", "GroupFormation", "Full"),
                      col = ModelColorSeq, lty = ModelLineSeq, bty = "n", cex = LegendSize, inset = -0.02)
               axis(2, las = 1, at = seq(-16, 0, by = 4), cex.axis = AxisSize)
               axis(2, at = seq(-16, 0, by = 1), labels = FALSE, tcl = -0.25)
          } else{
               axis(2, las = 1, at = seq(-16, 0, by = 4), cex.axis = AxisSize)
               axis(2, at = seq(-16, 0, by = 1), labels = FALSE, tcl = -0.25)
          }
     }
     # Add y and x axis labels
     mtext("Group size", side = 1, outer = TRUE, line = xLabLine, cex = LabelSize)
     mtext(yLabel, side = 2, outer = TRUE, line = yLabLine, cex = LabelSize)
dev.off()





