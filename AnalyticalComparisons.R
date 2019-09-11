# This script will check the deterministic model dynamics produced by the 
#    simulations against certain analytical constraints derived from the 
#    Ricker equation with group formation added in.
# The following lists the notation used in this script:
#    G_hat = number of groups at equilibrium
#    tau = group size
#    z = expected proportion of females (at both the population and group level)
#    R = intrinsic growth rate of the population
#    alpha = strength of density dependent competition (density calculated at population level)
#    K = number of indivduals produced each generation at equilibrium (pre group formation)=
# Equilibrium Constraint 1: 
#    G_hat <= -ln((tau*z)/(floor(z*tau)*R)) / (alpha*tau)
#    But only when: G_hat > 0, floor(z*tau) > 0, and (tau*z)/floor(tau*z) <= R
# Constraint 2:
#    K = G_hat*floor(z*tau)*(R/z)*exp(-alpha*G_hat*tau)
#    At z = 0.5, for even tau this reduces to: K = G_hat*tau*R*exp(-alpha*G_hat*tau)
#    At z = 0.5, for odd tau this reduces to: K < G_hat*tau*R*exp(-alpha*G_hat*tau)

library(schoolmath) # for functions is.even() and is.odd()

# Source the function script
source("ModelFunctions.R")
source("CommonParams.R")

# Set the parameter values for the upcoming simulations
NumGens <- 100
F_init <- 25
M_init <- 25
GroupSizeVec <- seq(2, DeterministicK, by = 1)
MateSys <- 2   # For now, just focus on mongamy because that seems most volatile
R <- R0
SurvProb <- 1

################################################################################
# First, test for constraint 1 by simulating model dynamics, calculating the
#    number of groups at equilibrium across group sizes, and plotting that against
#    the upper boundary produced by Constraint 1.
G_hat <- rep(NA, length(GroupSizeVec))
UprBound <- rep(NA, length(GroupSizeVec))
for(i in 1:length(GroupSizeVec)){
     tau <- GroupSizeVec[i]
     AbundVals <- Deterministic(ms = MateSys, R = R, alpha = alpha, z = z, 
                                GroupSize = tau, SurvProb = SurvProb, 
                                NumGens = NumGens, F_init = F_init, M_init = M_init,
                                ReturnN = TRUE)
     Gvals <- floor(AbundVals[(NumGens - 20):NumGens] / tau)
     G_hat[i] <- mean(Gvals)
     # This upper bound is only valid under the following condition
     if( ((tau*z) / floor(tau*z)) <= R ){
          UprBound[i] <- log( (floor(z*tau)*R) / (z * tau) ) / (alpha * tau)
     }
}

yRange <- c(0, 75)
AxisSize <- 1.5
LabelSize <- 1.5
LegendSize <- 1.5
pdf(file = "G_hat.pdf", width = 8, height = 5, onefile = FALSE, paper = "special")
     plot(x = NA, y = NA, xlim = range(GroupSizeVec), ylim = yRange, main = "",
          xlab = "Group size", ylab = "Number of groups", las = 1, cex.axis = AxisSize,
          cex.lab = LabelSize)
     axis(1, at = seq(2, DeterministicK, by = 2), tcl = -0.25, labels = FALSE)
     lines(x = GroupSizeVec[is.even(GroupSizeVec)], y = G_hat[is.even(GroupSizeVec)], col = "darkgreen")
     lines(x = GroupSizeVec[is.even(GroupSizeVec)], y = UprBound[is.even(GroupSizeVec)], lty = 2, col = "darkgreen")
     lines(x = GroupSizeVec[is.odd(GroupSizeVec)], y = G_hat[is.odd(GroupSizeVec)], col = "red")
     lines(x = GroupSizeVec[is.odd(GroupSizeVec)], y = UprBound[is.odd(GroupSizeVec)], lty = 2, col = "red")
     legend("topright", legend = c("Even", "Odd", "Analytical upper bound"), lty = c(1, 1, 2),
            col = c("darkgreen", "red", "black"), bty = "n", cex = LegendSize)
dev.off()

################################################################################
# Now, test for constraint 2 by simulating model dynamics, calculating the
#    population size pre group formation at equilibrium across group sizes, 
#    and plotting that against the expectations set by constraint 2.
Kvals <- rep(NA, length(GroupSizeVec))
Expectation <- rep(NA, length(GroupSizeVec))
for(i in 1:length(GroupSizeVec)){
     tau <- GroupSizeVec[i]
     AbundVals <- Deterministic(ms = MateSys, R = R, alpha = alpha, z = z, 
                                GroupSize = tau, SurvProb = SurvProb, 
                                NumGens = NumGens, F_init = F_init, M_init = M_init,
                                ReturnN = TRUE)
     Kvals[i] <- mean(AbundVals[(NumGens - 20):NumGens])
     Expectation[i] <- G_hat[i] * tau * R * exp(-1 * alpha * G_hat[i] * tau)
}

yRange <- range(c(Kvals, Expectation))
AxisSize <- 1.5
LabelSize <- 1.5
LegendSize <- 1.5
pdf(file = "K.pdf", width = 8, height = 5, onefile = FALSE, paper = "special")
     par(mar = c(5, 5, 4, 2) + 0.1)
     plot(x = NA, y = NA, xlim = range(GroupSizeVec), ylim = yRange, main = "",
          xlab = "Group size", ylab = "Equilibrium abundance", las = 1, cex.axis = AxisSize,
          cex.lab = LabelSize)
     axis(1, at = seq(2, DeterministicK, by = 2), tcl = -0.25, labels = FALSE)
     lines(x = GroupSizeVec[is.even(GroupSizeVec)], y = Kvals[is.even(GroupSizeVec)], col = "darkgreen")
     points(x = GroupSizeVec[is.even(GroupSizeVec)], y = Expectation[is.even(GroupSizeVec)], col = "darkgreen")
     lines(x = GroupSizeVec[is.odd(GroupSizeVec)], y = Kvals[is.odd(GroupSizeVec)], col = "red")
     lines(x = GroupSizeVec[is.odd(GroupSizeVec)], y = Expectation[is.odd(GroupSizeVec)], lty = 2, col = "red")
     abline(h = DeterministicK, lty = 3)
     legend("bottomright", legend = c("Even", "Odd", "Analytical expectation", "Analytical upper bound"), 
            lty = c(1, 1, 0, 2), pch = c(NA, NA, 1, NA), col = c("darkgreen", "red", "darkgreen", "red"), 
            bty = "n", cex = LegendSize)
dev.off()


# Make a 3 panel graph with everything put together
AxisSize <- 1.5
LabelSize <- 1.25
LegendSize <- 1.25
pdf(file = "AnalyticalExpectations.pdf", width = 8, height = 6, onefile = FALSE, paper = "special")
     par(mar = c(2, 6, 2, 2) + 0.1, mfrow = c(3,1), oma = c(2, 0, 2, 2))
     # First plot the number of groups with the expectations
     yRange <- c(0, 75)
     plot(x = NA, y = NA, xlim = range(GroupSizeVec), ylim = yRange, main = "",
          xlab = "Group size", ylab = "Number of groups", las = 1, cex.axis = AxisSize,
          cex.lab = LabelSize)
     axis(1, at = seq(2, DeterministicK, by = 2), tcl = -0.25, labels = FALSE)
     lines(x = GroupSizeVec[is.even(GroupSizeVec)], y = G_hat[is.even(GroupSizeVec)], col = "darkgreen")
     lines(x = GroupSizeVec[is.even(GroupSizeVec)], y = UprBound[is.even(GroupSizeVec)], lty = 2, col = "darkgreen")
     lines(x = GroupSizeVec[is.odd(GroupSizeVec)], y = G_hat[is.odd(GroupSizeVec)], col = "red")
     lines(x = GroupSizeVec[is.odd(GroupSizeVec)], y = UprBound[is.odd(GroupSizeVec)], lty = 2, col = "red")
     legend("topright", legend = c("Even", "Odd", "Analytical upper bound"), lty = c(1, 1, 2),
            col = c("darkgreen", "red", "black"), bty = "n", cex = LegendSize)
     abline(v = GroupSizeVec[25], lty = 3, col = "darkgreen")
     abline(v = GroupSizeVec[24], lty = 3, col = "red")
     abline(v = GroupSizeVec[16], lty = 3, col = "red")
     # Next the abundances at equilibrium with the expectations
     yRange <- c(0, DeterministicK)
     plot(x = NA, y = NA, xlim = range(GroupSizeVec), ylim = yRange, main = "",
          xlab = "Group size", ylab = "Equilibrium abundance", las = 1, cex.axis = AxisSize,
          cex.lab = LabelSize)
     axis(1, at = seq(2, DeterministicK, by = 2), tcl = -0.25, labels = FALSE)
     lines(x = GroupSizeVec[is.even(GroupSizeVec)], y = Kvals[is.even(GroupSizeVec)], col = "darkgreen")
     points(x = GroupSizeVec[is.even(GroupSizeVec)], y = Expectation[is.even(GroupSizeVec)], col = "darkgreen")
     lines(x = GroupSizeVec[is.odd(GroupSizeVec)], y = Kvals[is.odd(GroupSizeVec)], col = "red")
     lines(x = GroupSizeVec[is.odd(GroupSizeVec)], y = Expectation[is.odd(GroupSizeVec)], lty = 2, col = "red")
     abline(h = DeterministicK, lty = 3)
     abline(v = GroupSizeVec[25], lty = 3, col = "darkgreen")
     abline(v = GroupSizeVec[24], lty = 3, col = "red")
     abline(v = GroupSizeVec[16], lty = 3, col = "red")
     legend("bottomright", legend = c("Even", "Odd", "Analytical expectation", "Analytical upper bound"), 
            lty = c(1, 1, 0, 2), pch = c(NA, NA, 1, NA), col = c("darkgreen", "red", "darkgreen", "red"), 
            bty = "n", cex = LegendSize)
     # Finally, the relationsip between tau*z and floor(tau * z)
     EvenOddFactor <- (GroupSizeVec * z) / floor(GroupSizeVec * z)
     LabExpression <- expression(over(tau*z, (. %down% tau * z %down% .)))
     plot(NA, NA, xlim = range(GroupSizeVec), ylim = c(1, 1.5), main = "", xlab = "Group size",
          ylab = LabExpression, las = 1, cex.axis = AxisSize, cex.lab = LabelSize)
     axis(1, at = seq(2, DeterministicK, by = 2), tcl = -0.25, labels = FALSE)
     lines(x = GroupSizeVec[is.even(GroupSizeVec)], y = EvenOddFactor[is.even(GroupSizeVec)],
           col = "darkgreen")
     lines(x = GroupSizeVec[is.odd(GroupSizeVec)], y = EvenOddFactor[is.odd(GroupSizeVec)],
           col = "red")
     abline(v = GroupSizeVec[25], lty = 3, col = "darkgreen")
     abline(v = GroupSizeVec[24], lty = 3, col = "red")
     abline(v = GroupSizeVec[16], lty = 3, col = "red")
     legend("topright", legend = c("Even", "Odd"), lty = 1, col = c("darkgreen", "red"),
            bty = "n", cex = LegendSize)
dev.off()

