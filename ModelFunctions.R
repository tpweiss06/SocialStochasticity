# This script contains all the function definitions necessary to recreate the
#    models examining the effect of different stochastic processes on population
#    dynamics of social species under different mating systems

# ---------------------------------------------------------------- MatingFemales
# This function determines the number of mating females in a population
#    according to the mating system and group structure.
# Arguments:
#    ms -- The mating system of the population
#    Ng -- The number of groups in the population
#    GroupFemales/GroupMales -- Vectors with the number of females/males per group
MatingFemales <- function(ms, Ng, GroupFemales, GroupMales){
     MatingFemales <- rep(NA, Ng)
     # Loop through each group to calculate the number of mating females in each
     #    group
     for(i in 1:Ng){
          if(ms == 1){
               MatingFemales[i] <- ifelse((GroupFemales[i] > 0) & (GroupMales[i] > 0),
                                               1, 0)
          } else if(ms == 2){
               MatingFemales[i] <- min(GroupFemales[i], GroupMales[i])
          } else{
               MatingFemales[i] <- ifelse(GroupMales[i] > 0, GroupFemales[i], 0)
          }
     }
     return(sum(MatingFemales))
}

# ------------------------------------------------------------------- FormGroups
# This function forms groups of a given size from a full population consisting
#    of males and females. It works either by forming each group deterministically
#    according to the given sex ratio or by unbiased random sampling.
# Arguments:
#    Tau -- the size of groups to be formed
#    Females and Males -- The respective number of each in the total population
#    Stochastic -- A boolean indicator for whether group formation should be purely
#         stochastic (i.e. via unbiased random sampling) or biased towards forming
#         groups with a given sex ratio
#    z -- The sex ratio to work towards if Stochastic is set to FALSE
FormGroups <- function(Tau, Females, Males, Stochastic, z = NA){
     # Set the total population size
     N <- Females + Males
     if(Stochastic){
          # Determine the number of groups the population can form
          Ng <- floor(N / Tau)
          
          # Create vectors to track the number of females and males per group
          GroupFemales <- rep(0, Ng)
          GroupMales <- rep(0, Ng)
          
          # Now step through each group and randomly pull males and females from
          #    the full population to fill each group
          UngroupedFems <- Females
          UngroupedMales <- Males
          UngroupedPop <- c(rep(1, UngroupedFems), rep(0, UngroupedMales))
          for(i in 1:Ng){
               CurGroup <- sample(UngroupedPop, size = Tau, replace = FALSE)
               CurFemales <- sum(CurGroup)
               CurMales <- Tau - CurFemales
               GroupFemales[i] <- GroupFemales[i] + CurFemales
               GroupMales[i] <- GroupMales[i] + CurMales
               UngroupedFems <- UngroupedFems - CurFemales
               UngroupedMales <- UngroupedMales - CurMales
               UngroupedPop <- c(rep(1, UngroupedFems), rep(0, UngroupedMales))
          }
     }else{
          # Determine the expected number of females and males per group
          ExpectedFem <- floor(Tau * z)
          ExpectedMale <- Tau - ExpectedFem
          
          # Calculate the number of groups that can form from the female and
          #    male populations separately
          FemGroups <- floor(Females / ExpectedFem)
          MaleGroups <- floor(Males / ExpectedMale)
          
          # Set the total number of groups to be formed at exactly the exptected
          #    sex ratios based on the minimum between the two quantities
          Ng1 <- min(FemGroups, MaleGroups)
          
          # Calculate the remaining numbers of females and males after forming
          #    Ng1 groups at the expected sex ratio
          FemRemaining <- Females - Ng1 * ExpectedFem
          MaleRemaining <- Males - Ng1 * ExpectedMale
          
          # Now check to see if there are enough individuals to form one final 
          #    group with the remainders, regardless of sex ratio. All other 
          #    remaining individuals are disgarded after this step even if they
          #    could form a group since that group would be composed of only one
          #    sex and hence have no reproductive success.
          if( (FemRemaining > 0) & (MaleRemaining > 0) & ((FemRemaining + MaleRemaining) > Tau) ){
               Ng <- Ng1 + 1
               if(FemGroups > MaleGroups){
                    MaleFinal <- MaleRemaining
                    FemFinal <- Tau - MaleFinal
               } else{
                    FemFinal <- FemRemaining
                    MaleFinal <- Tau - FemFinal
               }
          } else{
               Ng <- Ng1
          }
          
          # Finally, make the vectors with the numbers of males and females per
          #    group
          GroupFemales <- rep(0, Ng)
          GroupMales <- rep(0, Ng)
          GroupFemales[1:Ng1] <- ExpectedFem
          GroupMales[1:Ng1] <- ExpectedMale
          if(Ng > Ng1){
               GroupFemales[Ng] <- FemFinal
               GroupMales[Ng] <- MaleFinal
          }
     }
     GroupedPop <- rbind(GroupFemales, GroupMales)
     return(GroupedPop)
}

# ---------------------------------------------------------------- Deterministic
# This function will model the population dynamics of populations with various
#    mating systems and group sizes fully deterministically.
# Arguments:
#    ms -- The mating system of the population
#    R -- The intrinsic growth rate of the population
#    alpha -- The strength of negative density dependence
#    z -- The expected sex ratio of the population
#    Tau -- The number of individuals in each group of the population
#    NumGens -- The number of generations to simulate
#    F_init -- The initial number of females
#    M_init -- The initial number of males
#    ReturnN -- A flag to indicate whether to return just the vector of population
#              sizes (TRUE) or the vectors of male and female abundances (FALSE)
Deterministic <- function(ms, R, alpha, z, Tau, NumGens, F_init, M_init, ReturnN){
     # Initialize population vectors
     N <- rep(0, NumGens)
     Fem <- rep(0, NumGens)
     M <- rep(0, NumGens)
     N[1] <- F_init + M_init
     Fem[1] <- F_init
     M[1] <- M_init
     
     # Create the generation counter and an extinction flag
     t <- 1
     Extinct <- FALSE
     
     # Now loop through the generations
     while((t < NumGens) & (!Extinct)){
          # First form groups "deterministically"
          GroupedPop <- FormGroups(Tau = Tau, Females = Fem[t],
                                   Males = M[t], Stochastic = FALSE, z = z)
          Ng <- dim(GroupedPop)[2]
          GroupFemales <- GroupedPop[1,]
          GroupMales <- GroupedPop[2,]
          AllPop <- sum(GroupFemales) + sum(GroupMales)
          
          # Determine the number of mating females in the population
          MatedFems <- MatingFemales(ms = ms, Ng = Ng, GroupMales = GroupMales, 
                                     GroupFemales = GroupFemales)
          
          # Calculate the expected population size in the next generation and
          #    use it to draw the realized population size from a Poisson 
          #    distribution
          ExpectedPopSize <- MatedFems * (R/z) * exp(-1 * alpha * AllPop)
          N[t+1] <- round(ExpectedPopSize)
          
          # Determine the stochastic sex ratio of next generation's population
          Fem[t+1] <- round(N[t+1] * z)
          M[t+1] <- N[t+1] - Fem[t+1]
          
          # Finally, check for extinction conditions and advance t
          if((N[t+1] < Tau) | (Fem[t+1] == 0) | (M[t+1] == 0)){
               Extinct <- TRUE
          }
          t <- t + 1
     }
     if(ReturnN){
          return(N)
     } else{
          PopMatrix <- cbind(Fem, M)
          return(PopMatrix)
     }
}

# ------------------------------------------------------------------ Demographic
# This function will model the population dynamics of populations with various
#    mating systems and group sizes, incorporating only demographic stochasticity.
# Arguments:
#    ms -- The mating system of the population
#    R -- The intrinsic growth rate of the population
#    alpha -- The strength of negative density dependence
#    z -- The expected sex ratio of the population
#    Tau -- The number of individuals in each group of the population
#    NumGens -- The number of generations to simulate
#    F_init -- The initial number of females
#    M_init -- The initial number of males
#    ReturnN -- A flag to indicate whether to return just the vector of population
#              sizes (TRUE) or the vectors of male and female abundances (FALSE)
Demographic <- function(ms, R, alpha, z, Tau, NumGens, F_init, M_init, ReturnN){
     # Initialize population vectors
     N <- rep(0, NumGens)
     Fem <- rep(0, NumGens)
     M <- rep(0, NumGens)
     N[1] <- F_init + M_init
     Fem[1] <- F_init
     M[1] <- M_init
     
     # Create the generation counter and an extinction flag
     t <- 1
     Extinct <- FALSE
     
     # Now loop through the generations
     while((t < NumGens) & (!Extinct)){
          # First form groups "deterministically"
          GroupedPop <- FormGroups(Tau = Tau, Females = Fem[t], Males = M[t], 
                                   Stochastic = FALSE, z = z)
          Ng <- dim(GroupedPop)[2]
          GroupFemales <- GroupedPop[1,]
          GroupMales <- GroupedPop[2,]
          AllPop <- sum(GroupFemales) + sum(GroupMales)
          
          # Determine the number of mating females in the population
          MatedFems <- MatingFemales(ms = ms, Ng = Ng, GroupMales = GroupMales, 
                                     GroupFemales = GroupFemales)
          
          # Calculate the expected population size in the next generation and
          #    use it to draw the realized population size from a Poisson 
          #    distribution
          ExpectedPopSize <- MatedFems * (R/z) * exp(-1 * alpha * AllPop)
          N[t+1] <- rpois(n = 1, lambda = ExpectedPopSize)
          
          # Determine the stochastic sex ratio of next generation's population
          Fem[t+1] <- round(N[t+1] * z)
          M[t+1] <- N[t+1] - Fem[t+1]
          
          # Finally, check for extinction conditions and advance t
          if((N[t+1] < Tau) | (Fem[t+1] == 0) | (M[t+1] == 0)){
               Extinct <- TRUE
          }
          t <- t + 1
     }
     if(ReturnN){
          return(N)
     } else{
          PopMatrix <- cbind(Fem, M)
          return(PopMatrix)
     }
}

# --------------------------------------------------------------------- SexRatio
# This function will model the population dynamics of populations with various
#    mating systems and group sizes, incorporating demographic stochasticity and
#    stochasticity in sex ratio.
# Arguments:
#    ms -- The mating system of the population
#    R -- The intrinsic growth rate of the population
#    alpha -- The strength of negative density dependence
#    z -- The expected sex ratio of the population
#    Tau -- The number of individuals in each group of the population
#    NumGens -- The number of generations to simulate
#    F_init -- The initial number of females
#    M_init -- The initial number of males
#    ReturnN -- A flag to indicate whether to return just the vector of population
#              sizes (TRUE) or the vectors of male and female abundances (FALSE)
SexRatio <- function(ms, R, alpha, z, Tau, NumGens, F_init, M_init, ReturnN){
     # Initialize population vectors
     N <- rep(0, NumGens)
     Fem <- rep(0, NumGens)
     M <- rep(0, NumGens)
     N[1] <- F_init + M_init
     Fem[1] <- F_init
     M[1] <- M_init
     
     # Create the generation counter and an extinction flag
     t <- 1
     Extinct <- FALSE
     
     # Now loop through the generations
     while((t < NumGens) & (!Extinct)){
          # First form groups "deterministically"
          GroupedPop <- FormGroups(Tau = Tau, Females = Fem[t], Males = M[t], 
                                   Stochastic = FALSE, z = z)
          Ng <- dim(GroupedPop)[2]
          GroupFemales <- GroupedPop[1,]
          GroupMales <- GroupedPop[2,]
          AllPop <- sum(GroupFemales) + sum(GroupMales)
          
          # Determine the number of mating females in the population
          MatedFems <- MatingFemales(ms = ms, Ng = Ng, GroupMales = GroupMales, 
                                     GroupFemales = GroupFemales)
          
          # Calculate the expected population size in the next generation and
          #    use it to draw the realized population size from a Poisson 
          #    distribution
          ExpectedPopSize <- MatedFems * (R/z) * exp(-1 * alpha * AllPop)
          N[t+1] <- rpois(n = 1, lambda = ExpectedPopSize)
          
          # Determine the stochastic sex ratio of next generation's population
          NewFemales <- rbinom(n = N[t+1], size = 1, prob = z)
          Fem[t+1] <- sum(NewFemales)
          M[t+1] <- N[t+1] - Fem[t+1]
          
          # Finally, check for extinction conditions and advance t
          if((N[t+1] < Tau) | (Fem[t+1] == 0) | (M[t+1] == 0)){
               Extinct <- TRUE
          }
          t <- t + 1
     }
     if(ReturnN){
          return(N)
     } else{
          PopMatrix <- cbind(Fem, M)
          return(PopMatrix)
     }
}

# --------------------------------------------------------------- GroupFormation
# This function will model the population dynamics of populations with various
#    mating systems and group sizes, incorporating demographic stochasticity and
#    stochasticity in group formation, but a deterministic sex ratio.
# Arguments:
#    All arguments are the same as in the SexRatio model
GroupFormation <- function(ms, R, alpha, z, Tau, NumGens, F_init, M_init, ReturnN){
     # Initialize population vectors
     N <- rep(0, NumGens)
     Fem <- rep(0, NumGens)
     M <- rep(0, NumGens)
     N[1] <- F_init + M_init
     Fem[1] <- F_init
     M[1] <- M_init
     
     # Create the generation counter and an extinction flag
     t <- 1
     Extinct <- FALSE
     
     # Now loop through the generations
     while((t < NumGens) & (!Extinct)){
          # First form groups stochasticly
          GroupedPop <- FormGroups(Tau = Tau, Females = Fem[t], Males = M[t], 
                                   Stochastic = TRUE)
          Ng <- dim(GroupedPop)[2]
          GroupFemales <- GroupedPop[1,]
          GroupMales <- GroupedPop[2,]
          AllPop <- sum(GroupFemales) + sum(GroupMales)
          
          # Determine the number of mating females in the population
          MatedFems <- MatingFemales(ms = ms, Ng = Ng, GroupMales = GroupMales, 
                                     GroupFemales = GroupFemales)
          
          # Calculate the expected population size in the next generation and
          #    use it to draw the realized population size from a Poisson 
          #    distribution
          ExpectedPopSize <- MatedFems * (R/z) * exp(-1 * alpha * AllPop)
          N[t+1] <- rpois(n = 1, lambda = ExpectedPopSize)
          
          # Determine the sex ratio of next generation's population
          Fem[t+1] <- round(N[t+1] * z)
          M[t+1] <- N[t+1] - Fem[t+1]
          
          # Finally, check for extinction conditions and advance t
          if((N[t+1] < Tau) | (Fem[t+1] == 0) | (M[t+1] == 0)){
               Extinct <- TRUE
          }
          t <- t + 1
     }
     if(ReturnN){
          return(N)
     } else{
          PopMatrix <- cbind(Fem, M)
          return(PopMatrix)
     }
}

# ------------------------------------------------------------------------- Full
# This function will model the population dynamics of populations with various
#    mating systems and group sizes, incorporating demographic stochasticity,
#    stochasticity in group formation, and a stochastic sex ratio.
# Arguments:
#    All arguments are the same as in the SexRatio model
Full <- function(ms, R, alpha, z, Tau, NumGens, F_init, M_init, ReturnN){
     # Initialize population vectors
     N <- rep(0, NumGens)
     Fem <- rep(0, NumGens)
     M <- rep(0, NumGens)
     N[1] <- F_init + M_init
     Fem[1] <- F_init
     M[1] <- M_init
     
     # Create the generation counter and an extinction flag
     t <- 1
     Extinct <- FALSE
     
     # Now loop through the generations
     while((t < NumGens) & (!Extinct)){
          # First form groups stochasticly
          GroupedPop <- FormGroups(Tau = Tau, Females = Fem[t], Males = M[t], 
                                   Stochastic = TRUE)
          Ng <- dim(GroupedPop)[2]
          GroupFemales <- GroupedPop[1,]
          GroupMales <- GroupedPop[2,]
          AllPop <- sum(GroupFemales) + sum(GroupMales)
          
          # Determine the number of mating females in the population
          MatedFems <- MatingFemales(ms = ms, Ng = Ng, GroupMales = GroupMales, 
                                     GroupFemales = GroupFemales)
          
          # Calculate the expected population size in the next generation and
          #    use it to draw the realized population size from a Poisson 
          #    distribution
          ExpectedPopSize <- MatedFems * (R/z) * exp(-1 * alpha * AllPop)
          N[t+1] <- rpois(n = 1, lambda = ExpectedPopSize)
          
          # Determine the sex ratio of next generation's population
          NewFemales <- rbinom(n = N[t+1], size = 1, prob = z)
          Fem[t+1] <- sum(NewFemales)
          M[t+1] <- N[t+1] - Fem[t+1]
          
          # Finally, check for extinction conditions and advance t
          if((N[t+1] < Tau) | (Fem[t+1] == 0) | (M[t+1] == 0)){
               Extinct <- TRUE
          }
          t <- t + 1
     }
     if(ReturnN){
          return(N)
     } else{
          PopMatrix <- cbind(Fem, M)
          return(PopMatrix)
     }
}
