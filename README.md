# SocialStochasticity

This file contains a brief overview of the model described in the manuscript "Stochasticity in social structure and mating system drive extinction risk" by Damon Leach, Allison Shaw, and Christopher Weiss-Lehman. In addition to describing the general structure of the model, this file contains descriptions of the R scripts contained within the repository and their purpose in the overall model.

The model implemented in this repository builds on the classic Ricker model (Ricker 1954) by adding social structure in the form of social groups and different mating systems. Additionally, the model selectively implements different combinations of stochasticity to assess the impact of different forms of stochasticity on extinction risk in social species. All together, the model implements three different mating systems (polygynandry, monogamy, and alpha pair monogamy in which only a single female reproduces per group). Stochasticity in the model arises from stochasticity in the number of offspring produced (demographic stochasticity), in the sex ratio of the population (stochastic sex determination), and in the group formation process (stochastic group formation). This means that for each mating system, we examine five different models over a range of group sizes: (1) a deterministic model, (2) a model with only demographic stochasticity, (3) a model with demographic stochasticity and stochastic sex determination, (4) a model with demographic stochasticity and stochastic group formation, and (5) the full model with all sources of stochasticity. Using these model suites, the scripts in this repository calculate the intrinsic mean time to extinction (Grimm and Wissel 2004) for all models across a range of group sizes, as well as the effect of stochastic group formation and stochastic sex determination specifically on expected population growth rates.

Below is a brief description of each of the scripts in this repository. Additionally, each script is itself thoroughly documented to promote understanding and reuse of the model implemented here.

-- ModelFunctions.R
This script contains all of the individual functions necessary to run the model simulations. In addition to a function for each of the 5 models incorporating different sources of stochasticity listed above, there are functions to perform the group formation process and determine the number of mated females in each group according to different mating systems. All functions are documented with descriptions of the necessary inputs and the resulting outputs.

--CommonParams.R
This script simply contains some common parameter values used across all subsequent simulations (e.g. the density-independent mean per capita growth rate for the Ricker model). This is included for convenience so that exploring other regions of parameter space can be accomplished by simply changing the values in this script, which is then subsequently run at the beginning of the following simulation scripts.

--ExampleAbundCurves.R
This script simply runs a single realization of each stochastic model for a given mating system and group size, producing a graph to illustrate the population fluctuations experienced in each model through time. This is meant simply as a useful visualization of the model dynamics and not as rigorous analysis since it only performs a single realization at a time. The following scripts perform the anaylese of the model described by Leach et al. (in revision). 

--ExtDynamics.R & ExtDynamics.sh
These scripts use the package Rmpi to run parralel simulations of the model on an external server to thoroughly characterize the extinction dynamics induced by the different combinations of mating systems, group sizes, and sources of stochasticity. The scripts perform 1000 simulations of each model scenario for 1000 generations to track whether and when each realization goes extinct. These data are then used in subsequent scripts to calculate and graph the intrinsic mean time to extinction for each scenario.

--ExtDynamicsPlots.R
This script processes the output from ExtDynamics.R to caculate and plot the intrinsic mean time to extinction. The calculation of the intrinsic mean time to extinction follows the procedure laid out in Grimm and Wissel (2004) for stochastic, simulation models.

--EffectOfGroupFormation.R
This script quantifies the effect of stochastic group formation by comparing expected population growth under deterministic versus stochastic group formation in populations at three different sex ratios. It runs 10000 simulations for each group size and mating system investigated in the manuscript, plotting the average effects of stochastic group formation (negative values indicate reduced growth compared to deterministic group formation) with error bars representing the standard deviations.

--EffectOfStochasticSexRatio.R
This script is similar to the above script, but it investigates and plots the effects of stochastic sex determination on expected population growth rates. It also performs 10000 simulations for each group size and mating system, comparing the expected population growth rates when sex determination is deterministic (i.e. exactly matching a given sex ratio) or stochastic. It also plots the results as the average effect with error bars representing the standard deviations.
