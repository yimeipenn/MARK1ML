
*## Libraries that are used
  library(xtable)
  library(MASS)
  library(geeM)


## Load the code for the analysis:
  source("h://Bender//20160610_Project2_Code.R")

## Read the data:
data2 = read.csv("h:/Bender/OneSimulatedData.csv")

## Run the analysis:

## Analysis of epilepsy data with Poisson and AR(1) assumptions
   PoisAR = EndResults(y ~ trt + base + age + period, "AR(1)", "Poisson", data2, data2$subject, data2$period, rep(0,6))
   PoisAR

## Analysis of epilepsy data with Poisson and Markov assumptions
  PoisMark = EndResults(y ~ trt + base + age + period, "Markov", "Poisson", data2, data2$subject, data2$period, rep(0,6))
  PoisMark