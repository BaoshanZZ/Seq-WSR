rm(list = ls())

library(xtable)
library(tidyr)
library(ggplot2)
library(reshape2)

source("FUNC_PMF.R")
source("FUNC_SeqWSR.R")
source("FUNC_Sim.R")

## Baseline Distribution Setting
pmf.X0m4 <- list(
  Increasing = c(0.05, 0.10, 0.20, 0.25, 0.40), 
  Uniform = rep(0.20, 5),
  Decreasing = c(0.40, 0.25, 0.20, 0.10, 0.05),
  U_shaped = c(0.30, 0.15, 0.10, 0.15, 0.30),
  Bell_shaped = c(0.1, 0.20, 0.40, 0.20, 0.1)
)
pmf.X0m5 <- list(
  Increasing = 1:6/21,
  Uniform = rep(1/6, 6),
  Decreasing = 6:1/21,
  U_shaped = c(0.3, 0.15, 0.05, 0.05, 0.15, 0.3),
  Bell_shaped = c(0.05, 0.15, 0.3, 0.3, 0.15, 0.05)
)
pmf.X0m6 <- list(
  Increasing = c(0.05, 0.075, 0.1, 0.125, 0.175, 0.225, 0.25),
  Uniform = rep(1/7, 7),
  Decreasing = c(0.25, 0.225, 0.175, 0.125, 0.1, 0.075, 0.05),
  U_shaped = c(0.25, 0.15, 0.075, 0.05, 0.075, 0.15, 0.25),
  Bell_shaped = c(0.05, 0.1, 0.15, 0.4, 0.15, 0.1, 0.05)
)

OR_values <- c(1, 1.5, 2)

pmf.X1m4 <- process_results_X1(pmf.X0m4, OR_values) 
pmf.Xm4 <- process_results_X(pmf.X0m4, OR_values) 

pmf.X1m5 <- process_results_X1(pmf.X0m5, OR_values) 
pmf.Xm5 <- process_results_X(pmf.X0m5, OR_values) 

pmf.X1m6 <- process_results_X1(pmf.X0m6, OR_values) 
pmf.Xm6 <- process_results_X(pmf.X0m6, OR_values) 

plot_pmf(pmf.Xm4)
plot_pmf(pmf.X1m4)
plot_pmf(pmf.Xm5)
plot_pmf(pmf.X1m5)
plot_pmf(pmf.Xm6)
plot_pmf(pmf.X1m6)



##### Simulation Part Begins
B = 50000
alpha = 0.05
N_values = c(60, 90, 120, 150)
### QUESITION
set.seed(0803) # WHERE SHOULD I SET THE SEEDS?
K <- 3

#### m = 4
SimResult_OR1m4 <- sim_by_OR(OddsR = 1, B = B, results_X = pmf.Xm4, 
                             K = K ,alpha = alpha, N_values = N_values)
SimResult_OR15m4 <- sim_by_OR(OddsR = 1.5, B = B, results_X = pmf.Xm4, 
                              K = K ,alpha = alpha, N_values = N_values)
SimResult_OR2m4 <- sim_by_OR(OddsR = 2, B = B, results_X = pmf.Xm4, 
                             K = K ,alpha = alpha, N_values = N_values)

#### m = 5
SimResult_OR1m5 <- sim_by_OR(OddsR = 1, B = B, results_X = pmf.Xm5, 
                             K = K ,alpha = alpha, N_values = N_values)
SimResult_OR15m5 <- sim_by_OR(OddsR = 1.5, B = B, results_X = pmf.Xm5, 
                              K = K ,alpha = alpha, N_values = N_values)
SimResult_OR2m5 <- sim_by_OR(OddsR = 2, B = B, results_X = pmf.Xm5, 
                             K = K ,alpha = alpha, N_values = N_values)

#### m = 6
SimResult_OR1m6 <- sim_by_OR(OddsR = 1, B = B, results_X = pmf.Xm6, 
                             K = K ,alpha = alpha, N_values = N_values)

SimResult_OR15m6 <- sim_by_OR(OddsR = 1.5, B = B, results_X = pmf.Xm6, 
                              K = K ,alpha = alpha, N_values = N_values)
SimResult_OR2m6 <- sim_by_OR(OddsR = 2, B = B, results_X = pmf.Xm6, 
                             K = K ,alpha = alpha, N_values = N_values)

###### Report Results simulation
#### OR = 1
SimPlot_AttainedAlpha(SimResult_OR1m4)
SimPlot_BiasRatio(SimResult_OR1m4)

SimPlot_BiasRatio(SimResult_OR1m5)
SimPlot_AttainedAlpha(SimResult_OR1m5)

SimPlot_BiasRatio(SimResult_OR1m6)
SimPlot_AttainedAlpha(SimResult_OR1m6)

### OR = 1.5
SimPlot_BiasRatio(SimResult_OR15m4)
SimPlot_AttainedPower(SimResult_OR15m4)

SimPlot_BiasRatio(SimResult_OR15m5)
SimPlot_AttainedPower(SimResult_OR15m5)

SimPlot_BiasRatio(SimResult_OR15m6)
SimPlot_AttainedPower(SimResult_OR15m6)


#### OR = 2
SimPlot_BiasRatio(SimResult_OR2m4)
SimPlot_AttainedPower(SimResult_OR2m4)

SimPlot_BiasRatio(SimResult_OR2m5)
SimPlot_AttainedPower(SimResult_OR2m5)

SimPlot_BiasRatio(SimResult_OR2m6)
SimPlot_AttainedPower(SimResult_OR2m6)

