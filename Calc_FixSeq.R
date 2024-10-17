rm(list = ls())

library(xtable)
library(tidyr)
library(ggplot2)
library(reshape2)
library(dplyr)
source("FUNC_PMF.R")
source("FUNC_SeqWSR.R")
source("FUNC_Sim.R")

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


alpha1 = 0.05
K = 3
N_values = c(60, 90, 120, 150)


### OR = 2 & m=4
pmf.Xm4OR2 <- pmf.Xm4 %>% filter(OR == 2)
power_tk(alpha1, K, pmf.Xm4OR2[1,3:11], N = 126)
power_tk(alpha1, 1, pmf.Xm4OR2[1,3:11], N = 123)

power_tk(alpha1, K, pmf.Xm4OR2[2,3:11], N = 106)
power_tk(alpha1, 1, pmf.Xm4OR2[2,3:11], N = 103)

power_tk(alpha1, K, pmf.Xm4OR2[3,3:11], N = 113)
power_tk(alpha1, 1, pmf.Xm4OR2[3,3:11], N = 110)

power_tk(alpha1, K, pmf.Xm4OR2[4,3:11], N = 94)
power_tk(alpha1, 1, pmf.Xm4OR2[4,3:11], N = 91)

power_tk(alpha1, K, pmf.Xm4OR2[5,3:11], N = 117)
power_tk(alpha1, 1, pmf.Xm4OR2[5,3:11], N = 114)

### OR = 2 & m = 5
pmf.Xm5OR2 <- pmf.Xm5 %>% filter(OR == 2)
power_tk(alpha1, K, pmf.Xm5OR2[1,3:13], N = 118)
power_tk(alpha1, 1, pmf.Xm5OR2[1,3:13], N = 115)

power_tk(alpha1, K, pmf.Xm5OR2[2,3:13], N = 107)
power_tk(alpha1, 1, pmf.Xm5OR2[2,3:13], N = 104)

power_tk(alpha1, K, pmf.Xm5OR2[3,3:13], N = 110)
power_tk(alpha1, 1, pmf.Xm5OR2[3,3:13], N = 107)

power_tk(alpha1, K, pmf.Xm5OR2[4,3:13], N = 95)
power_tk(alpha1, 1, pmf.Xm5OR2[4,3:13], N = 92)

power_tk(alpha1, K, pmf.Xm5OR2[5,3:13], N = 115)
power_tk(alpha1, 1, pmf.Xm5OR2[5,3:13], N = 112)

### OR = 2 & m = 6
pmf.Xm6OR2 <- pmf.Xm6 %>% filter(OR == 2)
power_tk(alpha1, K, pmf.Xm6OR2[1,3:15], N = 117)
power_tk(alpha1, 1, pmf.Xm6OR2[1,3:15], N = 114)

power_tk(alpha1, K, pmf.Xm6OR2[2,3:15], N = 107)
power_tk(alpha1, 1, pmf.Xm6OR2[2,3:15], N = 104)

power_tk(alpha1, K, pmf.Xm6OR2[3,3:15], N = 109)
power_tk(alpha1, 1, pmf.Xm6OR2[3,3:15], N = 106)

power_tk(alpha1, K, pmf.Xm6OR2[4,3:15], N = 102)
power_tk(alpha1, 1, pmf.Xm6OR2[4,3:15], N = 99)

power_tk(alpha1, K, pmf.Xm6OR2[5,3:15], N = 117)
power_tk(alpha1, 1, pmf.Xm6OR2[5,3:15], N = 114)

