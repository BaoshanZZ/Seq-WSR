## Brief Introduction of how to use this Seq-WSR document.

### FUNC_PMF.R:
- Generate PMF of $X^{(1)}$ and $X$ for different scenarios of baseline distribution setting of $X^{(0)}$;
- Plot the PMF figure and also generate the latex code for the PMF table.

### FUNC_SeqWSR.R
- Using for the estimated $\theta$ and its corresponding variance estimation
- Alpha spending function and critical value (or boundary)
- Variance-covariance matrix of seq-WSR
- Theoretical Power and Rejection Probability Calculation

### FUNC_Sim.R
- Simulation Generation and Results Plots
- Return Empirical prob for the given data set
- Simulate data and results test statistics
- Return the Attained Level of a simulation study in each interim analysis

### Calc_Simulation.R
- Use to Do the simulation study in the manuscript for different settings of $X^{(0)}$
- Generate the simulation results and plots for attained Power/Type I error level
  
### Calc_FixSeq.R
- Use to compare the fixed design vs sequential design for a given power
- generate the theoretical power level in each interim analysis
