library(ggplot2)
library(tidyr)
library(xtable)

######## Distribution of X^(1)
calc_probs <- function(p0, OR) {
  l <- length(p0)
  cum_p0 <- cumsum(p0)[1:l-1]
  cum_p1 <- 1/(OR*(1-cum_p0)/cum_p0 + 1)
  cum_p1 <- c(cum_p1, 1)
  p1 <- c(cum_p1[1], diff(cum_p1))
  return(p1)
}
###### Distribution of X
calc_diff_probs <- function(p0, p1) {
  m <- length(p0) -1 
  x_values <- -m : m
  probs <- numeric(length(x_values))
  for (i in seq_along(x_values)) {
    x_val <- x_values[i]
    for (j in 0 : m) {
      k <- j + x_val
      if (k >= 0 && k <= m) {
        probs[i] <- probs[i] + p0[j + 1] * p1[k + 1]
      }
    }
  }
  return(probs)
}

# Function to generate results_X1 and process data
process_results_X1 <- function(baseline_probs, OR_values) {
  results_X1 <- data.frame(stringsAsFactors = FALSE)
  
  for (scenario in names(baseline_probs)) {
    for (OR in OR_values) {
      p0 <- baseline_probs[[scenario]]
      p1 <- calc_probs(p0, OR)
      new_row <- c(Scenario = scenario, OR = OR, as.numeric(format(p1, digits = 3)))
      results_X1 <- rbind(results_X1, new_row)
    }
  }
  
  # Dynamically create column names
  colnames(results_X1) <- c("Scenario", "OR", paste0("p_{", seq_along(baseline_probs[[1]]) - 1, "}"))
  results_X1[-(1:2)] <- apply(results_X1[-(1:2)], 2, as.numeric)
  
  return(results_X1)
}

# Function to generate results_X and process data
process_results_X <- function(baseline_probs, OR_values) {
  results_X <- data.frame(stringsAsFactors = FALSE)
  
  for (scenario in names(baseline_probs)) {
    p0 <- baseline_probs[[scenario]]
    for (OR in OR_values) {
      p1 <- calc_probs(p0, OR)
      diff_probs <- calc_diff_probs(p0, p1)
      
      # Create a data frame for the new row with correct data types
      new_row <- data.frame(
        Scenario = scenario,
        OR = as.numeric(OR),
        t(diff_probs),
        stringsAsFactors = FALSE
      )
      
      results_X <- rbind(results_X, new_row)
    }
  }
  num_prob_columns <- length(diff_probs)
  prob_column_names <- paste0("p_{", seq(-(num_prob_columns - 1) / 2, (num_prob_columns - 1) / 2), "}")
  colnames(results_X) <- c("Scenario", "OR", prob_column_names)
  
  # Convert probability columns to numeric
  prob_columns <- prob_column_names
  results_X[prob_columns] <- lapply(results_X[prob_columns], as.numeric)
  
  return(results_X)
}

# Function to output LaTeX table from results
output_latex_table <- function(results) {
  latex_table <- xtable(results, digits = 3)
  print(latex_table, include.rownames = FALSE, hline.after = c(-1, 0, nrow(results)), 
        booktabs = TRUE, comment = FALSE)
}

# Plot PMF Function
plot_pmf <- function(results, x_label = "Level of X", y_label = "Probability", color_label = "Odds Ratio") {
  # Convert data to long format for plotting
  probability_columns <- grep("^p_", colnames(results), value = TRUE)
  long_results <- tidyr::pivot_longer(results, cols = all_of(probability_columns),
                                      names_to = "X_level", values_to = "Probability")
  
  # Extract numeric X_level from column names like "p_{-2}"
  long_results$X_level <- as.numeric(sub("p_\\{(-?\\d+)\\}", "\\1", long_results$X_level))
  
  # Convert Probability to numeric and handle missing values
  long_results$Probability <- suppressWarnings(as.numeric(long_results$Probability))
  long_results <- na.omit(long_results)
  
  # Convert OR and Scenario to factors
  long_results$OR <- as.factor(long_results$OR)
  if("Scenario" %in% names(long_results)){
    long_results$Scenario <- as.factor(long_results$Scenario)
  }
  
  # Plot PMF
  p <- ggplot(long_results, aes(x = X_level, y = Probability, group = OR, color = OR)) +
    geom_line() +
    theme_classic() +
    labs(x = x_label,
         y = y_label,
         color = color_label) +
    theme(legend.position = "bottom") +
    scale_color_discrete(name = color_label) +
    scale_y_continuous(breaks = function(x) pretty(x, n = 5))
  
  # Add faceting if Scenario is available
  if("Scenario" %in% names(long_results)){
    p <- p + facet_wrap(~Scenario, scales = "free")
  }
  
  print(p)
}
