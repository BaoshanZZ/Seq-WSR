# Simulation Generation and Results Plots #
sim_by_OR <- function(OddsR, B, results_X, K, alpha, N_values) {
  X_probOR <- subset(results_X, OR == OddsR)
  results_list <- list()
  result_id <- 1
  m <- (ncol(results_X) - 3)/2
  for (N in N_values) { 
    for (i in 1:nrow(X_probOR)) {
      p <- as.numeric(X_probOR[i, 3:ncol(results_X)])
      # Both Type I and power could be applied here. OR = 1 -> Type I. Other II
      if (OddsR == 1) {
        theoretical_level.seq <- alpha_tk(alpha, K, p, N)
        theoretical_level.fix <- alpha_tk(alpha, 1, p, N)
      } else {
        theoretical_level.seq <- power_tk(alpha, K, p, N)
        theoretical_level.fix <- power_tk(alpha, 1, p, N)}
      attained_level.seq <- Attained_level(B, N, p, K, alpha)
      bias.seq <- attained_level.seq - theoretical_level.seq
      theoretical_level.fix <- power_tk(alpha, 1, p, N)
      attained_level.fix <- Attained_level(B, N, p, 1, alpha)
      bias.fix <- attained_level.fix - theoretical_level.fix
      scenario_result <- c(Scenario = X_probOR$Scenario[i],
                           m,
                           OR = OddsR, 
                           N = N, 
                           K = K,
                           attained.Seq = attained_level.seq,
                           theoretical.Seq = theoretical_level.seq,
                           bias.Seq = bias.seq,
                           attained.Fix = attained_level.fix,
                           theoretical.Fix = theoretical_level.fix,
                           bias.Fix = bias.fix)
      names(scenario_result) <- c("Scenario", "Scale_Level", "OR", "N", "K",
                                  paste0("Attained_", 1:K),
                                  paste0("Theoretical_", 1:K),
                                  paste0("Bias_", 1:K),
                                  "Attained_Fix", "Theoretical_Fix", "Bias_Fix")
      results_list[[result_id]] <- scenario_result
      cat("OR:", OddsR, "N:", N, "Scenario:", X_probOR$Scenario[i], "\n")
      result_id <- result_id + 1
    }
  }
  results_df <- do.call(rbind, 
                        lapply(results_list, 
                               function(x) as.data.frame(t(x), stringsAsFactors = FALSE)))
  return(results_df)
}


#### Simulation Results Plot Function
SimPlot_Bias <- function(results_df) {
  # Extract unique N values and sort them (assuming numeric N)
  N_levels <- sort(unique(as.numeric(as.character(results_df$N))))
  N_levels <- as.character(N_levels)
  
  # Set N as a factor with specified levels
  results_df$N <- factor(results_df$N, levels = N_levels)
  
  # Specify the bias columns to include in the plot (exclude Bias_Fix)
  bias_columns <- grep("^Bias_[0-9]+$", names(results_df), value = TRUE)
  
  # Reshape the data from wide to long format
  long_df <- melt(results_df, 
                  id.vars = c("Scenario", "OR", "N", "K"),
                  measure.vars = bias_columns,
                  variable.name = "Stage", 
                  value.name = "Bias")
  
  # Extract the stage number from the variable names
  long_df$Stage <- as.numeric(gsub("Bias_", "", long_df$Stage))
  
  # Convert Bias to numeric for plotting
  long_df$Bias <- as.numeric(long_df$Bias)
  
  # Create the plot
  p <- ggplot(long_df, aes(x = Stage, y = Bias, group = Scenario, color = Scenario)) +
    geom_line() +
    geom_point() +
    facet_wrap(~N, scales = "free_x", labeller = labeller(N = function(x) paste("N =", x))) +
    labs(x = "Stage", y = "Bias (Attained - Theoretical)") +
    theme_classic() +
    theme(legend.position = "bottom",
          legend.title = element_blank(),
          strip.text.x = element_text(face = "bold"),
          plot.title = element_text(hjust = 0.5),
          axis.title = element_text(size = 12),
          axis.text = element_text(size = 10),
          legend.text = element_text(size = 10)) +
    scale_x_continuous(breaks = 1:max(long_df$Stage), labels = function(x) as.character(x))
  # +  scale_y_continuous(breaks = function(x) pretty(x, n = 5), labels = function(x) format(x, scientific = FALSE))
  
  return(p)
}


SimPlot_AttainedPower <- function(results_df) {
  # Extract unique N values and sort them (assuming numeric N)
  N_levels <- sort(unique(as.numeric(as.character(results_df$N))))
  
  # Convert N back to character if needed
  N_levels <- as.character(N_levels)
  
  # Set N as a factor with levels based on unique N values
  results_df$N <- factor(results_df$N, levels = N_levels)
  
  # Specify the columns to include in the plot
  # Exclude any columns like Attained_Fix if necessary
  attained_columns <- grep("^Attained_[0-9]+$", names(results_df), value = TRUE)
  
  # Reshape the data from wide to long format
  long_df <- melt(results_df, 
                  id.vars = c("Scenario", "OR", "N", "K"),
                  measure.vars = attained_columns,
                  variable.name = "Stage", 
                  value.name = "AttPower")
  
  # Extract the stage number from the variable names
  long_df$Stage <- as.numeric(gsub("Attained_", "", long_df$Stage))
  
  # Convert AttPower to numeric for plotting
  long_df$AttPower <- as.numeric(long_df$AttPower)
  
  # Create the plot
  p <- ggplot(long_df, aes(x = Stage, y = AttPower, group = Scenario, color = Scenario)) +
    geom_line() +
    geom_point() +
    facet_wrap(~N, scales = "free_x", labeller = labeller(N = function(x) paste("N =", x))) +
    labs(x = "Stage", y = "Attained Power Level") +
    theme_classic() +
    theme(legend.position = "bottom",
          legend.title = element_blank(),
          strip.text.x = element_text(face = "bold"),
          plot.title = element_text(hjust = 0.5),
          axis.title = element_text(size = 12),
          axis.text = element_text(size = 10),
          legend.text = element_text(size = 10)) +
    scale_x_continuous(breaks = sort(unique(long_df$Stage)), labels = function(x) as.character(x)) 
  
  return(p)
}


SimPlot_AttainedAlpha <- function(results_df) {
  # Extract unique N values and sort them (assuming numeric N)
  N_levels <- sort(unique(as.numeric(as.character(results_df$N))))
  
  # Convert N back to character if needed
  N_levels <- as.character(N_levels)
  
  # Set N as a factor with levels based on unique N values
  results_df$N <- factor(results_df$N, levels = N_levels)
  
  # Specify the columns to include in the plot
  stage_columns <- c("Attained_1", "Attained_2", "Attained_3")
  
  # Reshape the data from wide to long format, excluding Attained_Fix
  long_df <- melt(results_df, 
                  id.vars = c("Scenario", "OR", "N", "K"),
                  measure.vars = stage_columns,
                  variable.name = "Stage", 
                  value.name = "AttAlpha")
  
  # Extract the stage number from the variable names
  long_df$Stage <- as.numeric(gsub("Attained_", "", long_df$Stage))
  
  # Convert AttAlpha to numeric for plotting
  long_df$AttAlpha <- as.numeric(long_df$AttAlpha)
  
  # Create the plot
  p <- ggplot(long_df, aes(x = Stage, y = AttAlpha, group = Scenario, color = Scenario)) +
    geom_line() +
    geom_point() +
    facet_wrap(~N, scales = "free_x", labeller = labeller(N = function(x) paste("N =", x))) +
    labs(x = "Stage", y = "Attained Alpha Level") +
    theme_classic() +
    theme(legend.position = "bottom",
          legend.title = element_blank(),
          strip.text.x = element_text(face = "bold"),
          plot.title = element_text(hjust = 0.5),
          axis.title = element_text(size = 12),
          axis.text = element_text(size = 10),
          legend.text = element_text(size = 10)) +
    scale_x_continuous(breaks = 1:max(long_df$Stage), labels = function(x) as.character(x))
  
  return(p)
}

SimPlot_BiasRatio <- function(results_df) {
  # Load necessary libraries
  library(reshape2)
  library(ggplot2)
  
  # Extract unique N values and sort them (assuming numeric N)
  N_levels <- sort(unique(as.numeric(as.character(results_df$N))))
  N_levels <- as.character(N_levels)
  
  # Set N as a factor with specified levels
  results_df$N <- factor(results_df$N, levels = N_levels)
  
  # Identify Bias and Theoretical Power columns
  bias_columns <- grep("^Bias_[0-9]+$", names(results_df), value = TRUE)
  theoretical_columns <- grep("^Theoretical_[0-9]+$", names(results_df), value = TRUE)
  
  # Ensure that the necessary columns are numeric
  results_df[bias_columns] <- lapply(results_df[bias_columns], as.numeric)
  results_df[theoretical_columns] <- lapply(results_df[theoretical_columns], as.numeric)
  
  # Calculate BiasRatio for each stage
  for (bias_col in bias_columns) {
    # Extract stage number
    stage_num <- gsub("Bias_", "", bias_col)
    theoretical_col <- paste0("Theoretical_", stage_num)
    bias_ratio_col <- paste0("BiasRatio_", stage_num)
    
    # Check if the corresponding theoretical column exists
    if (theoretical_col %in% names(results_df)) {
      # Calculate BiasRatio
      results_df[[bias_ratio_col]] <- results_df[[bias_col]] / results_df[[theoretical_col]]
    } else {
      warning(paste("Theoretical column", theoretical_col, "not found."))
    }
  }
  
  # Identify all BiasRatio columns
  bias_ratio_columns <- grep("^BiasRatio_[0-9]+$", names(results_df), value = TRUE)
  
  # Reshape the data from wide to long format
  long_df <- melt(results_df, 
                  id.vars = c("Scenario", "OR", "N", "K"),
                  measure.vars = bias_ratio_columns,
                  variable.name = "Stage", 
                  value.name = "BiasRatio")
  
  # Extract the stage number from the variable names
  long_df$Stage <- as.numeric(gsub("BiasRatio_", "", long_df$Stage))
  
  # Convert BiasRatio to numeric for plotting
  long_df$BiasRatio <- as.numeric(long_df$BiasRatio)
  
  # Create the plot
  p <- ggplot(long_df, aes(x = Stage, y = BiasRatio, group = Scenario, color = Scenario)) +
    geom_line() +
    geom_point() +
    facet_wrap(~N, labeller = labeller(N = function(x) paste("N =", x))) +
    labs(x = "Stage", y = "Bias Ratio (Bias / Theoretical Power)") +
    theme_classic() +
    theme(legend.position = "bottom",
          legend.title = element_blank(),
          strip.text.x = element_text(face = "bold"),
          plot.title = element_text(hjust = 0.5),
          axis.title = element_text(size = 12),
          axis.text = element_text(size = 10),
          legend.text = element_text(size = 10)) +
    scale_x_continuous(breaks = sort(unique(long_df$Stage)), labels = function(x) as.character(x))
  
  return(p)
}
