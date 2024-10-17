library(dplyr)
library(mvtnorm)

### Theoretical Type I error and Power Analysis ###
# theta_fun is applied for both Type I and Power analysis 
theta_fun <- function(p) {
  m <- (length(p)-1)/2
  first_sum <- 0
  second_sum <- 0
  p <- as.numeric(p)
  for (l in -m:m) {
    inner_sum <- 0
    for (j in max(-l+1, -m):m) {
      # Check if the index is within the bounds of the vector
      index <- j + m + 1
      if (index >= 1 && index <= length(p)) {
        inner_sum <- inner_sum + p[index]
      }
    }
    l_index <- l + m + 1
    minus_l_index <- -l + m + 1
    if (l_index >= 1 && l_index <= length(p) && minus_l_index >= 1 && minus_l_index <= length(p)) {
      first_sum <- first_sum + inner_sum * p[l_index]
      second_sum <- second_sum + p[minus_l_index] * p[l_index]
    }
  }
  theta <- first_sum + 0.5 * second_sum
  return(theta)
}

#Estimated Theta
varSXiPXi_6 <- function(p) {
  m <- (length(p)-1)/2
  sum1 <- 0
  sum2 <- 0
  sum3 <- 0
  theta_value <- theta_fun(p)
  p <- as.numeric(p)
  for (l in -m:m) {
    inner_sum <- 0
    for (j in max(-l+1, -m):m) {
      index <- j + m + 1
      if (index >= 1 && index <= length(p)) {
        inner_sum <- inner_sum + p[index]
      }
    }
    l_index <- l + m + 1
    minus_l_index <- -l + m + 1
    if (l_index >= 1 && l_index <= length(p) && minus_l_index >= 1 && minus_l_index <= length(p)) {
      sum1 <- sum1 + inner_sum^2 * p[l_index]
      sum2 <- sum2 + 0.25 * p[minus_l_index]^2 * p[l_index]
      sum3 <- sum3 + inner_sum * p[minus_l_index] * p[l_index]
    }
  }
  var_value <- sum1 + sum2 + sum3 - theta_value^2
  return(var_value)
}

#Est theta
E_minusEsquare_7 <- function(p) {
  m <- (length(p) - 1) / 2
  sum1 <- 0
  sum2 <- 0
  theta_value <- theta_fun(p)
  p <- as.numeric(p)
  for (l in -m:m) {
    inner_sum <- 0
    for (j in max(-l+1, -m):m) {
      index <- j + m + 1
      if (index >= 1 && index <= length(p)) {
        inner_sum <- inner_sum + p[index]
      }
    }
    l_index <- l + m + 1
    minus_l_index <- -l + m + 1
    if (l_index >= 1 && l_index <= length(p)) {
      sum1 <- sum1 + inner_sum * p[l_index]
    }
    if (minus_l_index >= 1 && minus_l_index <= length(p)) {
      sum2 <- sum2 + 0.25 * p[minus_l_index] * p[l_index]
    }
  }
  result <- sum1 + sum2 - theta_value^2
  return(result)
}

TheVarUk <- function(K, p, N){
  k <- 1:K
  nk <- rev(N/k) # Sample Size in each Stage
  VarU <- 4* (nk-2)/(nk*(nk-1))*varSXiPXi_6(p) + 2/(nk*(nk-1))*E_minusEsquare_7(p)
  return(VarU)
}
#TheVarUk(3, p , 300)

tk <- function(K, p, N){
  VarUk <- TheVarUk(K, p, N)
  tk <- VarUk[K]/VarUk
  return(tk)
}

alpha_tk <- function(alpha, K, p, N) {
  tk_values <- tk(K ,p , N)
  alpha.tk <- alpha *tk_values^2
  return(alpha.tk)
}
#alpha_tk(0.05, 3, p, 300)

# Return Alpha Spending in Each Interim Analysis k
alpha_k <- function(alpha, K, p, N){
  alpha_raw <- alpha_tk(alpha, K, p, N)
  alpha_k <- c(alpha_raw[1], diff(alpha_raw))
  return(alpha_k)
}

Sigma_mat <- function(K, p, N){
  varVec <- TheVarUk(K, p, N)
  matrix <- diag(1, K, K)
  for (u in 1:(K-1)) {
    for (v in (u+1):K) {
      sqrtRatio <- sqrt(varVec[v] / varVec[u])
      matrix[u, v] <- sqrtRatio
      matrix[v, u] <- sqrtRatio
    }
  }
  return(matrix)
}

critical_values <- function(alpha, K, p, N) {
  if (K == 1) {
    return(qnorm(1 - alpha))
  }
  
  alpha_K <- alpha_k(alpha, K, p, N)
  Sigma <- Sigma_mat(K, p, N)
  if (!all(eigen(Sigma)$values > 0)) {
    stop("Sigma must be positive definite.")
  }
  c_values <- numeric(K)
  c_values[1] <- qnorm(1 - alpha_K[1])
  for (t in 2:K) {
    fc <- function(ct) {
      as.numeric(pmvnorm(upper = c_values[1:t-1], sigma = Sigma[1:t-1, 1:t-1])) - 
        as.numeric(pmvnorm(upper = c(c_values[1:t-1], ct), sigma = Sigma[1:t, 1:t])) - alpha_K[t]
    }
    c_values[t] <- uniroot(fc, c(0, 10))$root
  }
  return(c_values)
}

theo_RejectProb_k <- function(alpha, K, p, N) {
  if (K == 1) {
    c_value <- critical_values(alpha, K, p, N)
    mu <- (theta_fun(p) - 0.5) / sqrt(TheVarUk(K, p, N))
    RP_1 <- 1 - pnorm(c_value, mean = mu, sd = 1)
    return(RP_1)
  }
  c_values <- critical_values(alpha, K, p, N)
  mu <- (theta_fun(p) - 0.5) / sqrt(TheVarUk(K, p, N))
  Sigma <- Sigma_mat(K, p, N)
  RP_1 <- 1 - pmvnorm(upper = c_values[1], mean = mu[1], sigma = Sigma[1, 1])
  RP_k <- numeric(length = K)
  RP_k[1] <- RP_1
  for (k in 2:K) {
    lower_limits <- c(rep(-Inf, k-1), c_values[k]) 
    upper_limits <- c(c_values[1:(k-1)], Inf)  
    RP_k[k] <- pmvnorm(lower = lower_limits, upper = upper_limits, mean = mu[1:k], sigma = Sigma[1:k, 1:k])
  }
  return(RP_k)
}

# Return theoretical power for each interim. Def power[k] = sum RP_K[1:k].
power_tk <- function(alpha, K ,p ,N){
  RP_K <- theo_RejectProb_k(alpha, K ,p ,N)
  power <- numeric(length(RP_K))
  for (i in 1:length(RP_K)) {
    power[i] <- sum(RP_K[1:i])
  }
  return(power)
}

# Simulate Data for a given distribution
phi_XiXj <- function(x_i, x_j) {
  sum_x = x_i + x_j
  if (sum_x > 0) {
    return(1)
  } else if (sum_x == 0) {
    return(1/2)
  } else {
    return(0)}
}

U_statistics <- function(x) {
  n <- length(x); u_sum <- 0
  for (i in 1:(n-1)) {
    for (j in (i+1):n) 
    {u_sum <- u_sum + phi_XiXj(x[i], x[j])
    }}
  u_stat <- u_sum / choose(n, 2)
  return(u_stat)
}

sample_discrete <- function(N, prob) {
  l <- (length(prob) - 1)/2
  values <- -l : l
  simulated_value <- sample(values, size = N, replace = TRUE, prob = prob)
  return(simulated_value)
}

# Return Empirical prob for the given data set
Empirical_prob <- function(simulated_data_partial, l) {
  values <- -l:l
  simulated_data_partial_factor <- factor(simulated_data_partial, levels = values)
  table_data <- table(simulated_data_partial_factor)
  probabilities <- table_data / sum(table_data)
  probabilities_df <- as.data.frame(probabilities)
  names(probabilities_df) <- c("Value", "Probability")
  probabilities_df_formatted <- t(probabilities_df$Probability)
  colnames(probabilities_df_formatted) <- paste("X", probabilities_df$Value, sep="=")
  return(probabilities_df_formatted)
}

# Return a list of Empirical prob distribution for each interim analysis
Empirical_prob_k <- function(simulated_data, K, l) {
  N <- length(simulated_data)
  segment_length <- N %/% K
  values <- -l:l 
  cols <- paste("X", values, sep="=")
  empirical_probs_df <- data.frame(matrix(ncol = length(cols), nrow = K))
  colnames(empirical_probs_df) <- cols
  for (j in 1:K) {
    end_index <- min(j * segment_length, N)
    cumulative_data <- simulated_data[1:end_index]
    empirical_probs <- Empirical_prob(cumulative_data, l)
    current_probs <- setNames(rep(0, length(cols)), cols)
    if (ncol(empirical_probs) > 0) {
      probs_values <- colnames(empirical_probs)
      current_probs[probs_values] <- as.numeric(empirical_probs)
    }
    empirical_probs_df[j, ] <- current_probs
  }
  return(empirical_probs_df)
}

# Simulate data and results test statistics
Test_Stat <- function(N, p, K){
  l <- (length(p)-1)/2
  simulated_data <- sample_discrete(N, p)
  EmProb_k <- Empirical_prob_k(sample_discrete(N, p), K = K, l = l)
  U_stat <- NULL; T_stat <- NULL; Emp_VarUk <- NULL
  segment_length <- N %/% K
  for (k in 1:K) {
    end_index <- min(k * segment_length, N)
    cumulative_data <- simulated_data[1:end_index]
    U_stat[k] <- U_statistics(x = cumulative_data)
    Emp_VarUk[k] <- TheVarUk(K = K, p = EmProb_k[k,], N = N)[k] 
    T_stat[k] <- (U_stat[k]- 1/2) / sqrt(Emp_VarUk[k])
  }
  return(T_stat)
}

Attained_level <- function(B, N, p, K, alpha) {
  ctc_value <- critical_values( alpha, K, p, N)
  exceed_counts <- numeric(K)
  valid_B <- 0  
  pb <- txtProgressBar(min = 0, max = B, style = 3)
  for (b in 1:B) {
    test_stat <- Test_Stat(N, p, K)
    if(any(is.infinite(test_stat) | is.nan(test_stat))) {
      setTxtProgressBar(pb, b)
      next
    }
    count_B <- numeric(K)
    for (k in 1:K){
      count_B[k] <- sum((test_stat > ctc_value)[1:k]) > 0
    }
    exceed_counts <- exceed_counts + count_B
    valid_B <- valid_B + 1
    setTxtProgressBar(pb, b)
  }
  close(pb)
  if (valid_B > 0) {
    alpha_attained <- exceed_counts / valid_B 
  } else {
    alpha_attained <- rep(NA, K)
  }
  return(alpha_attained)
}
