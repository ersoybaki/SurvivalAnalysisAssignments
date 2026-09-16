set.seed(3131)

# Parameters
true_rate <- 0.05 # True lambda
n <- 500 
R <- 2000


T_cal <- rexp(n, rate = true_rate)



c_max_roots <- numeric(3)
pc_targets <- c(0, 0.1, 0.5)


count <- 1
for (pc in pc_targets) {
  # No censoring if pc = 0
  if (pc == 0) {
    c_max_roots[count] <- Inf
    count <- count + 1
    next
  }
  
  # Function to find c_max
  # ! This function was generated with Gemini 3.6 Flash (Ersoy). I will add it to the report ! 
  obj_fn <- function(c_max) {
    C_sim <- runif(n, 0, c_max)
    empirical_pc <- mean(C_sim < T_cal)
    return(empirical_pc - pc)
  }
  
  # ! End of AI Generated code !
  
  C_sim_fit <- uniroot(obj_fn, interval = c(1, 2000))
  c_max_roots[count] <- C_sim_fit$root
  count <- count + 1
}

# Sets
mean_Y <- numeric(R)
mle_est_x   <- numeric(R)
cens_prop_x <- numeric(R)

# For p_c = 0, 0.1, 0.5 creating variables

# To store estimations and censoring proportion
naive_est <- numeric(3)
mle_est <- numeric(3)
cens_prop <- numeric(3)
# Bias set
naive_bias <- numeric(3)
mle_bias   <- numeric(3)
#Variance set
naive_var <- numeric(3)
mle_var   <- numeric(3)

# To store MSE calculations
naive_est_mse <- numeric(3) 
mle_est_mse   <- numeric(3)


count <- 1

for (c_max in c_max_roots) {
  for (i in 1:R) {
    # Times 
    T_i <- rexp(n, rate=true_rate)
    
    # Right Censoring
    if (c_max == Inf) {
      C <- rep(Inf, n)
    } else {
      C <- runif(n, 0, c_max)
    }

    # Observed data
    Y <- pmin(T_i, C)
    
    # Mean Y calculation for Naive est
    mean_Y[i] <- sum(Y) / n
    
    # Delta calculation
    delta <- as.integer(T_i <= C)
    
    # Lamda_mle
    r <- sum(delta)
    if (r == 0 ) {
      mle_est_x[i] <- NA
    } else {
      mle_est_x[i] <- r / sum(Y)
    }
    
    
    # average censoring proportion (observed across replicates)
    cens_prop_x[i] <- 1 - mean(delta)
  }
  
  # Estimators
  naive_est[count] <- mean(1 / mean_Y)
  mle_est[count] <- mean(mle_est_x, na.rm = TRUE)
  
  #Bias
  naive_bias[count] <- naive_est[count] - true_rate
  mle_bias[count]   <- mle_est[count]   - true_rate
  
  # Variance
  naive_var[count] <- var(1 / mean_Y)
  mle_var[count]   <- var(mle_est_x, na.rm = TRUE)
  
  # MSE Calculations
  naive_est_mse[count] <- mean(((1 / mean_Y) - true_rate)^2)
  mle_est_mse[count] <- mean((mle_est_x - true_rate)^2, na.rm = TRUE)
  
  # Censoring proportion
  cens_prop[count] <- mean(cens_prop_x)
  
  count <- count + 1
  
}

# Print the findings as a Data Frame
data.frame(pc_target = c(0, 0.1, 0.5), c_max = c_max_roots,
           cens_prop, naive_est, mle_est, naive_est_mse, mle_est_mse, naive_bias, mle_bias, naive_var, mle_var)


