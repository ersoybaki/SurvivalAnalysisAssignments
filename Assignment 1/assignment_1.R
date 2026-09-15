library(survival)
set.seed(3131)

# Parameters
true_rate <- 0.05 # True lambda
n <- 2000 
target_pc <- 0.5 # CHANGE TO 0, 0.1, 0.5
R <- 2000

# Estimators
# naive_est <- numeric(1)

T_cal <- rexp(n, rate = true_rate)

# Function to find c_max
# ! This function was generated with Gemini 3.6 Flash (Ersoy). I will add it to the report ! 
obj_fn <- function(c_max) {
  C_sim <- runif(n, 0, c_max)
  empirical_pc <- mean(C_sim < T_cal)
  return(empirical_pc - target_pc)
}
C_sim_fit <- uniroot(obj_fn, interval = c(1, 2000))
# ! End of AI Generated code !

# Sets
mean_Y <- numeric(R)
mle_est   <- numeric(R)
cens_prop <- numeric(R)


for (i in 1:R) {
  # Times 
  T_i <- rexp(n, rate=true_rate)
  
  # Right Censoring
  C <- runif(n, 0, C_sim_fit$root)
  
  # Observed data
  Y <- pmin(T_i, C)
  
  # Mean Y calculation for Naive est
  mean_Y[i] <- sum(Y) / n
  
  # Delta calculation
  delta <- as.integer(T_i <= C)
  
  # Lamda_mle
  r <- sum(delta)
  mle_est[i] <- r / sum(Y)
  
  # average censoring proportion (observed across replicates)
  cens_prop[i] <- 1 - mean(delta)
}

naive_est <- 1 / mean_Y

