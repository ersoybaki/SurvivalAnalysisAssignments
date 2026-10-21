# install.packages("faux")
library(faux)
library(survival)

n <- 500
B <- 1000
C <- 12

# Distribution Parameters
age_mean <- 70
age_sd <- 10

ex_mean <- 3.5
ex_sd <- 3

hist_prob <- 0.3


# Beta Parameters
# Determined with LLMs. Prompted to make research and back parameters with real findings
# Possible future work
shape <- 0.8 # Weibull k=0.8 (Section 2 from report states k<1 but no specific value)
beta0 <- log(C / log(2)^(1 / shape)) - beta_age * age_mean - beta_ex * ex_mean # Deign choice: 50% readmitted by 12 months 
beta_age <- -log(1.98) / 10  / shape # HR 1.98 per 10 years (https://pmc.ncbi.nlm.nih.gov/articles/PMC2465638/table/tbl2/)
beta_ex   <-  -log(0.86) / 2.5 / shape    # RR 0.86 for 2.5 h/week (https://pubmed.ncbi.nlm.nih.gov/22039197/)
beta_hist <- -log(1.74) / shape    # HR 1.74 (https://pmc.ncbi.nlm.nih.gov/articles/PMC3773915/)

# Estimations and Coverages
est_true <- est_false <- est_naive_true <- est_naive_false <- numeric(B)
cov_true <- cov_false <- cov_naive_true <- cov_naive_false <- numeric(B)

# Rho (correlation between age and exercise)
# !!! According to RQ2, change rho and rerun the experiments !!!
rho=0.5

# Function to check 95% CI coverage
covers <- function(fit) {
  e  <- coef(fit)["exercise"]
  se <- sqrt(vcov(fit)["exercise", "exercise"])
  abs(e - beta_ex) <= 1.96 * se
}

# Step 1. Main loop to run experiments
for (b in 1:B) {
  
  # Step 2.a. Generate the distributed variables
  covs <- rnorm_multi(n=n, mu=c(age_mean, ex_mean), sd=c(age_sd, ex_sd), r=rho, varnames = c("age", "exercise"))
  covs$history <- rbinom(n, 1, hist_prob)
  
  # Scale via log linear link
  scale <- exp(beta0 + beta_age * covs$age + beta_ex * covs$exercise + beta_hist * covs$history)
  
  # Draw each true event time
  T <- rweibull(n, shape=shape, scale=scale)
  
  # Step 2.b. Application Censoring
  Y <- pmin(T, C)
  delta <- as.integer(T <= C) 
  
  # Step 2.c. Fitting Weibull Regressions (true vs false)
  d <- data.frame(Y, delta, covs)
  
  fit_true <- survreg(Surv(Y, delta) ~ age + exercise + history, data=d, dist="weibull")
  fit_false <- survreg(Surv(Y, delta) ~ exercise + history, data=d, dist="weibull")
  
  # Step 2.d. Save estimates and coverages from both models
  est_true[b] <- coef(fit_true)["exercise"]
  est_false[b] <- coef(fit_false)["exercise"]
  
  cov_true[b] <- covers(fit_true)
  cov_false[b] <- covers(fit_false)
  
  # ALTERNATIVE HYPOTHESIS: Naive approach (no censoring)
  d_naive <- d[d$delta==1,]
  
  fit_naive_true <- survreg(Surv(Y) ~ age + exercise + history, data=d_naive, dist="weibull")
  fit_naive_false <- survreg(Surv(Y) ~ exercise + history, data=d_naive, dist="weibull")
  
  est_naive_true[b] <- coef(fit_naive_true)["exercise"]
  est_naive_false[b] <- coef(fit_naive_false)["exercise"]
  
  cov_naive_true[b]  <- covers(fit_naive_true)
  cov_naive_false[b] <- covers(fit_naive_false)
  
}

# Step 3. Compute bias, variance and MSE
bias <- c(mean(est_true), mean(est_false)) - beta_ex
var <- c(var(est_true), var(est_false))
mse <- bias^2 + var


results <- data.frame(
  Model    = c("Correct (with Age)", "Misspecified (no Age)"),
  Mean_est = c(mean(est_true), mean(est_false)),
  Bias     = bias,
  Variance = var,
  MSE      = mse, 
  Coverage = c(mean(cov_true), mean(cov_false)) 
)

print(results, digits = 4, row.names = FALSE)

# ALTERNATIVE HYPOTHESIS: Computer bias, variance, and MSE
bias_naive <- c(mean(est_naive_true), mean(est_naive_false)) - beta_ex
var_naive <- c(var(est_naive_true), var(est_naive_false))
mse_naive <- bias_naive^2 + var_naive


results_naive <- data.frame(
  Model    = c("Correct (with Age)", "Misspecified (no Age)"),
  Mean_est = c(mean(est_naive_true), mean(est_naive_false)),
  Bias     = bias_naive,
  Variance = var_naive,
  MSE      = mse_naive, 
  Coverage = c(mean(cov_naive_true), mean(cov_naive_false)) 
)

cat("--- ALTERNATIVE HYPOTHESIS --- \n")
print(results_naive, digits = 4, row.names = FALSE)
