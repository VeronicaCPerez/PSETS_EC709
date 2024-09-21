library(MASS)  # For mvrnorm function
library(dplyr)
library(xtable)

# Define the Epanechnikov kernel as a vectorized function
epanechnikov_kernel <- function(u) {
  result <- numeric(length(u))
  result[abs(u) <= 1] <- 3/4 * (1 - u[abs(u) <= 1]^2)
  return(result)
}

# Define bandwidth values
bandwidths <- c(0.05, 0.10, 0.15, 0.20)
n <- 1000
x <- 0.5

# Define the number of simulations
num_simulations <- 1000

# Initialize storage for results
kernel_bias <- numeric(length(bandwidths))
kernel_variance <- numeric(length(bandwidths))
llr_bias <- numeric(length(bandwidths))
llr_variance <- numeric(length(bandwidths))

# Simulation loop
for (j in 1:length(bandwidths)) {
  h <- bandwidths[j]
  
  kernel_estimates <- numeric(num_simulations)
  llr_estimates <- numeric(num_simulations)
  
  for (sim in 1:num_simulations) {
    # Generate data
    X <- runif(n)
    epsilon <- rnorm(n)
    Y <- exp(X) * (1 + epsilon)
    
    # Kernel regression estimator
    kernel_estimates_sim <- numeric(n)
    weights <- epanechnikov_kernel((x - X) / h) / sum(epanechnikov_kernel((x - X) / h))
    kernel_estimates_sim <- sum(weights * Y)
    kernel_estimates[sim] <- kernel_estimates_sim
    
    # Locally linear regression estimator
    llr_model <- lm(Y ~ 1 + I(x - X), weights = epanechnikov_kernel((x - X) / h))
    llr_estimates[sim] <- coef(llr_model)[1]
  }
  
  # Compute bias and variance
  true_value <- exp(x)  # since E[Y | X=x] = exp(x)
  kernel_bias[j] <- mean(kernel_estimates) - true_value
  kernel_variance[j] <- var(kernel_estimates)
  llr_bias[j] <- mean(llr_estimates) - true_value
  llr_variance[j] <- var(llr_estimates)
}

# Print results
results <- data.frame(
  Bandwidth = bandwidths,
  Kernel_Bias = kernel_bias,
  Kernel_Variance = kernel_variance,
  LLR_Bias = llr_bias,
  LLR_Variance = llr_variance
)

# Convert to LaTeX format
latex_table <- xtable(results, caption = "Bias and Variance of Kernel and Locally Linear Regression Estimators")

print(latex_table)

# Save LaTeX code to a file
sink("/Users/veronica/Dropbox/Apps/Overleaf/EC_709_vcperez/PSET_1/Q1_2_table.tex")
print(latex_table, include.rownames = FALSE, booktabs = TRUE)
sink()