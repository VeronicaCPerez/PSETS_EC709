# 1_1 Generate a dataset with n observations ------------------------------------------------------------------

set.seed(123)

# Load necessary libraries
library(mvtnorm)
library(hdm)
library(glmnet)

# Set parameters
n <- 100
K <- 500
sigma2 <- 1
beta0 <- c(1, 1, 1/2, 1/3, 1/4, 1/5, rep(0, K - 6))
Sigma <- outer(1:K, 1:K, function(i, j) (1/2)^abs(i - j))
num_simulations <- 500
num_folds <- 5  # For cross-validation

# Function to generate data
generate_data <- function() {
  X <- rmvnorm(n, mean = rep(0, K), sigma = Sigma)
  epsilon <- rnorm(n, mean = 0, sd = sqrt(sigma2))
  Y <- X %*% beta0 + epsilon
  list(X = X, Y = Y)
}

# Function to perform cross-validated Lasso
cv_lasso <- function(X, Y) {
  # Perform 5-fold cross-validation
  fit <- cv.glmnet(X, Y, alpha = 1, nfolds = num_folds)
  return(coef(fit, s = "lambda.min")[-1])  # Exclude intercept
}

# Storage for results
results <- matrix(0, nrow = num_simulations, ncol = 4)  # 6 methods
bias <- matrix(0, nrow = num_simulations, ncol = 4)  # To store bias
coef_results <- matrix(0, nrow = num_simulations, ncol = 4)  # Store coefficients

colnames(coef_results) <- c("Lasso", "Post-Lasso","CV Lasso", "CV Post-Lasso")
colnames(results) <- c("Lasso", "Post-Lasso","CV Lasso", "CV Post-Lasso")
colnames(bias) <- c("Lasso", "Post-Lasso","CV Lasso", "CV Post-Lasso")

# Storage for coefficients


# Simulation
for (s in 1:num_simulations) {
  # Generate data
  data <- generate_data()
  X <- data$X
  Y <- data$Y
  
  # Lasso estimation
  lasso_fit <- rlasso(X, Y, post = FALSE)
  lasso_coefs <- coef(lasso_fit)[-1]  # Exclude intercept
  results[s, "Lasso"] <- mean((X %*% lasso_coefs - X %*% beta0)^2)
  bias[s, "Lasso"] <- mean(lasso_coefs - beta0)  # Compute bias
  coef_results[s, "Lasso"] <- mean(lasso_coefs)
  
  # Post-Lasso
  post_lasso_fit <- rlasso(X, Y, post = TRUE)
  post_lasso_coefs <- coef(post_lasso_fit)[-1]  # Exclude intercept
  results[s, "Post-Lasso"] <- mean((X %*% post_lasso_coefs - X %*% beta0)^2)
  bias[s, "Post-Lasso"] <- mean(post_lasso_coefs - beta0)  # Compute bias
  coef_results[s, "Post-Lasso"] <- mean(post_lasso_coefs)
  
  # CV Lasso
  cv_lasso_coefs <- cv_lasso(X, Y)
  results[s, "CV Lasso"] <- mean((X %*% cv_lasso_coefs - X %*% beta0)^2)
  bias[s, "CV Lasso"] <- mean(cv_lasso_coefs - beta0)  # Compute bias
  coef_results[s, "CV Lasso"] <- mean(cv_lasso_coefs)
  
  # CV Post-Lasso
  selected_cv <- which(cv_lasso_coefs != 0)
  post_cv_lasso_fit <- lm(Y ~ X[, selected_cv])
  post_cv_lasso_coefs <- rep(0, K)
  post_cv_lasso_coefs[selected_cv] <- coef(post_cv_lasso_fit)[-1]
  results[s, "CV Post-Lasso"] <- mean((X %*% post_cv_lasso_coefs - X %*% beta0)^2)
  bias[s, "CV Post-Lasso"] <- mean(post_cv_lasso_coefs - beta0)  # Compute bias
  coef_results[s, "CV Post-Lasso"] <- mean(post_cv_lasso_coefs)
  
  print(s)
}

# Average prediction errors and coefficients
avg_prediction_errors <- colMeans(results)
avg_bias <- colMeans(bias)
avg_coefs <- colMeans(coef_results)

# Tabulate results
tabulated_results <- data.frame(
  Estimator = round(avg_coefs,4),
  Average_Prediction_Error = round(avg_prediction_errors, 4),
  Bias = round(avg_bias, 4)
)

# Print the results
print(tabulated_results)

# Convert to LaTeX format
library(xtable)

# Create the LaTeX table with the caption
latex_table <- xtable(tabulated_results, caption = "Lasso, prediction error, bias, and average coefficients", digits=c(0,4,4,4))

# Save LaTeX code to a file with 4 decimal precision
sink("/Users/veronica/Dropbox/Apps/Overleaf/EC_709_vcperez/PSET_1/Q3_table.tex")
print(latex_table, include.rownames = FALSE, booktabs = TRUE)
sink()
