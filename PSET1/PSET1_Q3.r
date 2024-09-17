# 1_1 Generate a dataset with n observations ------------------------------------------------------------------

set.seed(123)
# Load necessary libraries
library(mvtnorm)
library(hdm)
library(glmnet)

# Set parameters
set.seed(123)  # For reproducibility
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
results <- matrix(0, nrow = num_simulations, ncol = 6)  # 6 methods
colnames(results) <- c("Lasso", "Post-Lasso", "Square-root Lasso", "Post-Square-root Lasso",
                       "CV Lasso", "CV Post-Lasso")

# Simulation
for (s in 1:num_simulations) {
  # Generate data
  data <- generate_data()
  X <- data$X
  Y <- data$Y
  
  # Lasso estimation
  lasso_fit <- rlasso(X, Y)
  lasso_coefs <- coef(lasso_fit)[-1]  # Exclude intercept
  results[s, "Lasso"] <- mean((X %*% lasso_coefs - X %*% beta0)^2)
  
  # Post-Lasso
  selected <- which(lasso_coefs != 0)
  post_lasso_fit <- lm(Y ~ X[, selected])
  post_lasso_coefs <- rep(0, K)
  post_lasso_coefs[selected] <- coef(post_lasso_fit)[-1]
  results[s, "Post-Lasso"] <- mean((X %*% post_lasso_coefs - X %*% beta0)^2)
  
  # Square-root Lasso
  sqrt_lasso_fit <- rlasso(X, Y, sqrt = TRUE)
  sqrt_lasso_coefs <- coef(sqrt_lasso_fit)[-1]
  results[s, "Square-root Lasso"] <- mean((X %*% sqrt_lasso_coefs - X %*% beta0)^2)
  
  # Post-Square-root Lasso
  selected_sqrt <- which(sqrt_lasso_coefs != 0)
  post_sqrt_lasso_fit <- lm(Y ~ X[, selected_sqrt])
  post_sqrt_lasso_coefs <- rep(0, K)
  post_sqrt_lasso_coefs[selected_sqrt] <- coef(post_sqrt_lasso_fit)[-1]
  results[s, "Post-Square-root Lasso"] <- mean((X %*% post_sqrt_lasso_coefs - X %*% beta0)^2)
  
  # CV Lasso
  cv_lasso_coefs <- cv_lasso(X, Y)
  results[s, "CV Lasso"] <- mean((X %*% cv_lasso_coefs - X %*% beta0)^2)
  
  # CV Post-Lasso
  selected_cv <- which(cv_lasso_coefs != 0)
  post_cv_lasso_fit <- lm(Y ~ X[, selected_cv])
  post_cv_lasso_coefs <- rep(0, K)
  post_cv_lasso_coefs[selected_cv] <- coef(post_cv_lasso_fit)[-1]
  results[s, "CV Post-Lasso"] <- mean((X %*% post_cv_lasso_coefs - X %*% beta0)^2)
}

# Average prediction errors
avg_prediction_errors <- colMeans(results)

# Tabulate results
tabulated_results <- data.frame(
  Estimator = colnames(results),
  Average_Prediction_Error = round(avg_prediction_errors, 4)
)

print(tabulated_results)
