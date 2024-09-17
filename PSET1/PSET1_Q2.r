##############################################################
##############################################################
#                       Question 2_1
##############################################################
##############################################################


# KERNEL ------------------------------------------------------------------

# Load necessary libraries
set.seed(123)
library(foreign)
library(haven)
library(KernSmooth)
library(gplm)
library(locpol)

# Load and process the data
data <- read.dta("/Users/veronica/Dropbox/PhD/2024_2/EC_709/PSETS_EC709/PSET1/yn.dta")
data$log_gas <- log(data$gas)
data$log_price <- log(data$price)

# Extract the variables
X <- data$log_price
Y <- data$log_gas

# Kernel Regression
bandwidths <- c(0.01, 0.05, 0.08, 0.1, 0.15, 0.2, 0.25, 0.3)

# Perform kernel regression and plot the results
png("/Users/veronica/Dropbox/Apps/Overleaf/EC_709_vcperez/PSET_1/Q2_1_kernel_regression_comparison.png")
plot(NULL, xlim = range(X), ylim = range(Y), 
     main = "Kernel Regression with Different Bandwidths",
     xlab = "Log Price", ylab = "Log Gasoline Consumption")

colors <- c("orange", "purple", "limegreen", "blue", "red", "yellow", "green", "cyan")
for (i in seq_along(bandwidths)) {
  Kreg <- ksmooth(x = X, y = Y, kernel = "normal", bandwidth = bandwidths[i])
  lines(Kreg, lwd = 2, col = colors[i])
}
legend("topright", legend = paste("h =", bandwidths), lwd = 2, col = colors)
dev.off()

# Cross-validation for bandwidth selection
n <- length(X)
h_seq <- seq(from = 0.01, to = 0.3, by = 0.01)
CV_err_h <- numeric(length(h_seq))

for (j in seq_along(h_seq)) {
  h_using <- h_seq[j]
  CV_err <- numeric(n)
  
  for (i in 1:n) {
    X_val <- X[i]
    Y_val <- Y[i]
    X_tr <- X[-i]
    Y_tr <- Y[-i]
    
    Y_val_predict <- ksmooth(x = X_tr, y = Y_tr, kernel = "normal", bandwidth = h_using, x.points = X_val)$y
    CV_err[i] <- (Y_val - Y_val_predict)^2
  }
  
  CV_err_h[j] <- mean(CV_err)
}

# Optimal bandwidth and plot
optimal_bandwidth <- h_seq[which.min(CV_err_h)]
png("/Users/veronica/Dropbox/Apps/Overleaf/EC_709_vcperez/PSET_1/Q2_1_bandwidth_selection_kernel.png")
plot(h_seq, CV_err_h, type = "b", lwd = 3, col = "blue",
     xlab = "Smoothing Bandwidth", ylab = "LOOCV Prediction Error",
     main = "Bandwidth Selection using Cross-Validation")
dev.off()

# Kernel regression with optimal bandwidth
Kreg_optimal <- ksmooth(x = X, y = Y, kernel = "normal", bandwidth = optimal_bandwidth)
png("/Users/veronica/Dropbox/Apps/Overleaf/EC_709_vcperez/PSET_1/Q2_1_optimal_kernel_regression.png")
plot(Kreg_optimal, col = "purple", lwd = 2, 
     main = "Optimal Kernel Regression of Log Gasoline Consumption on Log Price",
     xlab = "Log Price", ylab = "Log Gasoline Consumption")
dev.off()


# LOCALLY LINEAR ----------------------------------------------------------


library(KernSmooth)
library(ks)

library(locpol) # Ensure locpol is installed and loaded

# Cross-validation for Locally Linear Regression
CV_err_h_locpol <- numeric(length(h_seq))

for (j in seq_along(h_seq)) {
  h_using <- h_seq[j]
  CV_err <- numeric(n)
  
  for (i in 1:n) {
    # Leave-one-out sets
    X_val <- X[i]
    Y_val <- Y[i]
    X_tr <- X[-i]
    Y_tr <- Y[-i]
    
    # Locally linear regression on training set
    model <- locpoly(y = Y_tr, x = X_tr, kernel = "normal", bandwidth = h_using, degree = 1)
    
    # Predict on a grid of x-values
    x_grid <- seq(min(X_tr), max(X_tr), length.out = 100)
    fitted_vals <- locpoly(y = Y_tr, x = X_tr, kernel = "normal", bandwidth = h_using, degree = 1)
    
    # Extract fitted values from the model output
    fitted_vals <- fitted_vals$y  # Adjust according to the actual structure of `locpoly` output
    
    # Find the closest x_grid value to X_val
    closest_index <- which.min(abs(x_grid - X_val))
    Y_val_predict <- fitted_vals[closest_index]
    
    # Compute squared error
    CV_err[i] <- (Y_val - Y_val_predict)^2
  }
  
  # Average CV error for current bandwidth
  CV_err_h_locpol[j] <- mean(CV_err)
}

# Optimal bandwidth for locally linear regression
optimal_bandwidth_locpol <- h_seq[which.min(CV_err_h_locpol)]

# Plot the cross-validation error
png("/Users/veronica/Dropbox/Apps/Overleaf/EC_709_vcperez/PSET_1/Q2_1_locpol_bandwidth_selection.png")
plot(h_seq, CV_err_h_locpol, type = "b", lwd = 3, col = "blue",
     xlab = "Smoothing Bandwidth", ylab = "LOOCV Prediction Error",
     main = "Locally Linear Regression Bandwidth Selection")
dev.off()

# Plot with optimal bandwidth
model_optimal_locpol <- locpoly(y = Y, x = X, kernel = "normal", bandwidth = optimal_bandwidth_locpol, degree = 1)

png("/Users/veronica/Dropbox/Apps/Overleaf/EC_709_vcperez/PSET_1/Q2_1_locpol_optimal_regression.png")
plot(model_optimal_locpol, col = "purple", lwd = 2, 
     main = "Optimal Local Linear Regression of Log Gasoline Consumption on Log Price",
     xlab = "Log Price", ylab = "Log Gasoline Consumption")
dev.off()


# SERIES POWER ------------------------------------------------------------

# Cross-validation for Series Power Regression
degree_seq <- 1:10
CV_err_poly <- numeric(length(degree_seq))

for (d in degree_seq) {
  CV_err <- numeric(n)
  
  for (i in 1:n) {
    X_val <- X[i]
    Y_val <- Y[i]
    X_tr <- X[-i]
    Y_tr <- Y[-i]
    model <- lm(Y_tr ~ poly(X_tr, d))
    Y_val_predict <- predict(model, data.frame(X_tr = X_val))
    CV_err[i] <- (Y_val - Y_val_predict)^2
  }
  
  CV_err_poly[d] <- mean(CV_err)
}

# Optimal degree and plot
optimal_degree <- degree_seq[which.min(CV_err_poly)]
png("/Users/veronica/Dropbox/Apps/Overleaf/EC_709_vcperez/PSET_1/Q2_1_series_power_bandwidth_selection.png")
plot(degree_seq, CV_err_poly, type = "b", lwd = 3, col = "blue",
     xlab = "Degree of Polynomial", ylab = "LOOCV Prediction Error",
     main = "Series Power Regression Degree Selection")
dev.off()

# Series Power Regression with optimal degree
model_optimal_poly <- lm(Y ~ poly(X, optimal_degree))
x_grid <- seq(min(X), max(X), length.out = 100)
y_grid <- predict(model_optimal_poly, newdata = data.frame(X = x_grid))

png("/Users/veronica/Dropbox/Apps/Overleaf/EC_709_vcperez/PSET_1/Q2_2_series_power_optimal_regression.png")
plot(X, Y, type = "n", main = "Series Power Regression",
     xlab = "Log Price", ylab = "Log Gasoline Consumption")
lines(x_grid, y_grid, col = "purple", lwd = 2)
dev.off()



##############################################################
#                       B SPLINE
##############################################################

library(splines)

# Cross-validation for B-spline Regression
knots_seq <- 1:10
CV_err_bspline <- numeric(length(knots_seq))

for (k in knots_seq) {
  CV_err <- numeric(n)
  
  for (i in 1:n) {
    X_val <- X[i]
    Y_val <- Y[i]
    X_tr <- X[-i]
    Y_tr <- Y[-i]
    
    model <- lm(Y_tr ~ bs(X_tr, df = k))
    Y_val_predict <- predict(model, data.frame(X_tr = X_val))
    
    CV_err[i] <- (Y_val - Y_val_predict)^2
  }
  
  CV_err_bspline[k] <- mean(CV_err)
}

# Optimal number of knots and plot
optimal_knots <- knots_seq[which.min(CV_err_bspline)]
png("/Users/veronica/Dropbox/Apps/Overleaf/EC_709_vcperez/PSET_1/Q2_bspline_bandwidth_selection.png")
plot(knots_seq, CV_err_bspline, type = "b", lwd = 3, col = "blue",
     xlab = "Number of Knots", ylab = "LOOCV Prediction Error",
     main = "B-spline Regression Knots Selection")
dev.off

# B-spline regression with optimal number of knots
model_optimal_bspline <- lm(Y ~ bs(X, df = optimal_knots))
x_grid <- seq(min(X), max(X), length.out = 100)
y_grid <- predict(model_optimal_bspline, newdata = data.frame(X = x_grid))

png("/Users/veronica/Dropbox/Apps/Overleaf/EC_709_vcperez/PSET_1/Q2_bspline_optimal_regression.png")
plot(X, Y, type = "n", main = "B-spline Regression",
     xlab = "Log Price", ylab = "Log Gasoline Consumption")
lines(x_grid, y_grid, col = "purple", lwd = 2)
dev.off()

##############################################################
##############################################################
#                       Question 2_2  
##############################################################
##############################################################



# KERNEL REGRESSION -------------------------------------------------------

# Kernel regression with optimal bandwidth
Kreg_optimal <- ksmooth(x = X, y = Y, kernel = "normal", bandwidth = optimal_bandwidth)

# Create a function for interpolation
Kreg_fun <- approxfun(Kreg_optimal$x, Kreg_optimal$y, rule = 2)

# Estimate for price = 0.57
X_new <- log(0.57)
Y_estimate_kernel <- Kreg_fun(X_new)

# Standard error (assuming a simple approach)
residuals <- Y - Kreg_fun(X)
standard_error_kernel <- sd(residuals) / sqrt(length(X))

# Output results
cat("Kernel Regression Estimate at price 0.57:", Y_estimate_kernel, "\n")
cat("Kernel Regression Standard Error:", standard_error_kernel, "\n")



# LOCALLY LINEAR ----------------------------------------------------------

# Define the new data point
X_new <- log(0.57)

# Define x_grid for interpolation
x_grid <- seq(min(X), max(X), length.out = 401)

# Perform locally linear regression
model_optimal_locpol <- locpoly(y = Y, x = X, kernel = "normal", bandwidth = optimal_bandwidth_locpol, degree = 1)

# Extract fitted values for x_grid
fitted_vals <- approx(model_optimal_locpol$x, model_optimal_locpol$y, xout = x_grid)$y

# Interpolate to estimate the value at X_new
Y_estimate_locpol <- approx(x_grid, fitted_vals, xout = X_new)$y

# Compute standard error (assuming residuals from the whole dataset)
residuals <- Y - approx(x_grid, fitted_vals, xout = X)$y
standard_error_locpol <- sd(residuals) / sqrt(length(X))

# Output results
cat("Locally Linear Regression Estimate at price 0.57:", Y_estimate_locpol, "\n")
cat("Locally Linear Regression Standard Error:", standard_error_locpol, "\n")




# POWER REGRESSION --------------------------------------------------------

# Series Power Regression with optimal degree
model_optimal_poly <- lm(Y ~ poly(X, optimal_degree))

# Estimate for price = 0.57
X_new <- log(0.57)
Y_estimate_poly <- predict(model_optimal_poly, newdata = data.frame(X = X_new))

# Standard error estimation
# Calculating standard errors from the model
vcov_matrix <- vcov(model_optimal_poly)
standard_error_poly <- sqrt(diag(vcov_matrix))[2]  # Adjust index based on degree

# Output results
cat("Series Power Regression Estimate at price 0.57:", Y_estimate_poly, "\n")
cat("Series Power Regression Standard Error:", standard_error_poly, "\n")

# B SPLINE ----------------------------------------------------------------


# B-spline regression with optimal number of knots
model_optimal_bspline <- lm(Y ~ bs(X, df = optimal_knots))

# Estimate for price = 0.57
X_new <- log(0.57)
Y_estimate_bspline <- predict(model_optimal_bspline, newdata = data.frame(X = X_new))

# Standard error estimation
# Calculating standard errors from the model
vcov_matrix <- vcov(model_optimal_bspline)
standard_error_bspline <- sqrt(diag(vcov_matrix))[2]  # Adjust index based on df

# Output results
cat("B-Spline Regression Estimate at price 0.57:", Y_estimate_bspline, "\n")
cat("B-Spline Regression Standard Error:", standard_error_bspline, "\n")


# GRAPHS ------------------------------------------------------------------

library(ggplot2)

# Create a data frame with estimates and standard errors
results <- data.frame(
  Method = c("Kernel Regression", "Locally Linear Regression", "Series Power Regression", "B-Spline Regression"),
  Estimate = c(Y_estimate_kernel, Y_estimate_locpol, Y_estimate_poly, Y_estimate_bspline),
  StandardError = c(standard_error_kernel, standard_error_locpol, standard_error_poly, standard_error_bspline)
)

# Create a ggplot object
plot <- ggplot(results, aes(x = Method, y = Estimate, ymin = Estimate - StandardError, ymax = Estimate + StandardError)) +
  geom_point(size = 4) +
  geom_errorbar(width = 0.2) +
  labs(title = "Estimates with Standard Errors for Different Regression Methods",
       x = "Regression Method",
       y = "Estimate") +
  theme_minimal()

# Print the plot
print(plot)



###########################################################################
# QUESTION 2_3 ------------------------------------------------------------
###########################################################################

# Extract the variables
X1 <- log(data$price)
X2 <- log(data$income)
Y <- log(data$gas)

# Create the matrix for predictors
X <- cbind(X1, X2)


# Kernel Regression
bandwidths <- c(0.01, 0.05, 0.08, 0.1, 0.15, 0.2, 0.25, 0.3)

# Perform kernel regression and plot the results
png("/Users/veronica/Dropbox/Apps/Overleaf/EC_709_vcperez/PSET_1/Q2_3_kernel_regression_comparison.png")
plot(NULL, xlim = range(X), ylim = range(Y), 
     main = "Kernel Regression with Different Bandwidths",
     xlab = "Log Price", ylab = "Log Gasoline Consumption")

colors <- c("orange", "purple", "limegreen", "blue", "red", "yellow", "green", "cyan")
for (i in seq_along(bandwidths)) {
  Kreg <- ksmooth(x = X, y = Y, kernel = "normal", bandwidth = bandwidths[i])
  lines(Kreg, lwd = 2, col = colors[i])
}
legend("topright", legend = paste("h =", bandwidths), lwd = 2, col = colors)
dev.off()

# Cross-validation for bandwidth selection
n <- length(X)
h_seq <- seq(from = 0.01, to = 0.3, by = 0.01)
CV_err_h <- numeric(length(h_seq))

for (j in seq_along(h_seq)) {
  h_using <- h_seq[j]
  CV_err <- numeric(n)
  
  for (i in 1:n) {
    X_val <- X[i]
    Y_val <- Y[i]
    X_tr <- X[-i]
    Y_tr <- Y[-i]
    
    Y_val_predict <- ksmooth(x = X_tr, y = Y_tr, kernel = "normal", bandwidth = h_using, x.points = X_val)$y
    CV_err[i] <- (Y_val - Y_val_predict)^2
  }
  
  CV_err_h[j] <- mean(CV_err)
}

# Optimal bandwidth and plot
optimal_bandwidth <- h_seq[which.min(CV_err_h)]
png("/Users/veronica/Dropbox/Apps/Overleaf/EC_709_vcperez/PSET_1/Q2_3_bandwidth_selection_kernel.png")
plot(h_seq, CV_err_h, type = "b", lwd = 3, col = "blue",
     xlab = "Smoothing Bandwidth", ylab = "LOOCV Prediction Error",
     main = "Bandwidth Selection using Cross-Validation")
dev.off()

# Kernel regression with optimal bandwidth
Kreg_optimal <- ksmooth(x = X, y = Y, kernel = "normal", bandwidth = optimal_bandwidth)
png("/Users/veronica/Dropbox/Apps/Overleaf/EC_709_vcperez/PSET_1/Q2_3_optimal_kernel_regression.png")
plot(Kreg_optimal, col = "purple", lwd = 2, 
     main = "Optimal Kernel Regression of Log Gasoline Consumption on Log Price",
     xlab = "Log Price", ylab = "Log Gasoline Consumption")
dev.off()


# LOCALLY LINEAR ----------------------------------------------------------


library(KernSmooth)
library(ks)

library(locpol) # Ensure locpol is installed and loaded

# Cross-validation for Locally Linear Regression
CV_err_h_locpol <- numeric(length(h_seq))

for (j in seq_along(h_seq)) {
  h_using <- h_seq[j]
  CV_err <- numeric(n)
  
  for (i in 1:n) {
    # Leave-one-out sets
    X_val <- X[i]
    Y_val <- Y[i]
    X_tr <- X[-i]
    Y_tr <- Y[-i]
    
    # Locally linear regression on training set
    model <- locpoly(y = Y_tr, x = X_tr, kernel = "normal", bandwidth = h_using, degree = 1)
    
    # Predict on a grid of x-values
    x_grid <- seq(min(X_tr), max(X_tr), length.out = 100)
    fitted_vals <- locpoly(y = Y_tr, x = X_tr, kernel = "normal", bandwidth = h_using, degree = 1)
    
    # Extract fitted values from the model output
    fitted_vals <- fitted_vals$y  # Adjust according to the actual structure of `locpoly` output
    
    # Find the closest x_grid value to X_val
    closest_index <- which.min(abs(x_grid - X_val))
    Y_val_predict <- fitted_vals[closest_index]
    
    # Compute squared error
    CV_err[i] <- (Y_val - Y_val_predict)^2
  }
  
  # Average CV error for current bandwidth
  CV_err_h_locpol[j] <- mean(CV_err)
}

# Optimal bandwidth for locally linear regression
optimal_bandwidth_locpol <- h_seq[which.min(CV_err_h_locpol)]

# Plot the cross-validation error
png("/Users/veronica/Dropbox/Apps/Overleaf/EC_709_vcperez/PSET_1/Q2_3_locpol_bandwidth_selection.png")
plot(h_seq, CV_err_h_locpol, type = "b", lwd = 3, col = "blue",
     xlab = "Smoothing Bandwidth", ylab = "LOOCV Prediction Error",
     main = "Locally Linear Regression Bandwidth Selection")
dev.off()

# Plot with optimal bandwidth
model_optimal_locpol <- locpoly(y = Y, x = X, kernel = "normal", bandwidth = optimal_bandwidth_locpol, degree = 1)

png("/Users/veronica/Dropbox/Apps/Overleaf/EC_709_vcperez/PSET_1/Q2_3_locpol_optimal_regression.png")
plot(model_optimal_locpol, col = "purple", lwd = 2, 
     main = "Optimal Local Linear Regression of Log Gasoline Consumption on Log Price",
     xlab = "Log Price", ylab = "Log Gasoline Consumption")
dev.off()


# SERIES POWER ------------------------------------------------------------

# Cross-validation for Series Power Regression
degree_seq <- 1:10
CV_err_poly <- numeric(length(degree_seq))

for (d in degree_seq) {
  CV_err <- numeric(n)
  
  for (i in 1:n) {
    X_val <- X[i]
    Y_val <- Y[i]
    X_tr <- X[-i]
    Y_tr <- Y[-i]
    
    model <- lm(Y_tr ~ poly(X_tr, d))
    Y_val_predict <- predict(model, data.frame(X_tr = X_val))
    
    CV_err[i] <- (Y_val - Y_val_predict)^2
  }
  
  CV_err_poly[d] <- mean(CV_err)
}

# Optimal degree and plot
optimal_degree <- degree_seq[which.min(CV_err_poly)]
png("/Users/veronica/Dropbox/Apps/Overleaf/EC_709_vcperez/PSET_1/Q2_3_series_power_bandwidth_selection.png")
plot(degree_seq, CV_err_poly, type = "b", lwd = 3, col = "blue",
     xlab = "Degree of Polynomial", ylab = "LOOCV Prediction Error",
     main = "Series Power Regression Degree Selection")
dev.off()

# Series Power Regression with optimal degree
model_optimal_poly <- lm(Y ~ poly(X, optimal_degree))
x_grid <- seq(min(X), max(X), length.out = 100)
y_grid <- predict(model_optimal_poly, newdata = data.frame(X = x_grid))

png("/Users/veronica/Dropbox/Apps/Overleaf/EC_709_vcperez/PSET_1/Q2_3_series_power_optimal_regression.png")
plot(X, Y, type = "n", main = "Series Power Regression",
     xlab = "Log Price", ylab = "Log Gasoline Consumption")
lines(x_grid, y_grid, col = "purple", lwd = 2)
dev.off()



##############################################################
#                       B SPLINE
##############################################################

library(splines)

# Cross-validation for B-spline Regression
knots_seq <- 1:10
CV_err_bspline <- numeric(length(knots_seq))

for (k in knots_seq) {
  CV_err <- numeric(n)
  
  for (i in 1:n) {
    X_val <- X[i]
    Y_val <- Y[i]
    X_tr <- X[-i]
    Y_tr <- Y[-i]
    
    model <- lm(Y_tr ~ bs(X_tr, df = k))
    Y_val_predict <- predict(model, data.frame(X_tr = X_val))
    
    CV_err[i] <- (Y_val - Y_val_predict)^2
  }
  
  CV_err_bspline[k] <- mean(CV_err)
}

# Optimal number of knots and plot
optimal_knots <- knots_seq[which.min(CV_err_bspline)]
png("/Users/veronica/Dropbox/Apps/Overleaf/EC_709_vcperez/PSET_1/Q2_3_bspline_bandwidth_selection.png")
plot(knots_seq, CV_err_bspline, type = "b", lwd = 3, col = "blue",
     xlab = "Number of Knots", ylab = "LOOCV Prediction Error",
     main = "B-spline Regression Knots Selection")
dev.off

# B-spline regression with optimal number of knots
model_optimal_bspline <- lm(Y ~ bs(X, df = optimal_knots))
x_grid <- seq(min(X), max(X), length.out = 100)
y_grid <- predict(model_optimal_bspline, newdata = data.frame(X = x_grid))

png("/Users/veronica/Dropbox/Apps/Overleaf/EC_709_vcperez/PSET_1/Q2_3_bspline_optimal_regression.png")
plot(X, Y, type = "n", main = "B-spline Regression",
     xlab = "Log Price", ylab = "Log Gasoline Consumption")
lines(x_grid, y_grid, col = "purple", lwd = 2)
dev.off()



# QUESTION 2_4 ------------------------------------------------------------

# Load necessary libraries
library(hdm)

# Load and process the data
data <- read.dta("/Users/veronica/Dropbox/PhD/2024_2/EC_709/PSET1/yn.dta")

# Create log-transformed variables and interaction terms if needed
data$log_gas <- log(data$gas)
data$log_price <- log(data$price)
data$log_income <- log(data$income)

# Extract the variables
X <- data[, c("log_price", "log_income", "age", "driver", "hhsize", "month", 
              "fueltype", "urban", "prov", "year", "distance", "youngsingle")]
Y <- data$log_gas


# LASSO ALL VARIABLES -----------------------------------------------------

# Fit Lasso model with all variables
lasso_model_all = rlasso(data$log_gas ~ log_price + age + driver + hhsize + fueltype + urban + prov + year, data , post = FALSE)  # use lasso, not-Post-lasso
# lasso.reg = rlasso(X, Y, post=FALSE)
summary_lasso <- summary(lasso_model_all, all = FALSE)  # can also do print(lasso.reg, all=FALSE)

print(summary_lasso)

# LASSO 2 WAY INTERACTIONS ------------------------------------------------

# Prepare the data frame for interactions
data_for_interactions <- data[, c("log_price", "log_income", "age", "driver", "hhsize", 
                                  "month", "fueltype", "urban", "prov", "year", 
                                  "distance", "youngsingle")]
data_for_interactions$log_gas <- data$log_gas

# Create interaction terms
data_interactions <- model.matrix(~ .^2 - 1, data = data_for_interactions)

# Extract response
Y_interactions <- data_interactions[, "log_gas"]

# Remove the response column from the predictors matrix
X_interactions <- data_interactions[, -which(colnames(data_interactions) == "log_gas")]

# Fit Lasso model with interactions
lasso_model_interactions <- rlasso(X_interactions, Y_interactions, post = FALSE)

# Display the summary of the Lasso model with interactions
summary_lasso_interactions <- summary(lasso_model_interactions, all = FALSE)
print(summary_lasso_interactions)
