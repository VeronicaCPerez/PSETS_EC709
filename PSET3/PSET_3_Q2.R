library(haven)
library(quantreg)
library(boot)
library(reshape2)

# Load the custom IVQR function
source("/Users/veronica/Dropbox/PhD/2024_2/EC_709/PSETS_EC709/PSET3/IVQR.R")

# Load the data
data <- read_dta("/Users/veronica/Dropbox/PhD/2024_2/EC_709/PSETS_EC709/PSET3/jtpa.dta")

# Filter for men only
data_men <- subset(data, male == 1)


# Question 1 --------------------------------------------------------------

taus <- seq(0.1, 0.9, by = 0.1)

# Quantile regression without covariates
qreg_no_cov <- rq(earnings ~ trained, data = data_men, tau = taus)

# Plot results
plot(taus, coef(qreg_no_cov)[2,], type = "b", col = "blue", 
     xlab = "Quantiles", ylab = "JTPA Effect on Earnings",
     main = "Quantile Effects of JTPA without Covariates")

# Quantile regression with covariates (e.g., age, education, etc.)
qreg_with_cov <- rq(earnings ~ trained + age2225 + age2629 +
                    age3035 + age3644 + age4554 + 
                    black + hispanic + married + wkless13, 
                    data = data_men, tau = taus)

# Plot with covariates
points(taus, coef(qreg_with_cov)[2,], type = "b", col = "red")
legend("topright", legend = c("No Covariates", "With Covariates"),
       col = c("blue", "red"), lty = 1, pch = 1)


# Q2 ----------------------------------------------------------------------

# Check for missing values
print(sum(is.na(data_men)))  # Count missing values in the dataset

# Check for unique values in the 'offer' variable
print(unique(data_men$offer))  # Inspect unique values of the instrument

# Check for multicollinearity
cor_matrix <- cor(data_men[, c("trained", "earnings", "offer")], use = "pairwise.complete.obs")
print(cor_matrix)  # Print the correlation matrix

# Set quantiles to estimate
taus <- seq(0.1, 0.9, by = 0.1)

# Initialize vectors to store results
ivqr_results <- numeric(length(taus))

# Loop through each quantile
for (i in seq_along(taus)) {
  tau <- taus[i]
  
  # Step 1: Fit the quantile regression using the instrument
  # First stage: predict 'trained' from 'offer'
  first_stage <- rq(trained ~ offer, data = data_men, tau = tau)
  
  # Check if first stage is successful
  if (!is.null(first_stage)) {
    # Get fitted values from the first stage
    data_men$fitted_values <- fitted(first_stage)
    
    # Step 2: Fit the second stage regression using the fitted values
    # Regression of earnings on fitted values
    ivqr_model <- rq(earnings ~ fitted_values, data = data_men, tau = tau)
    
    # Extract the coefficient for the fitted values
    ivqr_results[i] <- coef(ivqr_model)["fitted_values"]
  } else {
    ivqr_results[i] <- NA  # If the first stage fails, store NA
  }
}

# Output the results
print(ivqr_results)

# Prepare to plot the results
plot(taus, ivqr_results, type = "b", col = "blue", 
     xlab = "Quantiles", ylab = "IVQR Effect on Earnings",
     main = "IVQR Estimates for JTPA on Earnings (Males)")


# Q3 ----------------------------------------------------------------------


# Define a function for IVQR estimation for bootstrap
ivqr_fn <- function(data, indices) {
  data_boot <- data[indices, ]
  ivqr_est <- ivqr(earnings ~ trained + age + education + exper + exper2 | offer, 
                   data = data_boot, taus = taus)
  return(coef(ivqr_est)[2,])
}

# Run bootstrap
boot_results <- boot(data_men, ivqr_fn, R = 200)

# Get 95% confidence intervals
conf_intervals <- boot.ci(boot_results, type = "perc")

# Add confidence intervals to the plot
plot(taus, coef(ivqr_with_cov)[2,], type = "b", col = "red", 
     ylim = range(c(conf_intervals$perc[,4], conf_intervals$perc[,5])),
     xlab = "Quantiles", ylab = "JTPA Effect on Earnings",
     main = "IVQR with Confidence Intervals")

arrows(taus, conf_intervals$perc[,4], taus, conf_intervals$perc[,5], 
       length = 0.05, angle = 90, code = 3, col = "gray")



# Q4 ----------------------------------------------------------------------


taus <- seq(0.1, 0.9, by = 0.1)

# Instrumental quantile regression for TOT (Local QTE)
lqte_no_cov <- sapply(taus, function(tau) {
  ivqr_fit <- rq(earnings ~ trained | offer, data = data_men, tau = tau)
  return(coef(ivqr_fit)[2])
})

# Plot LQTE estimates without covariates
plot(taus, lqte_no_cov, type = "b", col = "green", 
     xlab = "Quantiles", ylab = "LQTE of JTPA on Earnings",
     main = "LQTE Quantile Effects without Covariates")


# Q5 ----------------------------------------------------------------------

points(taus, coef(qreg_with_cov)[2,], type = "b", col = "red")

# Create a data frame for ggplot
plot(taus, lqte_no_cov, type = "b", col = "green", 
     xlab = "Quantiles", ylab = "LQTE of JTPA on Earnings",
     main = "LQTE Quantile Effects without Covariates")



# Add points or lines for other comparisons if needed
# For example, you can add the quantile regression without IV
points(taus, qr_with_cov, type = "b", col = "red", pch = 17)

# Add a legend
legend("topright", legend = c("IVQR (Instrumented)", "Quantile Regression"),
       col = c("blue", "red"), pch = c(19, 17), lty = c(1, 1))

