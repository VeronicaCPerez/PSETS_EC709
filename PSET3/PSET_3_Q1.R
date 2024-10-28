# Load necessary libraries
library(foreign)  
library(quantreg) 
library(boot)     
library(ggplot2)  


# paths
wp <- "/Users/veronica/Dropbox/PhD/2024_2/EC_709/PSETS_EC709/PSET3"

# Load the data
data <- read.dta(paste0(wp, "/census00.dta"))

# Inspect the data
head(data)


# Question 1 --------------------------------------------------------------

# Set quantile indices
taus <- seq(0.05, 0.95, by = 0.05)

# Run quantile regressions and store the results
quantile_results <- lapply(taus, function(tau) {
  rq(logwk ~ educ + race + exper + I(exper^2), tau = tau, data = data)
})

# Extract coefficients and standard errors for 'educ'
beta_hat <- sapply(quantile_results, function(model) coef(model)["educ"])
std_errors <- sapply(quantile_results, function(model) summary(model)$coefficients["educ", 2])

# Inspect the results
beta_hat
std_errors


# Question 2 --------------------------------------------------------------

# Load additional library for bootstrap
library(boot)

# Function to compute the maximal t-statistic for a bootstrap sample
bootstrap_max_t_stat <- function(data, indices) {
  # Create a bootstrap sample
  data_boot <- data[indices, ]
  
  # Run quantile regressions for the bootstrap sample
  quantile_results_boot <- lapply(taus, function(tau) {
    rq(logwk ~ educ + race + exper + I(exper^2), tau = tau, data = data_boot)
  })
  
  # Extract the bootstrap coefficients for 'educ'
  beta_boot <- sapply(quantile_results_boot, function(model) coef(model)["educ"])
  
  # Compute the t-statistic for each quantile
  t_stat <- abs((beta_boot - beta_hat) / std_errors)
  
  # Return the maximal t-statistic
  return(max(t_stat))
}

# Bootstrap the maximal t-statistic 200 times
set.seed(123)  # for reproducibility
n_boot <- 50

boot_results <- boot(data = data, statistic = bootstrap_max_t_stat, R = n_boot, parallel = "multicore")

# Estimate the 95th percentile of the maximal t-statistic
c_95 <- quantile(boot_results$t, 0.95)

# Inspect the result
c_95

# Confidence bands
conf_band_lower <- beta_hat - c_95 * std_errors
conf_band_upper <- beta_hat + c_95 * std_errors

# Inspect the confidence band
conf_band_lower
conf_band_upper


# Quesiton 3 --------------------------------------------------------------

# Load ggplot2 for plotting
library(ggplot2)

# Create a data frame for plotting
plot_data <- data.frame(
  tau = taus,
  beta = beta_hat,
  lower = conf_band_lower,
  upper = conf_band_upper
)

# Plot using ggplot2
ggplot(plot_data, aes(x = tau)) +
  geom_line(aes(y = beta), color = "blue") +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, fill = "blue") +
  labs(title = "Quantile Regression Coefficient of Years of Schooling",
       x = "Quantile",
       y = "Coefficient of Years of Schooling (educ)") +
  theme_minimal()


# Question 4 --------------------------------------------------------------

# Load the 1980 and 1990 data and repeat the same steps

# Load datasets
data_2000 <- read.dta(paste0(wp, "/census00.dta"))
data_1980 <- read.dta(paste0(wp, "/census80.dta"))
data_1990 <- read.dta(paste0(wp, "/census90.dta"))

# Set quantile indices
taus <- seq(0.05, 0.95, by = 0.05)

# Function to run quantile regressions and extract beta and standard errors for 'educ'
run_quantile_regressions <- function(data) {
  quantile_results <- lapply(taus, function(tau) {
    rq(logwk ~ educ + race + exper + I(exper^2), tau = tau, data = data)
  })
  beta_hat <- sapply(quantile_results, function(model) coef(model)["educ"])
  std_errors <- sapply(quantile_results, function(model) summary(model)$coefficients["educ", 2])
  list(beta_hat = beta_hat, std_errors = std_errors)
}

# Run quantile regressions for each year
results_2000 <- run_quantile_regressions(data_2000)
results_1980 <- run_quantile_regressions(data_1980)
results_1990 <- run_quantile_regressions(data_1990)

# Function to compute the maximal t-statistic for a bootstrap sample
bootstrap_max_t_stat <- function(data, beta_hat, std_errors) {
  function(data, indices) {
    data_boot <- data[indices, ]
    quantile_results_boot <- lapply(taus, function(tau) {
      rq(logwk ~ educ + race + exper + I(exper^2), tau = tau, data = data_boot)
    })
    beta_boot <- sapply(quantile_results_boot, function(model) coef(model)["educ"])
    t_stat <- abs((beta_boot - beta_hat) / std_errors)
    return(max(t_stat))
  }
}

# Function to bootstrap and compute the confidence interval for a year
bootstrap_conf_band <- function(data, beta_hat, std_errors) {
  n_boot <- 20
  boot_results <- boot(data = data, statistic = bootstrap_max_t_stat(data, beta_hat, std_errors), R = n_boot, parallel = "multicore")
  c_95 <- quantile(boot_results$t, 0.95)
  conf_band_lower <- beta_hat - c_95 * std_errors
  conf_band_upper <- beta_hat + c_95 * std_errors
  list(lower = conf_band_lower, upper = conf_band_upper)
}

# Bootstrap for confidence intervals for each year
conf_band_2000 <- bootstrap_conf_band(data_2000, results_2000$beta_hat, results_2000$std_errors)
conf_band_1980 <- bootstrap_conf_band(data_1980, results_1980$beta_hat, results_1980$std_errors)
conf_band_1990 <- bootstrap_conf_band(data_1990, results_1990$beta_hat, results_1990$std_errors)

# Create data frames for plotting
plot_data_2000 <- data.frame(
  tau = taus,
  beta = results_2000$beta_hat,
  lower = conf_band_2000$lower,
  upper = conf_band_2000$upper,
  year = "2000"
)

plot_data_1980 <- data.frame(
  tau = taus,
  beta = results_1980$beta_hat,
  lower = conf_band_1980$lower,
  upper = conf_band_1980$upper,
  year = "1980"
)

plot_data_1990 <- data.frame(
  tau = taus,
  beta = results_1990$beta_hat,
  lower = conf_band_1990$lower,
  upper = conf_band_1990$upper,
  year = "1990"
)

# Combine all data for plotting
plot_data <- rbind(plot_data_2000, plot_data_1980, plot_data_1990)

# Plot using ggplot2
ggplot(plot_data, aes(x = tau, y = beta, color = year)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = year), alpha = 0.2) +
  labs(title = "Quantile Regression Coefficient of Years of Schooling (1980, 1990, 2000)",
       x = "Quantile",
       y = "Coefficient of Years of Schooling (educ)") +
  theme_minimal() +
  scale_color_manual(values = c("blue", "green", "red")) +
  scale_fill_manual(values = c("blue", "green", "red"))

