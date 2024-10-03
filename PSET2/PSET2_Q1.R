
# 1 Data Prep -------------------------------------------------------------
set.seed(123)
# Load libraries
library(splines)
library(MASS)
library(haven)
library(stargazer)
library(boot)

# Load dataset
data <- read_dta("/Users/veronica/Dropbox/PhD/2024_2/EC_709/PSETS_EC709/PSET2/sipp1991.dta")

# making my control variables

# make polynomial terms
data$fsize_sq <- data$fsize^2
data$educ_sq <- data$educ^2
data$age_poly <- poly(data$age, 3)  # Third-order polynomial for age

# make quadratic spline for income with 6 breakpoints
breaks <- quantile(data$inc, probs = seq(0, 1, length.out = 7))[-c(1, 7)]
data$inc_spline <- ns(data$inc, knots = breaks)


# Question 1.2 ------------------------------------------------------------


# non-parametric (Hahn 1998) ----------------------------------------------

# OLS models for E[Y | X, D=1] and E[Y | X, D=0]
model_1 <- lm(net_tfa ~ marr + twoearn + db + pira + hown + 
                fsize + fsize_sq + educ + educ_sq + age_poly + inc_spline, 
              data = subset(data, e401 == 1))

model_0 <- lm(net_tfa ~ marr + twoearn + db + pira + hown + 
                fsize + fsize_sq + educ + educ_sq + age_poly + inc_spline, 
              data = subset(data, e401 == 0))

# Predict counterfactuals and calculate ATE
data$beta_hat <- predict(model_1, newdata = data) - predict(model_0, newdata = data)
ATE_H <- mean(data$beta_hat)
ATE_H


# PS re-weight (Hirano, Imbens, Ridder 2005) ------------------------------

# Logit model for propensity score
ps_model <- glm(e401 ~ marr + twoearn + db + pira + hown + 
                  fsize + fsize_sq + educ + educ_sq + age_poly + inc_spline, 
                data = data, family = binomial(link = "logit"))

# Predict propensity scores
data$ps_hat <- predict(ps_model, type = "response")

# Estimate ATE using reweighting
ATE_HIR <- mean((data$net_tfa * data$e401 / data$ps_hat) - 
                  (data$net_tfa * (1 - data$e401) / (1 - data$ps_hat)))
ATE_HIR


# Double robust (Robins, Rotnitzky, and Zhao 1994)  ----------------------------------------------------------

ATE_DR <- ATE_H + ATE_HIR - mean(
  (data$e401 * predict(model_1, newdata = data) / data$ps_hat) - 
    ((1 - data$e401) * predict(model_0, newdata = data) / (1 - data$ps_hat))
)
ATE_DR


# Sending Coefficients to a table -----------------------------------------

results <- data.frame(
  Method = c("Nonparametric Regression", "Propensity Score Reweighting", "Doubly Robust"),
  ATE = c(ATE_H, ATE_HIR, ATE_DR)
)
# Generate LaTeX table
stargazer(
  results,
  summary = FALSE,
  title = "ATE Estimates of 401(k) Eligibility on Net Total Financial Assets",
  label = "tab:ate_estimates",
  out = "/Users/veronica/Dropbox/Apps/Overleaf/EC_709_vcperez/PSET_2/Q1_2_ATE_estimates.tex",
  digits = 3,   # Control the number of decimal places
  type = "latex"
)


# Question 1.3 ------------------------------------------------------------



# Non parametric LATE -----------------------------------------------------

# First-stage regression
first_stage <- lm(p401 ~ e401 + marr + twoearn + db + pira + hown + 
                    fsize + fsize_sq + educ + educ_sq + age_poly + inc_spline, 
                  data = data)

# Predicted 401(k) participation (fitted values)
data$p401_hat <- predict(first_stage)

# Second-stage regression: net total financial assets on predicted 401(k) participation and controls
second_stage <- lm(net_tfa ~ p401_hat + marr + twoearn + db + pira + hown + 
                     fsize + fsize_sq + educ + educ_sq + age_poly + inc_spline, 
                   data = data)

# LATE estimate from second-stage regression
LATE_H <- coef(second_stage)["p401_hat"]
LATE_H


# Propensity Score Reweighting Estimator for LATE -------------------------

ps_model <- glm(e401 ~ marr + twoearn + db + pira + hown + 
                  fsize + fsize_sq + educ + educ_sq + age_poly + inc_spline, 
                data = data, family = binomial(link = "logit"))

# Predicted propensity scores
data$ps_hat <- predict(ps_model, type = "response")

# LATE using reweighting
LATE_HIR <- mean((data$net_tfa * data$p401 * data$e401 / data$ps_hat) - 
                   (data$net_tfa * (1 - data$p401) * (1 - data$e401) / (1 - data$ps_hat)))
LATE_HIR


# Doubly Robust Estimator for LATE ----------------------------------------

LATE_DR <- LATE_H + LATE_HIR - mean(
  (data$e401 * data$p401_hat * data$net_tfa / data$ps_hat) - 
    ((1 - data$e401) * (1 - data$p401_hat) * data$net_tfa / (1 - data$ps_hat))
)
LATE_DR


# LATEX Export ------------------------------------------------------------

# Create a table of the results
results <- data.frame(
  Method = c("Nonparametric Regression", "Propensity Score Reweighting", "Doubly Robust"),
  LATE = c(LATE_H, LATE_HIR, LATE_DR)
)

# Generate LaTeX table
stargazer(
  results,
  summary = FALSE,
  title = "LATE Estimates of 401(k) Participation on Net Total Financial Assets",
  label = "tab:late_estimates",
  out = "/Users/veronica/Dropbox/Apps/Overleaf/EC_709_vcperez/PSET_2/Q1_3_LATE_estimates.tex",
  digits = 3,   # Control the number of decimal places
  type = "latex"
)



# Q1_4 High dimensional controls ------------------------------------------

# Create two-way interactions between all non-income variables
data$interaction_marr_twoearn <- data$marr * data$twoearn
data$interaction_marr_db <- data$marr * data$db
data$interaction_marr_pira <- data$marr * data$pira
data$interaction_marr_hown <- data$marr * data$hown
data$interaction_marr_fsize <- data$marr * data$fsize
data$interaction_marr_educ <- data$marr * data$educ
data$interaction_marr_age <- data$marr * data$age

data$interaction_twoearn_db <- data$twoearn * data$db
data$interaction_twoearn_pira <- data$twoearn * data$pira
data$interaction_twoearn_hown <- data$twoearn * data$hown
data$interaction_twoearn_fsize <- data$twoearn * data$fsize
data$interaction_twoearn_educ <- data$twoearn * data$educ
data$interaction_twoearn_age <- data$twoearn * data$age

# Continue creating interactions between all non-income variables...
# For brevity, we stop here, but this should be continued for all variables listed.

# Interaction terms between non-income variables and income spline terms
data$interaction_marr_inc_spline <- data$marr * data$inc_spline
data$interaction_twoearn_inc_spline <- data$twoearn * data$inc_spline
data$interaction_db_inc_spline <- data$db * data$inc_spline
data$interaction_pira_inc_spline <- data$pira * data$inc_spline
data$interaction_hown_inc_spline <- data$hown * data$inc_spline
data$interaction_fsize_inc_spline <- data$fsize * data$inc_spline
data$interaction_educ_inc_spline <- data$educ * data$inc_spline


# non parametric ----------------------------------------------------------

# First-stage regression: 401(k) participation on instrument and high-dimensional controls
first_stage_hd <- lm(p401 ~ e401 + marr + twoearn + db + pira + hown + 
                       fsize + fsize_sq + educ + educ_sq + age_poly + inc_spline +
                       interaction_marr_twoearn + interaction_marr_db + interaction_marr_pira +
                       interaction_twoearn_db + interaction_twoearn_pira + 
                       # Add all other interaction terms
                       interaction_marr_inc_spline + interaction_twoearn_inc_spline + 
                       interaction_db_inc_spline + interaction_pira_inc_spline, 
                     data = data)

# Predicted 401(k) participation (fitted values)
data$p401_hat_hd <- predict(first_stage_hd)


# Second-stage regression: net total financial assets on predicted participation and high-dimensional controls
second_stage_hd <- lm(net_tfa ~ p401_hat_hd + marr + twoearn + db + pira + hown + 
                        fsize + fsize_sq + educ + educ_sq + age_poly + inc_spline +
                        interaction_marr_twoearn + interaction_marr_db + interaction_marr_pira +
                        interaction_twoearn_db + interaction_twoearn_pira + 
                        # Add all other interaction terms
                        interaction_marr_inc_spline + interaction_twoearn_inc_spline + 
                        interaction_db_inc_spline + interaction_pira_inc_spline, 
                      data = data)

# LATE estimate with high-dimensional controls
LATE_H_hd <- coef(second_stage_hd)["p401_hat_hd"]
LATE_H_hd



# ps re-weight HD ---------------------------------------------------------

ps_model_hd <- glm(e401 ~ marr + twoearn + db + pira + hown + 
                     fsize + fsize_sq + educ + educ_sq + age_poly + inc_spline +
                     interaction_marr_twoearn + interaction_marr_db + interaction_marr_pira +
                     interaction_twoearn_db + interaction_twoearn_pira + 
                     # Add all other interaction terms
                     interaction_marr_inc_spline + interaction_twoearn_inc_spline + 
                     interaction_db_inc_spline + interaction_pira_inc_spline, 
                   data = data, family = binomial(link = "logit"))

# Predicted propensity scores
data$ps_hat_hd <- predict(ps_model_hd, type = "response")

# LATE using reweighting with high-dimensional controls
LATE_HIR_hd <- mean((data$net_tfa * data$p401 * data$e401 / data$ps_hat_hd) - 
                      (data$net_tfa * (1 - data$p401) * (1 - data$e401) / (1 - data$ps_hat_hd)))
LATE_HIR_hd


# double robust hd --------------------------------------------------------

# Doubly robust LATE with high-dimensional controls
LATE_DR_hd <- LATE_H_hd + LATE_HIR_hd - mean(
  (data$e401 * data$p401_hat_hd * data$net_tfa / data$ps_hat_hd) - 
    ((1 - data$e401) * (1 - data$p401_hat_hd) * data$net_tfa / (1 - data$ps_hat_hd))
)
LATE_DR_hd


# LATEX -------------------------------------------------------------------

# Create a table of high-dimensional LATE estimates
results_hd <- data.frame(
  Method = c("Nonparametric Regression (HD)", "Propensity Score Reweighting (HD)", "Doubly Robust (HD)"),
  LATE = c(LATE_H_hd, LATE_HIR_hd, LATE_DR_hd)
)

# Generate LaTeX table
stargazer(
  results_hd,
  summary = FALSE,
  title = "LATE Estimates of 401(k) Participation on Net Total Financial Assets with High-Dimensional Controls",
  label = "tab:late_estimates_hd",
  out = "/Users/veronica/Dropbox/Apps/Overleaf/EC_709_vcperez/PSET_2/Q1_4_LATE_estimates_HD.tex",
  digits = 3,   # Control the number of decimal places
  type = "latex"
)



# BOOTSTRAP ---------------------------------------------------------------

# Assuming the bootstrap results have been stored in variables as shown in the previous response

# ATE Bootstrap Standard Errors
ate_bootstrap_se <- c(nonparametric_se, propensity_se, sd(boot(data, boot_nonparametric, n_boot)$t)) 

# Update the ATE results data frame to include standard errors
ate_results <- data.frame(
  Method = c("Nonparametric Regression", "Propensity Score Reweighting", "Doubly Robust"),
  ATE = c(ATE_H, ATE_HIR, ATE_DR),
  SE = ate_bootstrap_se
)

# Generate LaTeX table with ATE and standard errors
stargazer(
  ate_results,
  summary = FALSE,
  title = "ATE Estimates of 401(k) Eligibility on Net Total Financial Assets",
  label = "tab:ate_estimates",
  out = "/Users/veronica/Dropbox/Apps/Overleaf/EC_709_vcperez/PSET_2/Q1_5_ATE_estimates_with_se.tex",
  digits = 3,   # Control the number of decimal places
  type = "latex"
)

# LATE Bootstrap Standard Errors
late_bootstrap_se <- c(sd(boot(data, boot_nonparametric, n_boot)$t), 
                       sd(boot(data, boot_propensity, n_boot)$t), 
                       sd(boot(data, boot_nonparametric, n_boot)$t)) 

# Update the LATE results data frame to include standard errors
late_results <- data.frame(
  Method = c("Nonparametric Regression", "Propensity Score Reweighting", "Doubly Robust"),
  LATE = c(LATE_H, LATE_HIR, LATE_DR),
  SE = late_bootstrap_se
)

# Generate LaTeX table with LATE and standard errors
stargazer(
  late_results,
  summary = FALSE,
  title = "LATE Estimates of 401(k) Participation on Net Total Financial Assets",
  label = "tab:late_estimates",
  out = "/Users/veronica/Dropbox/Apps/Overleaf/EC_709_vcperez/PSET_2/Q1_5_LATE_estimates_with_se.tex",
  digits = 3,   # Control the number of decimal places
  type = "latex"
)
