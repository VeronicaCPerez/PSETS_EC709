library(dplyr)
library(tidyr)
library(kableExtra)
library(data.table)
library(xtable)
library(ggplot2)
library(fixest)
library(did)
library(haven)
library(stargazer)
library(modelsummary)
library(tinytable)

## Fucntion for stars

get_stars <- function(p_value) {
  if (p_value < 0.01) {
    return("***")
  } else if (p_value < 0.05) {
    return("**")
  } else if (p_value < 0.1) {
    return("*")
  } else {
    return("")
  }
}

# so i don't get scientific numbers

options(scipen = 999)

# outputs path

outputs_path <- "/Users/veronica/Dropbox/Apps/Overleaf/EC_709_vcperez/PSET_4"

# Load the data
data <- read_dta("/Users/veronica/Dropbox/PhD/2024_2/EC_709/PSETS_EC709/PSET4/castle.dta")

# use data.table
setDT(data)

# make Y_it ln(homicides per 100000)
data$ln_homicide_per_100k <- log(data$homicide)


# Question 1 --------------------------------------------------------------

# List of columns to aggregate
col_to_aggregate <- c('ln_homicide_per_100k', 'homicide', 'robbery', 'assault', 
                      'burglary', 'larceny', 'motor', 'murder', 'police', 
                      'unemployrt', 'income', 'blackm_15_24', 'whitem_15_24', 
                      'blackm_25_44', 'whitem_25_44', 'prisoner', 'poverty', 
                      'exp_subsidy', 'exp_pubwelfare', 'popwt')

# Create dummy variable (=1 if cdl > 0)
data$treatment <- ifelse(data$cdl >0, 1, 0)

# Corrected code
summary_table <- rbindlist(
  lapply(col_to_aggregate, function(var) {
    data[, list(
      variable = var,
      mean_D_it_0 = mean(get(var)[treatment == 0], na.rm = TRUE),
      sd_D_it_0 = sd(get(var)[treatment == 0], na.rm = TRUE),
      mean_D_it_1 = mean(get(var)[treatment == 1], na.rm = TRUE),
      sd_D_it_1 = sd(get(var)[treatment == 1], na.rm = TRUE)
    )]
  })
)

# table for overleaf
table_content_Q1 <- print(
  xtable(summary_table),
  include.rownames = FALSE,
  include.colnames = FALSE,
  only.contents = TRUE,
  print.results = FALSE
)

writeLines(table_content_Q1, paste0(outputs_path,"/Pset4_Q1.tex"))

# Q2 ----------------------------------------------------------------------

# make never treated and treated groups

# make maximum for the D_it variable by state
data[, max_treatment := max(treatment), by = state]

# make group "never_treated" = 1 if max_treatment == 0
data$never_treated <- ifelse(data$max_treatment == 0, 1, 0)

# making data for graph for ln_homicide_per_100k between 1 and 0 d
data_did_graph <- data %>%
  group_by(year, never_treated) %>%
  summarise(mean_homicide = weighted.mean(homicide, popwt, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(group = ifelse(never_treated == 1, "Never Treated", "Treated"))

# keep years < 2005
data_did_graph <- data_did_graph %>%
                  subset(year < 2005)

# Plotting the time trend by group
paralel_trend <- ggplot(data_did_graph, aes(x = year, y = mean_homicide, color = group, group = group)) +
  geom_line() +
  geom_point() +
  labs(title = "Pararel trend for homicides, weighted by population",
       x = "Year",
       y = "homicides per 100k",
       color = "Group") +
  theme_minimal()

ggsave(paste0(outputs_path,"/PSET4_Q2.png"), plot = paralel_trend, dpi = 300)


# Q3 ----------------------------------------------------------------------


# TWFE estimation

TWFE_basic <- feols(ln_homicide_per_100k ~ treatment | state + year,data,weights = ~ popwt)

TWFE_covariates <- feols(ln_homicide_per_100k ~ treatment + robbery + assault +
                           burglary + larceny + motor + murder + police + unemployrt +
                           income + blackm_15_24 + whitem_15_24 + blackm_25_44 +
                           whitem_25_44 + prisoner + poverty +
                           exp_subsidy + exp_pubwelfare 
                           | state + year,data,weights = ~ popwt)


TWFE_covariates_2 <- feols(ln_homicide_per_100k ~ larceny + poverty | state + year,data,weights = ~ popwt)
# Create LaTeX table for the TWFE Basic model

# Fit models (assuming your model objects TWFE_basic and TWFE_covariates are already created)

# Extract model summary details
summary_basic <- summary(TWFE_basic)

# Extract coefficients, standard errors, t-values, and p-values
coef_table_basic <- summary_basic$coeftable
coefs_basic <- coef_table_basic[,1]
std_error_basic <- coef_table_basic[,2]
p_value_basic <- coef_table_basic[,4]
  
# Additional model info
obs_basic <- summary_basic$nobs
adj_r2_basic <- r2(TWFE_basic,'ar2')


### Summary covariates model
# Extract model summary details
summary_covs <- summary(TWFE_covariates)

# Extract coefficients, standard errors, t-values, and p-values
coef_table_covs <- summary_covs$coeftable[1,]
coefs_covs <- coef_table_covs[1]
std_error_covs <- coef_table_covs[2]
p_value_covs <- coef_table_covs[4]

# Additional model info
obs_covs <- summary_covs$nobs
adj_r2_covs <- r2(TWFE_covariates,'ar2')


### Summary covariates model
# Extract model summary details
summary_covs_2 <- summary(TWFE_covariates_2)

# Extract coefficients, standard errors, t-values, and p-values
coef_table_covs_2 <- summary_covs_2$coeftable[1,]
coefs_covs_2 <- coef_table_covs_2[1]
std_error_covs_2 <- coef_table_covs_2[2]
p_value_covs_2 <- coef_table_covs_2[4]

# Additional model info
obs_covs_2 <- summary_covs_2$nobs
adj_r2_covs_2 <- r2(TWFE_covariates_2,'ar2')


# Add footer and model summary statistics
latex_table <- paste(
  "\\begin{tabular}{lccc} \n",
  "\\toprule \n",
  "Variables & \\multicolumn{3}{c}{$ln(homicides \\, per \\, 100 k)$} \\\\ \n",
  "\\midrule \n",
  "$D_{it}$ &", round(coefs_basic,4), get_stars(p_value_basic), "&", round(coefs_covs,4), get_stars(p_value_covs), "&",
  round(coefs_covs_2, 4),get_stars(p_value_covs_2), "\\\\ \n",
  " & (" , round(std_error_basic,4), ") & (", round(std_error_covs,4), ") & (", round(std_error_covs_2,4) ,") \\\\ \n",
  "\\midrule \n",
  "Covariates & No & All & reduced \\\\ \n",
  "Observations & ", obs_basic, "&", obs_covs, "&", obs_covs_2, "\\\\ \n",
  "State FE & Yes & Yes & Yes \\\\ \n",
  "Year FE & Yes & Yes & Yes \\\\ \n",
  "Weighed by popwt & Yes & Yes & Yes \\\\ \n",
  "Adjusted $R\\textsuperscript{2}$ & ", round(adj_r2, 4), "&", round(adj_r2_covs, 4) , "&",round(adj_r2_covs_2, 4)  , "\\\\ \n",
  "\\bottomrule",
  "\\multicolumn{4}{c}{\\footnotesize Note: Standard errors in parentheses. *** p<0.01, ** p<0.05, * p<0.1.} \\\\ \n",
  "\\multicolumn{4}{c}{\\footnotesize In Covariates 'All' includes controls for robbery, assault, burglary, larceny, motor, murder, police,} \\\\ \n",
  "\\multicolumn{4}{c}{\\footnotesize unemployment, income, number of black and white males, prisoners, poverty rates and public wellfare.} \\\\ \n",
  "\\multicolumn{4}{c}{\\footnotesize 'reduced' includes only larceny and poverty rates.} \\\\ \n",
  "\\end{tabular}"
)

# Write the LaTeX table to a .tex file
cat(latex_table, file = paste0(outputs_path,"/PSET4_Q3.tex"))

# Q4 ----------------------------------------------------------------------

# gname: The name of the variable in data that contains the first period when a particular 
# observation is treated. This should be a positive number for all observations in treated groups. 
# It defines which "group" a unit belongs to. It should be 0 for units in the untreated group.

# for each state make a variable that says the first treated period
data[, first_treated_period := ifelse(any(treatment == 1), min(year[treatment == 1]), 0), by = state]

# estimate group-time average treatment effects using att_gt method
results_attgt <- att_gt(yname = "ln_homicide_per_100k",
                        tname = "year",
                        idname = "sid",
                        gname = "first_treated_period",
                        data = data
)

results_attgt_covs <- att_gt(yname = "ln_homicide_per_100k",
                        tname = "year",
                        idname = "sid",
                        gname = "first_treated_period",
                        data = data,
                        xformla = ~larceny + poverty 
)

### Group-Time Average Treatment Effects
agg.es <- aggte(results_attgt, type = "dynamic")

agg.es.covs <- aggte(results_attgt_covs, type = "dynamic")

# did plot
ggdid(agg.es)

# did plot
ggdid(agg.es.covs)

latex_table <- paste(
  "\\begin{tabular}{lcc} \n",
  "\\toprule \n",
  "Variables & \\multicolumn{2}{c}{$ln(homicides \\, per \\, 100 k)$} \\\\ \n",
  "\\midrule \n",
  "$D_{it}$ &", round(agg.es$overall.att,4), "&", round(agg.es.covs$overall.att,4), "\\\\ \n",
  " & (" , round(agg.es$overall.se,4), ") & (", round(agg.es.covs$overall.se,4), ") \\\\ \n",
  "\\midrule \n",
  "Covariates & No & reduced \\\\ \n",
  "Observations & ", obs_basic, "&", obs_covs, "\\\\ \n",
  "State FE & Yes & Yes \\\\ \n",
  "Year FE & Yes & Yes \\\\ \n",
  "Weighed by popwt & Yes & Yes \\\\ \n",
  "Adjusted $R\\textsuperscript{2}$ & ", round(adj_r2, 4), "&", round(adj_r2_covs, 4), "\\\\ \n",
  "\\bottomrule",
  "\\multicolumn{3}{c}{\\footnotesize Note: Standard errors in parentheses.} \\\\ \n",
  "\\multicolumn{3}{c}{\\footnotesize In Covariates 'All' includes controls for robbery, assault, burglary, larceny, motor, murder, police,} \\\\ \n",
  "\\multicolumn{3}{c}{\\footnotesize unemployment, income, number of black and white males, prisoners, poverty rates and public wellfare.} \\\\ \n",
  "\\multicolumn{3}{c}{\\footnotesize 'reduced' includes only larceny and poverty rates.} \\\\ \n",
  "\\end{tabular}"
)


# Write the LaTeX table to a .tex file
cat(latex_table, file = paste0(outputs_path,"/PSET4_Q4.tex"))
