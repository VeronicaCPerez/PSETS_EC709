library(dplyr)
library(tidyr)
library(kableExtra)
library(data.table)
library(xtable)
library(ggplot2)
library(fixest)
library(did)
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
                           whitem_25_44 + prisoner + lagprisoner + poverty +
                           exp_subsidy + exp_pubwelfare 
                           | state + year,data,weights = ~ popwt)

# Q4 ----------------------------------------------------------------------

# gname: The name of the variable in data that contains the first period when a particular 
# observation is treated. This should be a positive number for all observations in treated groups. 
# It defines which "group" a unit belongs to. It should be 0 for units in the untreated group.

# for each state make a variable that says the first treated period
data[, first_treated_period := ifelse(any(treatment == 1), min(year[treatment == 1]), 0), by = state]

# estimate group-time average treatment effects using att_gt method
example_attgt <- att_gt(yname = "ln_homicide_per_100k",
                        tname = "year",
                        idname = "sid",
                        gname = "first_treated_period",
                        data = data
)

# summarize the results
summary(example_attgt)


### Aggregate results
agg.es <- aggte(example_attgt, type = "dynamic")
summary(agg.es)

# did plot
ggdid(agg.es)



