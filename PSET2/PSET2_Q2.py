import pandas as pd
from rdrobust import rdrobust
from rddensity import rddensity

#############################################################
## Data Loading and Preparation
#############################################################

# Load dataset
headstart_data = pd.read_stata('headstart.dta')

# Define cutoff point
cutoff_value = 59.1984

# Extracting covariates from the 1960 Census
covariates_1960 = headstart_data[['census1960_pop',
                                    'census1960_pctsch1417',
                                    'census1960_pctsch534', 
                                    'census1960_pctsch25plus',
                                    'census1960_pop1417',
                                    'census1960_pop534',
                                    'census1960_pop25plus',
                                    'census1960_pcturban',
                                    'census1960_pctblack']]

# Extracting covariates from the 1990 Census
covariates_1990 = headstart_data[['census1990_pop',
                                    'census1990_pop1824',
                                    'census1990_pop2534',
                                    'census1990_pop3554',
                                    'census1990_pop55plus',
                                    'census1990_pcturban',
                                    'census1990_pctblack',
                                    'census1990_percapinc']]

# Define the outcome variable, running variable, and treatment
outcome_variable = headstart_data['mort_age59_related_postHS']
raw_running_variable = headstart_data['povrate60']
running_variable = raw_running_variable - cutoff_value

treatment_indicator = (running_variable >= 0).astype(int)

# Placebo outcomes
placebo_outcomes = headstart_data[['mort_age59_injury_postHS', 
                                     'mort_age59_related_preHS']]

#############################################################
## Figure 1: Scatter Plot and RD Analysis
#############################################################

# Filter for valid data points
valid_indices = (outcome_variable < 20) & (outcome_variable.notna()) & (raw_running_variable.notna())

# Uncomment the following lines to visualize the RD plot
# rdplot(outcome_variable[valid_indices], raw_running_variable[valid_indices], c=cutoff_value, nbins=3000)

###################################################################
## Table 1: Binomial Tests
###################################################################

# Initialize results storage for binomial tests
binomial_results = []

for threshold in range(3, 14):  # Adjust thresholds as needed
    selected_indices = (abs(running_variable) <= threshold) & outcome_variable.notna() & raw_running_variable.notna()
    
    if selected_indices.sum() > 0:  # Ensure there are valid observations
        count_not_treated = sum(1 - treatment_indicator[selected_indices])
        count_treated = sum(treatment_indicator[selected_indices])
        
        binomial_test_result = sum(treatment_indicator[selected_indices])
        p_value = binom_test(binomial_test_result, count_treated + count_not_treated, p=0.5)
        
        binomial_results.append([threshold, count_not_treated, count_treated, p_value])

# Convert results to DataFrame for easy viewing
binomial_df = pd.DataFrame(binomial_results, columns=['Threshold', 'Not Treated', 'Treated', 'P-value'])
print(binomial_df.round(3))

###################################################################
## Table 2: Nonparametric Density Continuity Tests
###################################################################

# Initialize results storage for density tests
density_results = []

for bw_selection in [None, 'diff', 'restricted']:
    density_test_result = rddensity(raw_running_variable.dropna(), c=cutoff_value, bwselect=bw_selection)
    density_results.append([density_test_result.h[0], density_test_result.h[1],
                            density_test_result.N[4], density_test_result.N[5], 
                            density_test_result.test[4]])

# Convert results to DataFrame for easy viewing
density_df = pd.DataFrame(density_results, columns=['Bandwidth 1', 'Bandwidth 2', 'N 1', 'N 2', 'Test'])
print(density_df.round(3))

#############################################################
## Table 3: Flexible Parametric RD Methods
#############################################################

# Initialize storage for parametric RD results
parametric_results = []

for order in [1, 4]:  # 1 for linear, 4 for quartic
    for h in [9, 18, 20, 100]:  # Different bandwidths
        selected_indices = abs(running_variable) <= h
        if selected_indices.sum() > 0:
            polynomial_terms = pd.DataFrame({'Poly': np.poly(running_variable[selected_indices], order)})
            
            # Fit linear model
            model = sm.OLS(outcome_variable[selected_indices], 
                           sm.add_constant(pd.concat([treatment_indicator[selected_indices], polynomial_terms], axis=1))).fit()
            parametric_results.append([order, h, model.params[1], model.conf_int().iloc[1]])
    
# Convert results to DataFrame
parametric_df = pd.DataFrame(parametric_results, columns=['Order', 'Bandwidth', 'Coefficient', 'Confidence Interval'])
print(parametric_df.round(3))

#############################################################
## Table 4: Robust Nonparametric Local Polynomial Methods
#############################################################

# Initialize storage for local polynomial results
local_poly_results = []

for p in [0, 1]:  # Local constant and linear
    for h in [None, 9]:  # Default bandwidth and custom bandwidth
        result = rdrobust(outcome_variable, running_variable, p=p, h=h)
        local_poly_results.append([result.p, result.bws[0], result.coef[0], result.ci[0], result.ci[1], result.pv[0], result.N_h[0], result.N_h[1]])

# Convert results to DataFrame
local_poly_df = pd.DataFrame(local_poly_results, columns=['P-value', 'Bandwidth', 'Coefficient', 'CI Lower', 'CI Upper', 'PV', 'N Left', 'N Right'])
print(local_poly_df.round(3))

#############################################################
# Figure 2: Window Selection and Outcome Visualization
#############################################################

# Prepare data for window selection
X_window = pd.concat([headstart_data['mort_age59_related_preHS'], covariates_1960], axis=1)

# Select optimal window
window_selection = rdwinselect(running_variable, X_window, reps=1000, statistic="ksmirnov", wmin=0.3, wstep=0.2, level=0.2)

# Plot P-values
rdwinselect(running_variable, X_window, reps=1000, statistic="ksmirnov", wmin=0.3, wstep=0.2, level=0.2, nwindows=40, plot=True, quietly=True)

# Uncomment the following lines for scatter plot
valid_window_indices = abs(running_variable) <= 1.1
plt.scatter(running_variable[valid_window_indices], outcome_variable[valid_window_indices])
plt.title('Outcome by Running Variable')
plt.xlabel('Running Variable')
plt.ylabel('Outcome Variable')
plt.show()
