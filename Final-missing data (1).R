library(tidyverse)
library(readxl)
library(imputeTS)
library(VIM)
library(dplyr)
library(forecast)
library(ggplot2)

water_data = read_excel("Water-Data.xlsx")

# Relevant columns
water_data = water_data[c("Date", "pH", "TURB", "temperature", "ChemCostper1000")]
head(water_data)

# Simulate MCAR missingness
simulate_MCAR = function(data, missing_rate) {
  data_mcar = data
  for (var in colnames(data)[-1]) {  # Exclude 'Date' column
    missing_indices <- sample(1:nrow(data), size = floor(missing_rate * nrow(data)))
    data_mcar[missing_indices, var] <- NA
  }
  return(data_mcar)
}

# Simulate MAR missingness
simulate_MAR = function(data, missing_rate) {
  data_mar = data
  for (var in colnames(data)[-1]) {  # Exclude 'Date' column
    if (var == "pH") {
      missing_indices = which(runif(nrow(data)) < missing_rate & !is.na(data$TURB))
      data_mar[missing_indices, var] <- NA
    } else {
      missing_indices = sample(1:nrow(data), size = floor(missing_rate * nrow(data)))
      data_mar[missing_indices, var] = NA
    }
  }
  return(data_mar)
}

# Simulate MNAR missingness
simulate_MNAR = function(data, missing_rate) {
  data_mnar <- data
  for (var in colnames(data)[-1]) {  # Exclude 'Date' column
    if (var == "pH") {
      missing_indices = which(runif(nrow(data)) < missing_rate * (data$pH / max(data$pH, na.rm = TRUE)))
      data_mnar[missing_indices, var] = NA
    } else {
      missing_indices = sample(1:nrow(data), size = floor(missing_rate * nrow(data)))
      data_mnar[missing_indices, var] = NA
    }
  }
  return(data_mnar)
}

# 1000 simulations 
set.seed(123)  
num_simulations = 50
missing_rate = 0.1 

# List to store simulated datasets
mcar_simulations = vector("list", num_simulations)
mar_simulations = vector("list", num_simulations)
mnar_simulations = vector("list", num_simulations)

for (i in 1:num_simulations) {
  mcar_simulations[[i]] = simulate_MCAR(water_data, missing_rate)
  mar_simulations[[i]] = simulate_MAR(water_data, missing_rate)
  mnar_simulations[[i]] = simulate_MNAR(water_data, missing_rate)
}

head(mcar_simulations[[1]])
head(mar_simulations[[1]])
head(mnar_simulations[[1]])


# Visualize missing data patterns 
visualize_missing_patterns = function(data_list, mechanism) {
  combined_data = do.call(rbind, lapply(data_list, function(data) {
    data %>%
      mutate(Mechanism = mechanism)
  }))
  return(combined_data)
}

mcar_combined = visualize_missing_patterns(mcar_simulations, "MCAR")
mar_combined = visualize_missing_patterns(mar_simulations, "MAR")
mnar_combined = visualize_missing_patterns(mnar_simulations, "MNAR")

all_combined = bind_rows(mcar_combined, mar_combined, mnar_combined)
aggr(all_combined[, -1], col = c('navyblue', 'red'), numbers = TRUE, sortVars = TRUE, labels = names(all_combined)[-1], cex.axis = 0.7, gap = 3, ylab = c("Missing data", "Pattern"))

# Kalman smoothing for imputation
apply_kalman = function(data) {
  data_imputed = data
  for (col in colnames(data)[-1]) {  # Exclude 'Date' column
    data_imputed[[col]] = na_kalman(data[[col]], model = "auto.arima")
  }
  return(data_imputed)
}

mcar_imputed = lapply(mcar_simulations, apply_kalman)
mar_imputed = lapply(mar_simulations, apply_kalman)
mnar_imputed = lapply(mnar_simulations, apply_kalman)

head(mcar_imputed[[1]])
head(mar_imputed[[1]])
head(mnar_imputed[[1]])

# Function to calculate RMSE between original and imputed values
calculate_rmse = function(original, imputed, var) {
  original_values = original[[var]]
  imputed_values = imputed[[var]]
  missing_indices = which(is.na(original_values))
  valid_indices = which(!is.na(original_values) & !is.na(imputed_values))
  rmse = sqrt(mean((original_values[valid_indices] - imputed_values[valid_indices])^2, na.rm = TRUE))
  return(rmse)
}

# Summarize RMSE across simulations and variables
summarize_rmse = function(original_data, imputed_data_list, mechanism) {
  results = data.frame()
  for (i in seq_along(imputed_data_list)) {
    imputed_data = imputed_data_list[[i]]
    for (var in colnames(original_data)[-1]) {
      rmse = calculate_rmse(original_data, imputed_data, var)
      results = rbind(results, data.frame(
        simulation = i,
        variable = var,
        mechanism = mechanism,
        rmse = rmse
      ))
    }
  }
  return(results)
}

# Calculate RMSE summaries for all mechanisms
mcar_rmse = summarize_rmse(water_data, mcar_imputed, "MCAR")
mar_rmse = summarize_rmse(water_data, mar_imputed, "MAR")
mnar_rmse = summarize_rmse(water_data, mnar_imputed, "MNAR")

rmse_summary = bind_rows(mcar_rmse, mar_rmse, mnar_rmse)

# Calculate overall RMSE averages
overall_rmse = rmse_summary %>%
  group_by(mechanism) %>%
  summarize(mean_rmse = mean(rmse), .groups = "drop")

overall_rmse

# Combine RMSE results
rmse_summary = bind_rows(mcar_rmse, mar_rmse, mnar_rmse)

# Visualize RMSE 
ggplot(rmse_summary, aes(x = variable, y = rmse, fill = mechanism)) +
  geom_boxplot() +
  labs(
    title = "RMSE of Kalman Imputation by Missingness Mechanism",
    x = "Variable",
    y = "RMSE",
    fill = "Mechanism"
  ) +
  theme_minimal()

# Function to calculate MAE between original and imputed values
calculate_mae = function(original, imputed, var) {
  original_values = original[[var]]
  imputed_values = imputed[[var]]
  missing_indices = which(is.na(original_values))
  valid_indices = which(!is.na(original_values) & !is.na(imputed_values))
  mae = mean(abs(original_values[valid_indices] - imputed_values[valid_indices]), na.rm = TRUE)
  return(mae)
}

# Summarize MAE across simulations and variables
summarize_mae = function(original_data, imputed_data_list, mechanism) {
  results = data.frame()
  for (i in seq_along(imputed_data_list)) {
    imputed_data = imputed_data_list[[i]]
    for (var in colnames(original_data)[-1]) {
      mae = calculate_mae(original_data, imputed_data, var)
      results = rbind(results, data.frame(
        simulation = i,
        variable = var,
        mechanism = mechanism,
        mae = mae
      ))
    }
  }
  return(results)
}

# Calculate MAE summaries for all mechanisms
mcar_mae = summarize_mae(water_data, mcar_imputed, "MCAR")
mar_mae = summarize_mae(water_data, mar_imputed, "MAR")
mnar_mae = summarize_mae(water_data, mnar_imputed, "MNAR")

# Combine MAE results
mae_summary = bind_rows(mcar_mae, mar_mae, mnar_mae)

# Visualize MAE across mechanisms and variables
ggplot(mae_summary, aes(x = variable, y = mae, fill = mechanism)) +
  geom_boxplot() +
  labs(
    title = "MAE of Kalman Imputation by Missingness Mechanism",
    x = "Variable",
    y = "MAE",
    fill = "Mechanism"
  ) +
  theme_minimal()

# Compare MAE values
overall_mae =mae_summary %>%
  group_by(mechanism) %>%
  summarize(mean_mae = mean(mae, na.rm = TRUE))
overall_mae

# Function to calculate bias, variance, and confidence intervals
calculate_bias_variance_ci = function(original, imputed, var) {
  original_values = original[[var]]
  imputed_values = imputed[[var]]
  missing_indices = which(is.na(original_values))
  valid_indices = which(!is.na(original_values) & !is.na(imputed_values))
  
  bias = mean(imputed_values[valid_indices] - original_values[valid_indices], na.rm = TRUE)
  variance = var(imputed_values[valid_indices], na.rm = TRUE)
  ci = t.test(imputed_values[valid_indices])$conf.int
  
  return(list(bias = bias, variance = variance, ci = ci))
}

# Summarize bias, variance, and confidence intervals across simulations and variables
summarize_bias_variance_ci = function(original_data, imputed_data_list, mechanism) {
  results = data.frame()
  for (i in seq_along(imputed_data_list)) {
    imputed_data = imputed_data_list[[i]]
    for (var in colnames(original_data)[-1]) {
      bvci = calculate_bias_variance_ci(original_data, imputed_data, var)
      results = rbind(results, data.frame(
        simulation = i,
        variable = var,
        mechanism = mechanism,
        bias = bvci$bias,
        variance = bvci$variance,
        ci_lower = bvci$ci[1],
        ci_upper = bvci$ci[2]
      ))
    }
  }
  return(results)
}

# Calculate bias, variance, and confidence intervals for all mechanisms
mcar_bvci = summarize_bias_variance_ci(water_data, mcar_imputed, "MCAR")
mar_bvci = summarize_bias_variance_ci(water_data, mar_imputed, "MAR")
mnar_bvci = summarize_bias_variance_ci(water_data, mnar_imputed, "MNAR")

# Combine results
bvci_summary = bind_rows(mcar_bvci, mar_bvci, mnar_bvci)

# Visualize bias, variance, and confidence intervals
ggplot(bvci_summary, aes(x = variable, y = bias, fill = mechanism)) +
  geom_boxplot() +
  labs(
    title = "Bias of Kalman Imputation by Missingness Mechanism",
    x = "Variable",
    y = "Bias",
    fill = "Mechanism"
  ) +
  theme_minimal()

ggplot(bvci_summary, aes(x = variable, y = variance, fill = mechanism)) +
  geom_boxplot() +
  labs(
    title = "Variance of Kalman Imputation by Missingness Mechanism",
    x = "Variable",
    y = "Variance",
    fill = "Mechanism"
  ) +
  theme_minimal()

ggplot(bvci_summary, aes(x = variable, y = ci_upper - ci_lower, fill = mechanism)) +
  geom_boxplot() +
  labs(
    title = "Confidence Interval Width of Kalman  by Missingness",
    x = "Variable",
    y = "CI Width",
    fill = "Mechanism"
  ) +
  theme_minimal()

# Compare bias, variance, and confidence intervals
summary_stats = bvci_summary %>%
  group_by(mechanism) %>%
  summarize(
    mean_bias = mean(bias, na.rm = TRUE),
    mean_variance = mean(variance, na.rm = TRUE),
    mean_ci_width = mean(ci_upper - ci_lower, na.rm = TRUE)
  )



# table
mechanism_summary = tibble::tibble(
  mechanism = c("MAR", "MCAR", "MNAR"),
  mean_bias = c(0.000504, -0.00162, -0.000793),
  mean_variance = c(29.5, 29.4, 29.4),
  mean_ci_width = c(0.345, 0.345, 0.345)
)

mechanism_summary

#Comparison of Kalman and Exponential
# Function to apply exponential smoothing for imputation
apply_exponential_smoothing = function(data) {
  data_imputed <- data
  for (col in colnames(data)[-1]) {  # Exclude 'Date' column if applicable
    ts_data = ts(data[[col]], frequency = 1)
    if (any(is.na(ts_data))) {
      # Interpolate missing values first
      ts_data = na.interp(ts_data)
    }
    ets_model = ets(ts_data)
    # Use the fitted values for imputation
    data_imputed[[col]] <- fitted(ets_model)
  }
  return(data_imputed)
}

# Impute missing data using exponential smoothing
mar_imputed_ets = lapply(mar_simulations, apply_exponential_smoothing)

# RSME for expo
mar_rmse_exp = summarize_rmse(water_data, mar_imputed_ets, "MAR")

# MAE for expo
mar_mae_exp = summarize_mae(water_data, mar_imputed_ets, "MAR")

#BVCI
mar_bvci_exp = summarize_bias_variance_ci(water_data, mar_imputed_ets, "MAR")

# Mean of RMSE

str(mar_rmse_exp)
str(mar_rmse)

# Extract the RMSE column and compute the mean
mean_rmse_exp = mean(mar_rmse_exp$rmse, na.rm = TRUE)
mean_rmse_kalman = mean(mar_rmse$rmse, na.rm = TRUE)

mean_rmse_exp
mean_rmse_kalman

# Mean of MAE
mean_mae_exp = mean(mar_mae_exp$mae, na.rm = TRUE)
mean_mae_kalman = mean(mar_mae$mae, na.rm = TRUE)

mean_mae_exp
mean_mae_kalman

# Mean of BVCI
mean_bvci_exp = mean(mar_bvci_exp$bias, na.rm = TRUE)
mean_bvci_kalman <- mean(mar_bvci$bias, na.rm = TRUE)

mean_bvci_exp
mean_bvci_kalman

# Mean of BVCI1
mean_bvci1_exp <- mean(mar_bvci_exp$variance, na.rm = TRUE)
mean_bvci1_kalman <- mean(mar_bvci$variance, na.rm = TRUE)


# Create a comparison table
comparison_results <- data.frame(
  Metric = c("RMSE", "MAE", "Bias", "Variance"),
  Exponential_Smoothing = c(mean_rmse_exp, mean_mae_exp, mean_bvci_exp,mean_bvci1_exp ),
  Kalman_SMoothing = c(mean_rmse_kalman, mean_mae_kalman, mean_bvci_kalman,mean_bvci1_kalman )
)

comparison_results

# Add a method column to distinguish between the two datasets
mar_rmse_exp$method <- "Exponential Smoothing"
mar_rmse$method <- "Kalman_smoothing"

# Combine the two data frames
rmse_combined <- rbind(mar_rmse_exp, mar_rmse)

# Plot the RMSE values
ggplot(rmse_combined, aes(x = variable, y = rmse, fill = method)) +
  geom_boxplot() +
  labs(
    title = "Comparison of RMSE by Variable and Method",
    x = "Variable",
    y = "RMSE",
    fill = "Imputation Method"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# Add a method column to distinguish between the two datasets
mar_mae_exp$method <- "Exponential Smoothing"
mar_mae$method <- "Kalman_smoothing"

# Combine the two data frames
mae_combined <- rbind(mar_mae_exp, mar_mae)

# Plot the MAE values
ggplot(mae_combined, aes(x = variable, y = mae, fill = method)) +
  geom_boxplot() +
  labs(
    title = "Comparison of MAE by Variable and Method",
    x = "Variable",
    y = "MAE",
    fill = "Imputation Method"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# Add a method column to distinguish between the two datasets
mar_bvci_exp$method <- "Exponential Smoothing"
mar_bvci$method <- "Kalman Smoothing"

# Combine the two data frames
bvci_combined <- rbind(mar_bvci_exp, mar_bvci)

#plots
# Reshape the data to long format
bvci_long <- bvci_combined %>%
  pivot_longer(cols = c(bias, variance, ci_lower, ci_upper),  # Replace with your actual column names for variance and CI
               names_to = "metric",
               values_to = "value")


# Plot BIAS values
ggplot(bvci_long %>% filter(metric == "bias"), aes(x = variable, y = value, fill = method)) +
  geom_boxplot() +
  labs(
    title = "Comparison of BIAS by Variable and Method",
    x = "Variable",
    y = "BIAS",
    fill = "Imputation Method"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Plot Variance values
ggplot(bvci_long %>% filter(metric == "variance"), aes(x = variable, y = value, fill = method)) +
  geom_boxplot() +
  labs(
    title = "Comparison of Variance by Variable and Method",
    x = "Variable",
    y = "Variance",
    fill = "Imputation Method"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Reshape the data to long format for ci_lower and ci_upper
bvci_long_ci <- bvci_combined %>%
  pivot_longer(cols = c(ci_lower, ci_upper),  # Pivot for ci_lower and ci_upper
               names_to = "CI_bound",
               values_to = "value")

# Plot both ci_lower and ci_upper
ggplot(bvci_long_ci, aes(x = variable, y = value, fill = CI_bound)) +
  geom_boxplot() +
  labs(
    title = "Comparison of Confidence Interval (Lower and Upper Bounds)",
    x = "Variable",
    y = "Confidence Interval Value",
    fill = "CI Bound"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
