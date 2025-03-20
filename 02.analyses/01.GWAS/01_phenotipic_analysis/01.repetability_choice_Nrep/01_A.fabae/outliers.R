# Load necessary libraries
library(emmeans)
library(lme4)
##loop ###########################################################################
# Define a vector of x values to loop through
x_values <- seq(1, 2, by = 0.05)

# Initialize an empty data frame to store results
results <- data.frame(x = numeric(0), AIC = numeric(0), BIC = numeric(0), RMSE = numeric(0), 
                      n_before = numeric(0), n_after = numeric(0), n_diff = numeric(0), stringsAsFactors = FALSE)

# Loop through different values of x
for (x in x_values) {
  # Calculate IQR and bounds for outlier removal
  IQR_value <- IQR(data250_fabae$larvae_total_number)
  lower_bound <- quantile(data250_fabae$larvae_total_number, 0.25) - x * IQR_value
  upper_bound <- quantile(data250_fabae$larvae_total_number, 0.75) + x * IQR_value
  
  # Number of rows before removing outliers
  n_before <- nrow(data250_fabae)
  
  # Filter out outliers
  data250_fabae_no_outliers <- data250_fabae[data250_fabae$larvae_total_number >= lower_bound & data250_fabae$larvae_total_number <= upper_bound, ]
  
  # Number of rows after removing outliers
  n_after <- nrow(data250_fabae_no_outliers)
  
  # Calculate the difference in number of observations
  n_diff <- n_before - n_after
  
  # Fit the model
  mod <- lmer(larvae_total_number_shifted ~ ivip + (1|session/block), data = data250_fabae_no_outliers)
  
  # Extract fitted values and residuals
  fitted_values <- fitted(mod)
  residuals_values <- residuals(mod)
  
  # Calculate AIC and BIC
  model_AIC <- AIC(mod)
  model_BIC <- BIC(mod)
  
  # Calculate RMSE (Root Mean Squared Error)
  RMSE <- sqrt(mean(residuals_values^2))
  
  # Store the results in the data frame
  results <- rbind(results, data.frame(x = x, AIC = model_AIC, BIC = model_BIC, RMSE = RMSE, 
                                       n_before = n_before, n_after = n_after, n_diff = n_diff))
  
  # Plot residual diagnostics (optional for each iteration)
  # Residuals vs Fitted
  plot(fitted_values, residuals_values, 
       main=paste("Residuals vs Fitted (x =", x, ")"), 
       xlab="Fitted values", 
       ylab="Residuals")
  abline(h=0, col="red")
  
  # Histogram of residuals
  hist(residuals_values, main=paste("Histogram of Residuals (x =", x, ")"), xlab="Residuals", breaks=20)
  
  # Q-Q plot for residuals
  qqnorm(residuals_values)
  qqline(residuals_values, col = "red")
  
  # Residuals over time
  plot(residuals(mod), type="o", col="blue", main=paste("Residuals over Time (x =", x, ")"), xlab="Index", ylab="Residuals")
}

# Print the results for comparison
print(results)


library(emmeans)
library(lme4)
library(ggplot2)
library(tidyr) 



# Reshape the results data frame to a long format for plotting
results_long <- results %>%
  pivot_longer(cols = c(AIC, BIC, RMSE), 
               names_to = "Metric", 
               values_to = "Value")

# Create the plots
ggplot(results_long, aes(x = x, y = Value, color = Metric)) +
  geom_line() + 
  geom_point() + 
  facet_wrap(~ Metric, scales = "free_y") +
  labs(title = "Model Comparison by Metric and Observations Removed", 
       x = "x (IQR Multiplier)",
       y = "Metric Value") +
  theme_minimal() +
  theme(legend.position = "bottom") +
  scale_color_manual(values = c("red", "blue", "green"))  # Customize colors for each metric

# Add a plot for n_diff (difference in observations removed)
ggplot(results, aes(x = x, y = n_diff)) + 
  geom_line(color = "purple") + 
  geom_point(color = "purple") +
  labs(title = "Difference in Observations Removed", 
       x = "x (IQR Multiplier)", 
       y = "Number of Observations Removed") +
  theme_minimal()


##what is filtered? ###########################################################################


x<-1.45
IQR_value <- IQR(data250_fabae$larvae_total_number)
lower_bound <- quantile(data250_fabae$larvae_total_number, 0.25) - x * IQR_value
upper_bound <- quantile(data250_fabae$larvae_total_number, 0.75) + x * IQR_value

# Filter out outliers
data250_fabae_no_outliers <- data250_fabae[data250_fabae$larvae_total_number >= lower_bound & data250_fabae$larvae_total_number <= upper_bound, ]

# Identify the rows that were removed (outliers)
removed_outliers <- data250_fabae[data250_fabae$larvae_total_number < lower_bound | data250_fabae$larvae_total_number > upper_bound, ]

# Display the DataFrame with the removed outliers
removed_outliers


# Check the Q-Q plot again without outliers
qqnorm(data250_fabae_no_outliers$larvae_total_number)
qqline(data250_fabae_no_outliers$larvae_total_number, col = "red")

# Number of rows before removing outliers
n_before <- nrow(data250_fabae)
print(n_before)

# Number of rows after removing outliers
n_after <- nrow(data250_fabae_no_outliers)
print(n_after)


mod<-lmer(larvae_total_number_shifted ~ ivip + (1|session/block),data = data250_fabae_no_outliers)


fitted_values<-fitted(mod)

# Extract fitted values and residuals
residuals_values <- residuals(mod)


#homosedasticity
# Plot residuals vs fitted values
plot(fitted_values, residuals_values, 
     main="Residuals vs Fitted", 
     xlab="Fitted values", 
     ylab="Residuals")
abline(h=0, col="red")

#residual normality
hist(residuals_values, main="Histogram of Residuals", xlab="Residuals", breaks=20)
# 2. Q-Q Plot for Residuals
qqnorm(residuals_values)
qqline(residuals_values, col = "red")

#residual independance
# Plot residuals over time (for time-series data)
plot(residuals(mod), type="o", col="blue", main="Residuals over Time", xlab="Index", ylab="Residuals")




emmeans_model <- emmeans(mod, ~ ivip)  # ~ ivip gives marginal means for the levels of ivip

coef_df <-as.data.frame(emmeans_model)


hist(coef_df$emmean)
# Define the file path
file_path <- "~/Documents/gitrepo/faba_specifics/02.analyses/01.GWAS/01_phenotipic_analysis/02.ajusted_means/BLUEs_fabae.csv"

# Save the data frame to the specified location
write.csv(coef_df, file = file_path, row.names = FALSE)


