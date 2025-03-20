

data250_afabae_00 <- read.csv("01.data/2023_screen_specifics_251genotypes/02.aphis_fabae_on_faba_collection/aphis_fabae_fecundity_on_faba_collection_2024-04-02.csv", stringsAsFactors = FALSE, sep = ",") %>%  
  mutate(ivip = as.factor(ivip)) %>% 
  mutate(session = as.factor(session)) 



data250_afabae_01 <- data250_afabae_00 %>% filter(!is.na(larvae_total_number))
data250_afabae_01$ivip <- droplevels(data250_afabae_01$ivip)

data250_afabae_02<- data250_afabae_01 %>% filter(adult_alive != "0") #remove if no adult
data250_afabae_02$ivip <- droplevels(data250_afabae_02$ivip)

data250_afabae_02$ivip <-as.factor(data250_afabae_02$ivip)


x_values <- seq(1, 2, by = 0.05)

# Initialize an empty data frame to store results
results <- data.frame(x = numeric(0), i = numeric(0), AIC = numeric(0), BIC = numeric(0), RMSE = numeric(0), 
                      n_before = numeric(0), n_after = numeric(0), n_diff = numeric(0), H2=numeric(0),shapiro=numeric(0),NGEN=numeric(0), stringsAsFactors = FALSE)

# Loop through different values of x
for (x in x_values) {
  print(x)
  for (i in 3:7) {
    print(i)
    # Calculate IQR and bounds for outlier removal
    IQR_value <- IQR(data250_afabae_02$larvae_total_number)
    lower_bound <- quantile(data250_afabae_02$larvae_total_number, 0.25) - x * IQR_value
    upper_bound <- quantile(data250_afabae_02$larvae_total_number, 0.75) + x * IQR_value
    
    # Filter out outliers
    data <- data250_afabae_02[data250_afabae_02$larvae_total_number >= lower_bound & data250_afabae_02$larvae_total_number <= upper_bound, ]
    
    data <- data %>%
      group_by(ivip) %>%           # Group by genotype
      filter(n() >= i) %>%             # Keep only groups with i or more repetitions
      ungroup() 
    
    
    data$ivip <- droplevels(data$ivip)
    
    NGEN<-nlevels(data$ivip)
    
   
    # Number of rows before removing outliers
    n_before <- nrow(data250_afabae_02)
    
    
    # Number of rows after removing outliers
    n_after <- nrow(data)
    
    # Calculate the difference in number of observations
    n_diff <- n_before - n_after
    
    # Fit the model
    mod <- lmer(larvae_total_number ~ ivip + (1|session/block), data = data)
    
    # Extract fitted values and residuals
    fitted_values <- fitted(mod)
    residuals_values <- residuals(mod)
    
    # Perform the Shapiro-Wilk test
    shapiro_result <- shapiro.test(residuals_values)
    
    # Print the p-value (this is your quantification of normality)
    shapiro <- shapiro_result$p.value
    
    # Calculate AIC and BIC
    model_AIC <- AIC(mod)
    model_BIC <- BIC(mod)
    
    # Calculate RMSE (Root Mean Squared Error)
    RMSE <- sqrt(mean(residuals_values^2))
    
    # Fit the model with a fixed formula
    mod_h2 <- lmer(larvae_total_number ~ (1 | ivip) + (1 | session), data = data)
    
    # Get the variance components
    estimates <- as.data.frame(VarCorr(mod_h2))
    
    # Extract the relevant variance components
    genetic_variance <- estimates$vcov[1]
    residual_variance <- estimates$vcov[3]  # Adjust based on the placement of residual variance
    
    nbE <-nlevels(data$session)
    
    # Calculate heritability (h2)
    H2 <- genetic_variance / (genetic_variance + (residual_variance / nbE))
    
    
    
    # Store the results in the data frame
    results <- rbind(results, data.frame(x = x, i = i, AIC = model_AIC, BIC = model_BIC, RMSE = RMSE, 
                                         n_before = n_before, n_after = n_after, n_diff = n_diff, H2 = H2, shapiro = shapiro, NGEN=NGEN))
  }
    
    
}

library(ggplot2)

ggplot(results, aes(x = n_after, y = shapiro, color = i)) + 
  geom_point() +  # Scatter plot with colored points based on 'x'
  labs(title = "Shapiro-Wilk Test p-values vs n_after", 
       x = "n_after", 
       y = "Shapiro-Wilk p-value", 
       color = "X Variable") +  # Color legend label
  theme_minimal()


ggplot(results, aes(x = NGEN, y = AIC, color = x)) + 
  geom_point() +  # Scatter plot with colored points based on 'x'
  labs(title = "Shapiro-Wilk Test p-values vs n_after", 
       x = "NGEN", 
       y = "AIC", 
       color = "X Variable") +  # Color legend label
  theme_minimal()

ggplot(results, aes(x = NGEN, y = RMSE, color = i)) + 
  geom_point() +  # Scatter plot with colored points based on 'x'
  labs(title = "Shapiro-Wilk Test p-values vs n_after", 
       x = "n_after", 
       y = "Shapiro-Wilk p-value", 
       color = "X Variable") +  # Color legend label
  theme_minimal()


ggplot(results, aes(x = NGEN, y = H2, color = x)) + 
  geom_point() +  # Scatter plot with colored points based on 'x'
  labs(title = "Shapiro-Wilk Test p-values vs n_after", 
       x = "NGEN", 
       y = "H2", 
       color = "X Variable") +  # Color legend label
  theme_minimal()


ggplot(results, aes(x = AIC, y = H2, color = i)) + 
  geom_point() +  # Scatter plot with colored points based on 'x'
  labs(title = "Shapiro-Wilk Test p-values vs n_after", 
       x = "n_after", 
       y = "Shapiro-Wilk p-value", 
       color = "X Variable") +  # Color legend label
  theme_minimal()

ggplot(results, aes(x = RMSE, y = AIC, color = i)) + 
  geom_point() +  # Scatter plot with colored points based on 'x'
  labs(title = "Shapiro-Wilk Test p-values vs n_after", 
       x = "n_after", 
       y = "Shapiro-Wilk p-value", 
       color = "X Variable") +  # Color legend label
  theme_minimal()
