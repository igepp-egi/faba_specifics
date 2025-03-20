library(gamlss)

data250_pisum<-read.csv("02.analyses/01.GWAS/01_phenotipic_analysis/acyrthosiphon_pisum_filtered.csv")
data250_pisum<-na.omit(data250_pisum)

data250_fabae<-read.csv("02.analyses/01.GWAS/01_phenotipic_analysis/aphis_fabae_data_filtered.csv")
data250_fabae<-na.omit(data250_fabae)

summary(data250_fabae)


# List of all distribution families to compare
families <- c("NO","BE", "BEOI", "BEZI", "BEINF", "BCCG", "BCPE", "BCT", "EXP", "exGAUS", "EGB2", 
              "GA", "GB1", "GB2", "GG", "GIG", "GT", "GU", "IG", "JSU", "JSUo", "LO", "LOGNO", 
              "LNO", "NET", "NO", "NOF", "PE", "RG", "SEP1", "SEP2", "SEP3", "SEP4", "SHASH", 
              "ST1", "ST2", "ST3", "ST4", "ST5", "TF", "WEI", "WEI2", "WEI3", "ZAIG", "BB", "BI", 
              "DEL", "NBI", "NBII", "PO", "PIG", "SI", "SICHEL", "ZIP", "ZIP2")

# Initialize a data frame to store AIC, BIC, and Family
results <- data.frame(Family = character(0), AIC = numeric(0), BIC = numeric(0))

# Loop over each family and fit the model
for(fam in families) {
  # Fit the model
  model <- tryCatch({
    gamlss(larvae_total_number ~ ivip + session/block, data = data250_pisum, family = get(fam)())
  }, error = function(e) {
    return(NULL)  # If the model doesn't converge or fails, return NULL
  })
  
  # If model is successfully fitted, extract AIC and BIC
  if (!is.null(model)) {
    aic_val <- AIC(model)
    bic_val <- BIC(model)
    
    # Store the results
    results <- rbind(results, data.frame(Family = fam, AIC = aic_val, BIC = bic_val))
  }
}

# Sort the results by AIC or BIC to find the best model
results_sorted_by_AIC <- results[order(results$AIC), ]
results_sorted_by_BIC <- results[order(results$BIC), ]

# View the sorted results
print("Sorted by AIC:")
print(results_sorted_by_AIC)

print("Sorted by BIC:")
print(results_sorted_by_BIC)

modbis<-mod


data250_pisum$ivip<-as.factor(data250_pisum$ivip)
data250_pisum$block<-as.factor(data250_pisum$block)
data250_pisum$session<-as.factor(data250_pisum$session)

hist(data250_fabae$larvae_total_number)
hist(data250_pisum$larvae_total_number)

#mod<-gamlss(larvae_total_number ~ ivip + session/block, data = data250_pisum, family = SEP2)
#mod<-gamlss(larvae_total_number ~ ivip + session/block, data = data250_pisum, family = SE)
for (i in 1) {
  
mod0<-lm(larvae_total_number ~ ivip + session/block, data = data250_pisum)
print("0")
#mod1<-gamlss(larvae_total_number ~ ivip + session, data = data250_pisum, family = ST3)
print("1")
mod2<-gamlss(larvae_total_number ~ ivip + random(session), data = data250_fabae, family = ST3)
print("2")
#mod3<-gamlss(larvae_total_number ~ ivip + session/block, data = data250_pisum, family = ST3)
print("3")
mod4<-gamlss(larvae_total_number ~ ivip, random = ~ session/block, data = data250_pisum, family = ST3)
print("4")
}

summary(data250_pisum)
summary(data250_pisum$ivip)
table(data250_pisum$ivip)  # Check for lack of variability

compare_performance(mod0, mod2, mod4, rank = FALSE, verbose = TRUE)

# Load necessary package
library(lme4)

data250_pisum$larvae_total_number_shifted <- data250_pisum$larvae_total_number + 1

# Fit Gamma GLMM with identity link
model_identity <- glmer(larvae_total_number_shifted ~ ivip + (1|session), 
                        family = Gamma(link = "identity"), 
                        data = data250_pisum)


# Fit Gamma GLMM with log link
model_log <- glmer(larvae_total_number_shifted ~ ivip + (1|session), 
                   family = Gamma(link = "log"), 
                   data = data250_pisum)

# Fit Gamma GLMM with inverse link
model_inverse <- glmer(larvae_total_number_shifted ~ ivip + (1|session), 
                       family = Gamma(link = "inverse"), 
                       data = data250_pisum)



# Check summaries of each model
summary(model_identity)
summary(model_log)
summary(model_inverse)




# Compare AIC and BIC for model selection
AIC(model_identity, model_log, model_inverse)
BIC(model_identity, model_log, model_inverse)

# Q-Q plot for normality check
qqnorm(data250_fabae$larvae_total_number)
qqline(data250_fabae$larvae_total_number, col = "red")  # Add reference line



# Fit Poisson GLMM
model_poisson <- glmer(larvae_total_number ~ ivip + (1|session), 
                       family = poisson(link = "log"), 
                       data = data250_fabae)



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












coef_df <- as.data.frame(coef(mod4))
# Convert row names to a new column
coef_df <- cbind(ivip = rownames(coef_df), coef_df)

# Now remove the original row names (optional)
rownames(coef_df) <- NULL

# 1. Remove "ivip" prefix from the 'ivip' column
coef_df$ivip <- sub("^ivip", "", coef_df$ivip)

# 2. Rename the columns
colnames(coef_df) <- c("ivip", "BLUE")

# 3. Remove the last row from the dataframe
coef_df <- coef_df[-nrow(coef_df), ]

# 4. Add the value of the first 'BLUE' row to all other rows
first_blue_value <- coef_df$BLUE[1]
coef_df$BLUE <- coef_df$BLUE + first_blue_value

# 5. Remove the first row after the addition
coef_df <- coef_df[-1, ]

# View the modified data
print(coef_df)

hist(coef_df$BLUE)

# Define the file path
file_path <- "~/Documents/gitrepo/faba_specifics/02.analyses/01.GWAS/01_phenotipic_analysis/02.ajusted_means/BLUEs_pisum.csv"

# Save the data frame to the specified location
write.csv(coef_df, file = file_path, row.names = FALSE)







# Install and load 'e1071' package to calculate skewness
install.packages("e1071")
library(e1071)

skewness(data250_pisum$larvae_total_number)
skewness(data250_fabae$larvae_total_number)

# Visualize with a histogram or density plot
hist(data250_pisum$larvae_total_number, main = "Histogram of larvae_total_number (Pisum)", xlab = "larvae_total_number")


# Apply log transformation to reduce left skewness
data250_fabae$larvae_log <- log(data250_fabae$larvae_total_number + 1)  # Adding 1 to avoid log(0)
qqnorm(data250_fabae$larvae_log)
qqline(data250_fabae$larvae_log, col = "red")

# Apply square root transformation
data250_fabae$larvae_sqrt <- sqrt(data250_fabae$larvae_total_number)
qqnorm(data250_fabae$larvae_sqrt)
qqline(data250_fabae$larvae_sqrt, col = "red")



# Calculate IQR
x<-1.75
IQR_value <- IQR(data250_fabae$larvae_total_number)
lower_bound <- quantile(data250_fabae$larvae_total_number, 0.25) - x * IQR_value
upper_bound <- quantile(data250_fabae$larvae_total_number, 0.75) + x * IQR_value

# Filter out outliers
data250_fabae_no_outliers <- data250_fabae[data250_fabae$larvae_total_number >= lower_bound & data250_fabae$larvae_total_number <= upper_bound, ]

# Check the Q-Q plot again without outliers
qqnorm(data250_fabae_no_outliers$larvae_total_number)
qqline(data250_fabae_no_outliers$larvae_total_number, col = "red")

# Number of rows before removing outliers
n_before <- nrow(data250_fabae)
print(n_before)

# Number of rows after removing outliers
n_after <- nrow(data250_fabae_no_outliers)
print(n_after)

library(emmeans)

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


