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


#mod<-gamlss(larvae_total_number ~ ivip + session/block, data = data250_pisum, family = SEP2)
#mod<-gamlss(larvae_total_number ~ ivip + session/block, data = data250_pisum, family = SE)
for (i in 1) {
  
mod0<-lm(larvae_total_number ~ ivip + session/block, data = data250_pisum)
print("0")
#mod1<-gamlss(larvae_total_number ~ ivip + session, data = data250_pisum, family = ST3)
print("1")
mod2<-gamlss(larvae_total_number ~ ivip + random(session), data = data250_pisum, family = ST3)
print("2")
#mod3<-gamlss(larvae_total_number ~ ivip + session/block, data = data250_pisum, family = ST3)
print("3")
mod4<-gamlss(larvae_total_number ~ ivip, random = ~ session/block, data = data250_pisum, family = ST3)
print("4")
}

compare_performance(mod0, mod2, mod4, rank = FALSE, verbose = TRUE)




fitted_values<-fitted(mod2)

# Extract fitted values and residuals
residuals_values <- residuals(mod2)

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







# Extract the fitted values (adjusted means)
fitted_values <- fitted(model)

# Extract the random effects
random_effects <- ranef(model)








coef_df <- as.data.frame(coef(mod2))
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









