require(lme4)

# Define the function
calculate_broad_heritability <- function(data) {
  
  # Fit the model with a fixed formula
  mod_h2 <- lme4::lmer(larvae_total_number ~ (1 | ivip) + (1 | session), data = data)
  
  # Get the variance components
  estimates <- as.data.frame(VarCorr(mod_h2))
  
  # Extract the relevant variance components
  genetic_variance <- estimates$vcov[1]
  residual_variance <- estimates$vcov[3]  # Adjust based on the placement of residual variance
  
  # Get the number of levels for the 'session' random effect
  nbE <- nlevels(data$session)
  
  # Calculate heritability (h2)
  h2 <- genetic_variance / (genetic_variance + (residual_variance / nbE))
  
  
  nlvl<-nlevels(data$ivip)
  
  nnotation<-nrow(data)
  
  # Return a list with heritability and variance components
  result <- list(
    heritability = h2,
    number_of_genotypes = nlvl,
    number_of_notation = nnotation)
  
}