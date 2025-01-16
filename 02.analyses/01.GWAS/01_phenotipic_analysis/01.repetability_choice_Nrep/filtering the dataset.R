

# Define the function
analyse_outpout <- function(data) {
  
  # Fit the model with a fixed formula
  mod_h2 <- lmer(larvae_total_number ~ (1 | short_code) + (1 | session), data = data)
  
  # Get the variance components
  estimates <- as.data.frame(VarCorr(mod_h2))
  
  # Extract the relevant variance components
  genetic_variance <- estimates$vcov[1]
  residual_variance <- estimates$vcov[3]  # Adjust based on the placement of residual variance
  
  # Get the number of levels for the 'session' random effect
  nbE <- nlevels(data$session)
  
  # Calculate heritability (h2)
  h2 <- genetic_variance / (genetic_variance + (residual_variance / nbE))
  

  nlvl<-nlevels(data$genotype)
  
  # Return a list with heritability and variance components
  result <- list(
    heritability = h2,
    number_of_genotypes = nlvl)
  
  return(result)
}

# Data preparation ----
data250pisum_00 <- read.csv("acyrthospipon_pisum_fecundity_on_faba_collection.csv") %>% 
  mutate(genotype = as.factor(genotype)) %>% mutate(session = as.factor(session)) %>% mutate(short_code = as.factor(short_code))
analyse_outpout(data = data250pisum_00)


data250pisum_01 <- data250pisum_00 %>% filter(!is.na(larvae_total_number))
data250pisum_01$genotype <- droplevels(data250pisum_01$genotype)
analyse_outpout(data = data250pisum_01)


data250pisum_02<- data250pisum_01 %>% filter(adult_alive != "0") #remove if no adult
data250pisum_02$genotype <- droplevels(data250pisum_02$genotype)
analyse_outpout(data = data250pisum_02)


# Step 0: Filter the data if n replicates < 6 ----
data250pisum_03 <- data250pisum_02 %>%
  group_by(short_code) %>%           # Group by genotype
  filter(n() >= 6) %>%             # Keep only groups with 6 or more repetitions
  ungroup()  
data250pisum_03$genotype <- droplevels(data250pisum_03$genotype)
analyse_outpout(data = data250pisum_03)
# we lose genotypes therfore statistical power, 
#this does not lead to an increase in heritability (almost the opposite)
#less repetitions is not a problem because the snps are the one we want replicated 
#if we have repetition for the different haplotype we don't need repetition from
#the exact same plant. GWAS was developped for animals and not plant 
#with animals we dont have many times the same animal or even human but it is not a problem