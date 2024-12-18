setwd("~/Documents/annalyses/screening_choice_genotypes/A.Fabae")
data250fabae<-read.csv("aphis_fabae_fecundity_on_faba_collection_2024-04-02.csv", sep =",")

# Remove rows with missing values
data250fabae <- data250fabae %>% filter(!is.na(larvae_total_number))


data250fabae$short_code<-as.factor(data250fabae$short_code)
data250fabae$session<-as.factor(data250fabae$session)
data_filtered_250fabae <- data250fabae %>% filter(adult_alive != "0")
summary(data_filtered_250fabae)

data_filtered_250fabae <- data_filtered_250fabae %>%
  group_by(short_code) %>%           # Group by genotype
  filter(n() >= 6) %>%             # Keep only groups with 6 or more repetitions
  ungroup()                        # Ungroup to return to original data structure







#heritability----
mod_h2<-lmer(larvae_total_number~ (1|short_code)+(1|session), data=data_filtered_250fabae)
test <-anova(mod_h2)
summary(test)



estimates <- as.data.frame(VarCorr(mod_h2)) # get the Var. components

# extract the relevant data and put into df
genetic_variance <- estimates$vcov[1]
residual_variance <- estimates$vcov[3]#changer le numero en fonction du plaement de Residual

nbE<-nlevels(data_filtered_250fabae$session)

h2<-genetic_variance/(genetic_variance+(residual_variance/nbE))
print(h2)


