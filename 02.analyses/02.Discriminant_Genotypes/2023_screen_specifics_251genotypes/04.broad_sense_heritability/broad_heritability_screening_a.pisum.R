#Packages ----
library(dplyr)
library(tidyr)
library(ggplot2)
library(lme4)
library(emmeans)
library(pbkrtest)
library(randtests)
#data preparation ----
setwd("~/Documents/annalyses/screening_choice_genotypes/A.Pissum")
data250pisum<-read.csv("acyrthospipon_pisum_fecundity_on_faba_collection.csv", sep =",")

data250pisum$short_code<-as.factor(data250pisum$short_code)
data250pisum$session<-as.factor(data250pisum$session)
data250pisum <- data250pisum %>% filter(!is.na(larvae_total_number))
data_filtered_250_pisum <- data250pisum %>% filter(adult_alive != "0")
summary(data_filtered_250_pisum)

data_filtered_250_pisum <- data_filtered_250_pisum %>%
  group_by(short_code) %>%           # Group by genotype
  filter(n() >= 6) %>%             # Keep only groups with 6 or more repetitions
  ungroup()                        # Ungroup to return to original data structure


#heritability----
mod_h2<-lmer(larvae_total_number~ (1|short_code)+(1|session), data=data_filtered_250_pisum)



estimates <- as.data.frame(VarCorr(mod_h2)) # get the Var. components

# extract the relevant data and put into df
genetic_variance <- estimates$vcov[1]
residual_variance <- estimates$vcov[3]#changer le numero en fonction du plaement de Residual

nbE<-nlevels(data_filtered_250_pisum$session)

h2<-genetic_variance/(genetic_variance+(residual_variance/nbE))
print(h2)
heritability <- c(heritability,h2 )

var_components <- VarCorr(mod_h2)
