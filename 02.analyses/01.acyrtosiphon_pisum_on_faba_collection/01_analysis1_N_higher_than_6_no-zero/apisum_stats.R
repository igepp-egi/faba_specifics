library("lme4")
library("MASS")
library("glmmTMB")

source("02.analyses/01.acyrtosiphon_pisum_on_faba_collection/01_analysis1_N_higher_than_6_no-zero/analysis01_n_higher_than_6_apisum_.R")

########################################################
# Genotype levels: ordered from highest larvae to lowest
########################################################

genotype_order <- 
  apisum_fecundity_n6 %>% 
  group_by(genotype) %>% 
  summarise(median_larvae = median(larvae_total_number, na.rm=TRUE)) %>% 
  arrange(desc(median_larvae)) %>% 
  pull(genotype) %>% 
  unique()

apisum_fecundity_n6$genotype = factor(apisum_fecundity_n6$genotype, levels = genotype_order)
  
#######
### Stat
########
model0 <- glm.nb(formula = larvae_total_number ~ 1, data = apisum_fecundity_n6)
model1 <- MASS::glm.nb(formula = larvae_total_number ~ genotype, data = apisum_fecundity_n6) # sans effet aleatoire
#model2 <- lme4::glmer.nb(formula = larvae_total_number ~ genotype + (1|genotype/replicate), data=apisum_fecundity_n6)
model3 <- glmmTMB::glmmTMB(formula = larvae_total_number ~ genotype + (1|genotype/replicate), data=apisum_fecundity_n6)
              
              
         
                 
apisum_fecundity_n6
