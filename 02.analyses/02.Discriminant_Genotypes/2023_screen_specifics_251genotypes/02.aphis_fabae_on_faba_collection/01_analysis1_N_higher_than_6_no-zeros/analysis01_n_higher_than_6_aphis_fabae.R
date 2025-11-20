rm(list = ls()) # Clear the environment
suppressPackageStartupMessages(library("tidyverse"))
suppressPackageStartupMessages(library("MASS"))
suppressPackageStartupMessages(library("lme4"))
suppressPackageStartupMessages(library("emmeans"))

##############
## Data import
##############

# 248 unique genotypes
# Keep only data if N>6
afabae_fecundity <- read.csv("01.data/2023_screen_specifics_251genotypes/02.aphis_fabae_on_faba_collection/aphis_fabae_fecundity_on_faba_collection_2024-04-02.csv", 
                            stringsAsFactors = F) %>% 
  mutate(type = gsub(pattern = "Equina", x = type, replacement = "Vicia faba var. equina")) %>% 
  mutate(type = gsub(pattern = "Minor", x = type, replacement = "Vicia faba var. minor")) %>% 
  mutate(type = gsub(pattern = "Major", x = type, replacement = "Vicia faba var. major")) %>% 
  mutate(type = gsub(pattern = "Paucijuga", x = type, replacement = "Vicia paucijuga")) %>% 
  #filter(adult_alive == 1) %>% # adult has to be alive for larvae number results to be kept
  filter(larvae_total_number != 0) 

genotypes_nrep_higher_than_6 = 
  afabae_fecundity %>% 
  dplyr::select(genotype, larvae_total_number) %>% 
  group_by(genotype) %>% 
  tally() %>% 
  filter(n >= 6) %>% # keep only if N plants => 6 
  pull(genotype)
afabae_fecundity_n6 <- afabae_fecundity %>% filter(genotype %in% genotypes_nrep_higher_than_6)
afabae_fecundity_n6 %>% pull(genotype) %>% unique() %>% length() # 193 genotypes kept out of 248

# Count number of biological replicate per genotype and plot
afabae_fecundity_n6 %>% 
  dplyr::select(genotype, larvae_total_number) %>% 
  group_by(genotype) %>% 
  tally() %>%
  ggplot(., aes(x = n)) +
  geom_histogram(fill = "lightgrey", color = "black", bins = 6) +
  ggtitle("Number of biological replicates per genotype (filtered, N>=6)") +
  labs(x = "Number of biological replicates", y = "Number of genotypes") + 
  scale_x_continuous(breaks = seq(6, 12, by = 1))
ggsave(filename = "02.analyses/2023_screen_specifics_251genotypes/02.aphis_fabae_on_faba_collection/01_analysis1_N_higher_than_6_no-zeros/plot00_n_reps_per_genotype.pdf", width = 8, height = 6) 

write.csv(afabae_fecundity_n6, 
          file = "02.analyses/2023_screen_specifics_251genotypes/02.aphis_fabae_on_faba_collection/01_analysis1_N_higher_than_6_no-zeros/aphis_fabae_fecundity_on_faba_collection_N_higher_than_6_no_zeros.csv", 
          row.names = F)
rm(afabae_fecundity)

###################################################
## Stat models: genotype effect on number of larvae
###################################################

# Continuous normally distributed data of the total number of larvae ?
afabae_fecundity_n6 %>% 
  ggplot(., aes(x = larvae_total_number)) +
  geom_histogram(aes(y = ..density..), fill = "lightgrey", color = "black", bins = 30) +
  geom_density() +
  ggtitle("Distribution of the number of A. fabae larvae (all genotypes)")

ggsave(filename = "02.analyses/2023_screen_specifics_251genotypes/02.aphis_fabae_on_faba_collection/01_analysis1_N_higher_than_6_no-zeros/plot01_distribution_larvae_number.pdf") 
ggsave(filename = "02.analyses/2023_screen_specifics_251genotypes/02.aphis_fabae_on_faba_collection/01_analysis1_N_higher_than_6_no-zeros/plot01_distribution_larvae_number.png") 


#############################
# Barplot of fecundity levels
#############################

median_larvae_counts_afabae = 
  afabae_fecundity_n6 %>% 
  group_by(ivip) %>% 
  summarise(median_larvae = median(larvae_total_number)) %>% 
  arrange(median_larvae) 

order_afabae <- median_larvae_counts_afabae %>% pull(ivip) %>% unique()

p1_afabae <- median_larvae_counts_afabae %>% 
  mutate(ivip = factor(ivip, levels = order_afabae)) %>% 
  ggplot(aes(x = ivip, y = median_larvae)) + 
  geom_bar(stat = "identity", fill="grey", alpha = 0.8, color = "black") +
  theme_bw() +
  theme(axis.text.x = element_blank(), 
        axis.title.x = element_text(size = 12), 
        axis.title.y = element_text(size = 12)) + 
  labs(x = "Vicia faba genotype", y = "Total number of larvae\n(from one single L1 larvae)") + 
  ggtitle("A. fabae fecundity on Vicia faba collection");p1_afabae

ggsave(plot = p1_afabae, width = 12, height = 8, 
       filename="02.analyses/2023_screen_specifics_251genotypes/02.aphis_fabae_on_faba_collection/01_analysis1_N_higher_than_6_no-zeros/plot02_afabae_barplot.pdf")
ggsave(plot = p1_afabae, width = 12, height = 8, units = "cm",
       filename="02.analyses/2023_screen_specifics_251genotypes/02.aphis_fabae_on_faba_collection/01_analysis1_N_higher_than_6_no-zeros/plot02_afabae_barplot.png")

#################
# 10 most extreme
#################
ten_most_resistant_genotypes_afabae <- 
  afabae_fecundity_n6 %>% 
  group_by(ivip) %>% 
  summarise(median_afabae = median(larvae_total_number)) %>% 
  arrange(median_afabae) %>% 
  head(n = 10) %>% 
  mutate(phenotype = "Ten most resistant")

ten_most_sensitive_genotypes_afabae <- 
  afabae_fecundity_n6 %>% 
  group_by(ivip) %>% 
  summarise(median_afabae = median(larvae_total_number)) %>% 
  arrange(median_afabae) %>% 
  tail(n = 10) %>%
  mutate(phenotype = "Ten most sensitive") 

# Number of biological replicate per genotype (to add to the graph)
n_per_genotype <- 
  afabae_fecundity_n6 %>% 
  dplyr::select(ivip, replicate) %>%
  group_by(ivip) %>% 
  tally() %>% 
  dplyr::rename("n_rep" = "n")

most_extreme_genotype_afabae <- bind_rows(ten_most_resistant_genotypes_afabae, ten_most_sensitive_genotypes_afabae) %>% 
  inner_join(n_per_genotype, by = "ivip") %>%
  mutate(n_rep = paste("n=",n_rep, sep = "")) %>% 
  dplyr::select(ivip, phenotype, n_rep)
  
p_most_contrasted_fabae <- 
  afabae_fecundity_n6 %>% 
  inner_join(., most_extreme_genotype_afabae) %>%
  mutate(ivip = factor(ivip, levels = order_afabae)) %>% 
  ggplot(., aes(x = ivip, y = larvae_total_number, fill = phenotype, label = paste(n_rep, sep = " "))) +
  geom_boxplot() +
  geom_point(alpha=0.3) +
  scale_fill_brewer(type = 'qual', palette = 1) +
  ggtitle("Most contrasted genotypes to Aphis fabae\n(N>=6 and zeros removed)") +
  coord_flip() +
  labs(x = "Vicia faba genotype", y = "Total number of Aphis fabae larvae (from one L1)") +
  scale_y_continuous(limits = c(0,90)) +
  geom_text(y = 75, hjust=0) + 
  theme(axis.text = element_text(size = 10));p_most_contrasted_fabae

ggsave(plot = p_most_contrasted_fabae, 
       filename = "02.analyses/2023_screen_specifics_251genotypes/02.aphis_fabae_on_faba_collection/01_analysis1_N_higher_than_6_no-zeros/plot03_afabae_most_contrasted_genotypes.pdf", 
       width = 10, height = 8)
ggsave(plot = p_most_contrasted_fabae, width = 10, height = 8, units = "cm",
       filename = "02.analyses/2023_screen_specifics_251genotypes/02.aphis_fabae_on_faba_collection/01_analysis1_N_higher_than_6_no-zeros/plot03_afabae_most_contrasted_genotypes.png")

#####################
## List of genotypes
#####################

afabae_fecundity_n6 %>% 
  inner_join(., most_extreme_genotype_afabae) %>% 
  dplyr::select(ivip, larvae_total_number, phenotype) %>% 
  group_by(ivip, phenotype) %>% 
  summarise(median_larvae = median(larvae_total_number, na.rm=T)) %>% 
  arrange(median_larvae, ) %>% 
  knitr::kable()

afabae_fecundity_n6 %>% 
  inner_join(., most_extreme_genotype_afabae) %>% 
  dplyr::select(ivip, larvae_total_number, phenotype) %>% 
  group_by(ivip, phenotype) %>% 
  summarise(median_larvae = median(larvae_total_number, na.rm=T)) %>% 
  arrange(median_larvae) %>% 
  write.csv(file = "02.analyses/2023_screen_specifics_251genotypes/02.aphis_fabae_on_faba_collection/01_analysis1_N_higher_than_6_no-zeros/most_extreme_genotypes_aphis_fabae.csv", 
            quote = F, row.names = F)

########
## Stats
########

# reverse a vector of genotype names
reverse_order_afabae <- rev(order_afabae) %>% intersect(most_extreme_genotype_afabae$ivip)

# Data of fecundity of the 20 most contrasted genotypes
df_fabae <- 
  afabae_fecundity_n6 %>% 
  inner_join(., most_extreme_genotype_afabae) %>% 
  dplyr::select(ivip, larvae_total_number, replicate) 
df_fabae$ivip = factor(x = df_fabae$ivip, levels = reverse_order_afabae)

# Fit a null and genotype GLM negative binomial model
fit0 = glm.nb(formula = larvae_total_number ~ 1, data = df_fabae)
fit1 = glm.nb(formula = larvae_total_number ~ 1 + ivip, data = df_fabae)
anova(fit0, fit1)
# Likelihood ratio tests of Negative Binomial Models
# 
# Response: larvae_total_number
# Model    theta Resid. df    2 x log-lik.   Test    df LR stat.      Pr(Chi)
# 1        1 4.794820       164       -1379.271                                   
# 2 1 + ivip 7.793472       145       -1315.745 1 vs 2    19 63.52587 1.057775e-06

# Fit a null (only replicate) and a genotype GLMM negative binomial model 
fit2 = lme4::glmer.nb(larvae_total_number ~ 1 + (1|ivip), data = df_fabae)
fit3 = lme4::glmer.nb(larvae_total_number ~ 1 + ivip + (1|replicate), data = df_fabae)
anova(fit2, fit3)
# Data: df_fabae
# Models:
#   fit2: larvae_total_number ~ 1 + (1 | ivip)
# fit3: larvae_total_number ~ 1 + ivip + (1 | replicate)
#         npar    AIC    BIC  logLik deviance  Chisq Df Pr(>Chisq)    
# fit2    3       1390.2 1399.6 -692.12   1384.2                         
# fit3    22      1379.6 1448.4 -667.82   1335.6 48.603 19  0.0002098 ***


AIC(fit0, fit1, fit2, fit3)
#       df      AIC
# fit0  2       1407.944
# fit1  21      1378.359
# fit2  3       1390.239
# fit3  22      1379.636

pairs(emmeans(fit3, "ivip")) %>% 
  broom::tidy() %>% 
  # filter rows that contains string 'Castell'
  dplyr::filter(., grepl("IVIP100336",contrast)) %>%
  dplyr::filter(adj.p.value < 0.05) %>% 
  knitr::kable()

# |term     |contrast        | null.value|  estimate| std.error|  df| statistic| adj.p.value|
# |:--------|:---------------|----------:|---------:|---------:|---:|---------:|-----------:|
# |ivip |JF2120 - JF2244 |          0| 0.7383705| 0.2065503| Inf|  3.574774|   0.0451584|
# |ivip |JF2120 - KF2221 |          0| 0.7349489| 0.1992506| Inf|  3.688566|   0.0306791|
# |ivip |JF2120 - KF0190 |          0| 0.7565193| 0.1945160| Inf|  3.889239|   0.0148198|
# |ivip |JF2120 - JF2236 |          0| 0.6869415| 0.1885369| Inf|  3.643539|   0.0358329|