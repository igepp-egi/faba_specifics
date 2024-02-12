suppressPackageStartupMessages(library("tidyverse"))
library(multcompView)    
library("ggdist")

##############
## Data import
##############

afabae_fecundity <- read.csv("01.data/aphis_fabae_fecundity_on_faba_collection.csv", 
                            stringsAsFactors = F) %>% 
  mutate(type = gsub(pattern = "Equina", x = type, replacement = "Vicia faba var. equina")) %>% 
  mutate(type = gsub(pattern = "Minor", x = type, replacement = "Vicia faba var. minor")) %>% 
  mutate(type = gsub(pattern = "Major", x = type, replacement = "Vicia faba var. major")) %>% 
  mutate(type = gsub(pattern = "Paucijuga", x = type, replacement = "Vicia paucijuga")) %>% 
  filter(larvae_total_number != 0) %>% # Larvae probably died
  drop_na()

#############################################
# Number of biological replicate per genotype
#############################################
n_per_genotype_afabae <- 
  afabae_fecundity %>% 
  drop_na() %>% 
  filter(adult_alive != 0) %>%       # Keep only data with adult alive (from first L1 larvae deposited)
  select(genotype, replicate) %>%
  group_by(genotype) %>% 
  tally() %>% 
  dplyr::rename("n_rep" = "n")

aphis_fabae_n_rep_per_genotype <- 
  ggplot(n_per_genotype_afabae, aes(x = n_rep)) +
  geom_bar() + 
  scale_color_brewer(type = 'seq', palette = 1) +
  labs(x = "Number of biological replicates (plants)", y = "Genotype occurences") +
  ggtitle("Aphis fabae screen") +
  scale_x_continuous(breaks = seq(0,14,1))

ggsave(plot = aphis_fabae_n_rep_per_genotype, 
       filename = "02.analyses/01.acyrtosiphon_pisum_on_faba_collection/plot01_nb_rep_per_genotype.pdf")  

##############################
## Distribution of larvae data
##############################

afabae_fecundity %>% 
  ggplot(., aes(x = larvae_total_number)) +
  geom_histogram(aes(y = ..density..), fill = "lightgrey", color = "black", bins = 30) +
  geom_density() +
  ggtitle("Distribution of the number of Aphis fabae larvae across all genotypes")

###################
# Line plot with CI
###################

order_afabae <- afabae_fecundity %>% 
  group_by(genotype) %>% 
  summarise(median_afabae = median(larvae_total_number)) %>% 
  arrange(median_afabae) %>% 
  pull(genotype) %>% 
  unique()


afabae_fecundity %>% 
  mutate(genotype = factor(genotype, levels = order_afabae)) %>% 
  ggplot(aes(x = genotype, y = larvae_total_number)) + 
  stat_lineribbon(.width = c(0.50, 0.90, 0.95)) +
  scale_fill_brewer(palette = "Oranges") + 
  coord_flip() +
  theme(axis.text.y = element_text(size = 4)) +
  labs(x = "Vicia faba genotype", y = "Number of A. fabae larvae obtained from one L1 (Median -/+ CI)")

ggsave(filename = "02.analyses/02.aphis_fabae_on_faba_collection/plot01_aphis_fabae_ribbonplot.pdf")

##############################
# Median nb of A. fabae larvae
##############################

stat_res_table_aphis_fabae <- 
  afabae_fecundity %>% 
  select(genotype, replicate, larvae_total_number) %>% 
  group_by(genotype) %>% 
  summarise(median_afabae_larvae = median(larvae_total_number, na.rm = T), 
            sd_afabae_larvae = sd(larvae_total_number, na.rm = T)) %>% 
  arrange(median_afabae_larvae)

#########################################
# 10 most resistant / sensitive genotypes
#########################################

ten_most_res_genotypes <- 
  stat_res_table_aphis_fabae %>% 
  tail(n = 10) %>% 
  mutate(phenotype = "Ten most resistant")

ten_most_sensitive_genotypes <- 
  stat_res_table_aphis_fabae %>% 
  head(n = 10) %>% 
  mutate(phenotype = "Ten most sensitive") 

most_extreme_genotype_afabae <- bind_rows(ten_most_res_genotypes, ten_most_sensitive_genotypes) %>% 
  select(genotype, phenotype)

afabae_fecundity %>% 
  inner_join(., most_extreme_genotype_afabae) %>%
  mutate(genotype = factor(genotype, levels = order_afabae)) %>% 
  ggplot(., aes(x = genotype, y = larvae_total_number, fill = phenotype)) +
  geom_boxplot() +
  geom_point(alpha=0.3) +
  scale_fill_brewer(type = 'qual', palette = 2) +
  ggtitle("Most contrasted genotypes\n (Aphis fabae screen)") +
  coord_flip() +
  labs(x = "Vicia faba genotype", y = "Total number of Aphis fabae larvae (from one L1)") +
  scale_y_continuous(limits = c(0,80))

ggsave(filename = "02.analyses/02.aphis_fabae_on_faba_collection/plot02_aphis_fabae_10_most_contrasted_genotypes.pdf")
