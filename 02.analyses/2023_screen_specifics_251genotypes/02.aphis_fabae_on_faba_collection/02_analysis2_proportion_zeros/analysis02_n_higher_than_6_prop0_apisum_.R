suppressPackageStartupMessages(library("tidyverse"))
#library(multcompView)    
#library("ggdist")

##############
## Data import
##############

### Genetic structure info (Jonathan, UMR Agroécologie/ECP)
igepp_to_ivip <- read.csv("00.info/02.genetic_structure_jonathan_dijon/faba_bean_ivip_ids.csv") %>% 
  rename("ivip" = "seedbank_id") %>% 
  rename("genotype" = "igepp_id")

genetic_structure <- read.csv("00.info/02.genetic_structure_jonathan_dijon/genetic_groups_K_equal_to_6.csv") %>% 
  rename("ivip" = "genotype")

igepp_ivip_genetic_group <- inner_join(igepp_to_ivip, genetic_structure)
rm(igepp_to_ivip, genetic_structure)

### Fecundity data
# 222 unique genotypes
# Keep only data if N>6
afabae_fecundity <- read.csv("01.data/2023_screen_specifics_251genotypes/02.aphis_fabae_on_faba_collection/aphis_fabae_fecundity_on_faba_collection_2024-04-02.csv", 
                            stringsAsFactors = F) %>% 
  mutate(type = gsub(pattern = "Equina", x = type, replacement = "Vicia faba var. equina")) %>% 
  mutate(type = gsub(pattern = "Minor", x = type, replacement = "Vicia faba var. minor")) %>% 
  mutate(type = gsub(pattern = "Major", x = type, replacement = "Vicia faba var. major")) %>% 
  mutate(type = gsub(pattern = "Paucijuga", x = type, replacement = "Vicia paucijuga")) 

genotypes_nrep_higher_than_6 = 
  afabae_fecundity %>% 
  drop_na() %>% 
  count(genotype) %>% 
  filter(n >= 6) %>% # keep only if N plants => 6 >> results = 207 genotypes out of 222
  pull(genotype)

afabae_fecundity_n6 <- afabae_fecundity %>% filter(genotype %in% genotypes_nrep_higher_than_6)
rm(afabae_fecundity, genotypes_nrep_higher_than_6)

write.csv(x = afabae_fecundity_n6, 
          file = "01.data/2023_screen_specifics_251genotypes/02.aphis_fabae_on_faba_collection/aphis_fabae_fecundity_on_faba_collection_n_higher_6.csv", 
          row.names = FALSE)

###################################################
# Proportion of 0 from total number of measurements
###################################################

# calculate proportions of 0 per genotype
afabae_fecundity_n6_prop_zeros <- 
  afabae_fecundity_n6 %>% 
  select(genotype, larvae_total_number) %>%
  drop_na() %>% 
  group_by(genotype) %>% 
  summarise(n_zeros = sum(larvae_total_number == 0, na.rm = TRUE), 
            n_pos = sum(larvae_total_number > 0, na.rm = TRUE)) %>% 
  ungroup() %>% 
  mutate(freq_zeros = n_zeros / (n_zeros + n_pos) * 100) %>% 
  mutate(n_values = n_zeros + n_pos)

summary(afabae_fecundity_n6_prop_zeros$freq_zeros)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.00   16.67   25.00   23.82   33.33   72.73 

p_distri_zeros <- 
  afabae_fecundity_n6_prop_zeros %>% 
  drop_na() %>% 
  ggplot(., aes(x = freq_zeros)) +
  geom_histogram(aes(y = after_stat(density)), colour = 1, fill = "white", bins = 20) +
  geom_density() +
  scale_x_continuous(breaks = c(0,10,20,30,40,50,60,70,80,90,100)) +
  ggtitle("Distribution of the number of true zeros of A. fabae larvae across all genotypes") +
  labs(x = "Proportion of true zeros among valid values (%)", y = "Number of genotypes (total: 207 genotypes)")
p_distri_zeros

ggsave(plot = p_distri_zeros, 
       filename = "02.analyses/2023_screen_specifics_251genotypes/02.aphis_fabae_on_faba_collection/02_analysis2_proportion_zeros/plot01_distribution_zeros.pdf") 
ggsave(plot = p_distri_zeros, 
       filename = "02.analyses/2023_screen_specifics_251genotypes/02.aphis_fabae_on_faba_collection/02_analysis2_proportion_zeros/plot01_distribution_zeros.png") 

###################################
# Plot genotypes of the 90th percentile 
###################################

# 90th percentile
selected_treshold = quantile(afabae_fecundity_n6_prop_zeros$freq_zeros, probs = c(0.9)) %>% as.numeric()

order_genotypes <- afabae_fecundity_n6_prop_zeros %>% 
  filter(freq_zeros > selected_treshold) %>% 
  arrange(desc(freq_zeros)) %>% 
  pull(genotype)

p_zeros = 
  afabae_fecundity_n6_prop_zeros %>% 
  inner_join(., igepp_ivip_genetic_group) %>% 
  filter(freq_zeros > selected_treshold) %>% 
  mutate(cluster = as.factor(cluster)) %>% 
  mutate(genotype = factor(genotype, levels = order_genotypes)) %>% 
  ggplot(aes(x = genotype, y = freq_zeros, fill = cluster, label = paste("N_zeros:", n_zeros,"/",n_values))) + 
  geom_bar(stat = "identity") + 
  coord_flip() +
  scale_y_continuous(limits = c(0,100)) +
  geom_text(aes(x = genotype, y = 90), size = 3) +
  labs(x = "Vicia faba genotype", y = "Proportion of true zeros (%)") +
  ggtitle("Proportions of true zeros (A. fabae AF001)") +
  scale_fill_brewer(type = "qual", palette = 2) +
  theme(axis.text = element_text(size = 10, color = "black")) 
p_zeros

ggsave(plot = p_zeros, 
       filename = "02.analyses/2023_screen_specifics_251genotypes/02.aphis_fabae_on_faba_collection/02_analysis2_proportion_zeros/plot02_number_of_true_zeros.pdf") 
ggsave(plot = p_zeros, 
       filename = "02.analyses/2023_screen_specifics_251genotypes/02.aphis_fabae_on_faba_collection/02_analysis2_proportion_zeros/plot02_number_of_true_zeros.png") 

############################
# List of genotypes to keep
###########################
afabae_fecundity_n6_prop_zeros %>% 
  inner_join(., igepp_ivip_genetic_group) %>% 
  filter(freq_zeros > selected_treshold) %>% 
  arrange(freq_zeros) %>% 
  select(- cluster, - species, - n_values, - lieu, - ivip, - clusterName) %>% 
  knitr::kable()

#################################################
# List of genotypes to use as controls (no zeros)
#################################################

afabae_fecundity_n6_prop_zeros %>% 
  filter(freq_zeros == 0) %>% 
  filter(n_values == 12) %>% 
  knitr::kable()
