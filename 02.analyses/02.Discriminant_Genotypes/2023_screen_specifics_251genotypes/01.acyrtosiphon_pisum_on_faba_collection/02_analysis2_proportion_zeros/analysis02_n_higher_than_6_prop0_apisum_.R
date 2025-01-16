suppressPackageStartupMessages(library("tidyverse"))
library(multcompView)    
library("ggdist")

##############
## Data import
##############

### Genetic structure info (Jonathan, UMR Agroécologie/ECP)
igepp_to_ivip <- read.csv("01.data/03.genetic_structure_jonathan_dijon/faba_bean_ivip_ids.csv") %>% 
  rename("ivip" = "seedbank_id") %>% 
  rename("genotype" = "igepp_id")

genetic_structure <- read.csv("01.data/03.genetic_structure_jonathan_dijon/genetic_groups_K_equal_to_6.csv") %>% 
  rename("ivip" = "genotype")

igepp_ivip_genetic_group <- inner_join(igepp_to_ivip, genetic_structure)
rm(igepp_to_ivip, genetic_structure)

### Fecundity data
# 248 unique genotypes
# Keep only data if N>6
apisum_fecundity <- read.csv("01.data/01.acyrthosiphon_pisum_P123amp_on_faba_collection/acyrthospipon_pisum_fecundity_on_faba_collection.csv", 
                            stringsAsFactors = F) %>% 
  mutate(type = gsub(pattern = "Equina", x = type, replacement = "Vicia faba var. equina")) %>% 
  mutate(type = gsub(pattern = "Minor", x = type, replacement = "Vicia faba var. minor")) %>% 
  mutate(type = gsub(pattern = "Major", x = type, replacement = "Vicia faba var. major")) %>% 
  mutate(type = gsub(pattern = "Paucijuga", x = type, replacement = "Vicia paucijuga")) 

genotypes_nrep_higher_than_6 = 
  apisum_fecundity %>% 
  drop_na() %>% 
  count(genotype) %>% 
  filter(n >= 6) %>% # keep only if N plants => 6 >> results 223 genotypes out of 248
  pull(genotype)

apisum_fecundity_n6 <- apisum_fecundity %>% filter(genotype %in% genotypes_nrep_higher_than_6)
rm(apisum_fecundity)

###################################################
# Proportion of 0 from total number of measurements
###################################################

apisum_fecundity_n6_prop_zeros <- 
  apisum_fecundity_n6 %>% 
  select(genotype, larvae_total_number) %>%
  drop_na() %>% 
  mutate(condition = ifelse(larvae_total_number > 0, "higher_than_zeros","null")) %>% 
  select(genotype, condition) %>% 
  group_by(genotype) %>% 
  summarise(n_positive = length(condition == "higher_than_zeros"), n_values = n())

apisum_fecundity_n6_prop_zeros <- 
  apisum_fecundity_n6 %>% 
  select(genotype, larvae_total_number) %>%
  drop_na() %>% 
  group_by(genotype) %>% 
  summarise(n_zeros = sum(larvae_total_number == 0, na.rm = TRUE), 
            n_pos = sum(larvae_total_number > 0, na.rm = TRUE)) %>% 
  ungroup() %>% 
  mutate(freq_zeros = n_zeros / (n_zeros + n_pos) * 100) %>% 
  mutate(n_values = n_zeros + n_pos)

summary(apisum_fecundity_n6_prop_zeros$freq_zeros)
# Min.    1st Qu.  Median    Mean    3rd Qu.    Max. 
# 0.000   0.000    8.333     10.837  16.667     41.667 

p_distri_zeros <- 
  apisum_fecundity_n6_prop_zeros %>% 
  drop_na() %>% 
  ggplot(., aes(x = freq_zeros)) +
  geom_histogram(color = "black", binwidth = 5) +
  ggtitle("Distribution of the number of true zeros of A. pisum larvae across all genotypes") +
  labs(x = "Proportion of true zeros among valid values (%)", y = "Number of genotypes (total: 217 genotypes)")

p_distri_zeros

ggsave(plot = p_distri_zeros, 
       filename = "02.analyses/01.acyrtosiphon_pisum_on_faba_collection/02_analysis2_proportion_zeros/plot01_distribution_zeros.pdf") 
ggsave(plot = p_distri_zeros, 
       filename = "02.analyses/01.acyrtosiphon_pisum_on_faba_collection/02_analysis2_proportion_zeros/plot01_distribution_zeros.png") 

###################################
# Plot genotypes of the 90th percentile 
###################################

# 90th percentile
selected_trehsold = quantile(apisum_fecundity_n6_prop_zeros$freq_zeros, probs = c(0.9)) %>% as.numeric()

order_genotypes <- apisum_fecundity_n6_prop_zeros %>% 
  filter(freq_zeros > selected_trehsold) %>% 
  arrange(desc(freq_zeros)) %>% 
  pull(genotype)

p_zeros = 
  apisum_fecundity_n6_prop_zeros %>% 
  inner_join(., igepp_ivip_genetic_group) %>% 
  filter(freq_zeros > selected_trehsold) %>% 
  arrange(desc(freq_zeros)) %>% 
  mutate(genotype = factor(genotype, levels = order_genotypes)) %>% 
  ggplot(aes(x = genotype, y = freq_zeros, fill = clusterName, label = paste("N_zeros:", n_zeros,"/",n_values))) + 
  geom_bar(stat = "identity") + 
  coord_flip() +
  scale_y_continuous(limits = c(0,50)) +
  geom_text(aes(x = genotype, y = 45), size = 3) +
  labs(x = "Vicia faba genotype", y = "Proportion of true zeros (%)") +
  ggtitle("Proportions of true zeros (A. pisum screen)") +
  scale_fill_brewer(type = "qual", palette = 2) +
  theme(axis.text = element_text(size = 10, color = "black")) 
p_zeros

ggsave(plot = p_zeros, 
       filename = "02.analyses/01.acyrtosiphon_pisum_on_faba_collection/02_analysis2_proportion_zeros/plot02_number_of_true_zeros.pdf") 
ggsave(plot = p_zeros, 
       filename = "02.analyses/01.acyrtosiphon_pisum_on_faba_collection/02_analysis2_proportion_zeros/plot02_number_of_true_zeros.png") 

############################
# List of genotypes to keep
###########################
apisum_fecundity_n6_prop_zeros %>% 
  inner_join(., igepp_ivip_genetic_group) %>% 
  filter(freq_zeros > selected_trehsold) %>% 
  select(- cluster, - species, - n_values, - lieu, - ivip, - clusterName) %>% 
  knitr::kable()
