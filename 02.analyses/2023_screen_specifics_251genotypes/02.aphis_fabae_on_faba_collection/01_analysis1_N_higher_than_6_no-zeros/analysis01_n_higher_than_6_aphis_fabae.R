suppressPackageStartupMessages(library("tidyverse"))

##############
## Data import
##############

###########################################################
### Genetic structure info (Jonathan, UMR Agroécologie/ECP)
###########################################################@
igepp_to_ivip <- read.csv("00.info/02.genetic_structure_jonathan_dijon/faba_bean_ivip_ids.csv") %>% 
  rename("ivip" = "seedbank_id") %>% 
  rename("genotype" = "igepp_id")

genetic_structure <- read.csv("00.info/02.genetic_structure_jonathan_dijon/genetic_groups_K_equal_to_6.csv") %>% 
  rename("ivip" = "genotype")

igepp_ivip_genetic_group <- inner_join(igepp_to_ivip, genetic_structure)
rm(igepp_to_ivip, genetic_structure)

###################
### Fecundity data
##################
# 248 unique genotypes
# Keep only data if N>6
afabae_fecundity <- read.csv("01.data/2023_screen_specifics_251genotypes/02.aphis_fabae_on_faba_collection/aphis_fabae_fecundity_on_faba_collection_2024-04-02.csv", 
                            stringsAsFactors = F) %>% 
  mutate(type = gsub(pattern = "Equina", x = type, replacement = "Vicia faba var. equina")) %>% 
  mutate(type = gsub(pattern = "Minor", x = type, replacement = "Vicia faba var. minor")) %>% 
  mutate(type = gsub(pattern = "Major", x = type, replacement = "Vicia faba var. major")) %>% 
  mutate(type = gsub(pattern = "Paucijuga", x = type, replacement = "Vicia paucijuga")) %>% 
  filter(larvae_total_number != 0) # Larvae probably died

genotypes_nrep_higher_than_6 = 
  afabae_fecundity %>% 
  dplyr::select(genotype, larvae_total_number) %>% 
  group_by(genotype) %>% 
  tally() %>% 
  filter(n >= 6) %>% # keep only if N plants => 6 >> results 223 genotypes out of 248
  pull(genotype)

afabae_fecundity_n6 <- afabae_fecundity %>% filter(genotype %in% genotypes_nrep_higher_than_6)
rm(afabae_fecundity)

###################################################
## Stat models: genotype effect on number of larvae
###################################################

# Continuous normally distributed data of the total number of larvae ?
afabae_fecundity_n6 %>% 
  ggplot(., aes(x = larvae_total_number)) +
  geom_histogram(aes(y = ..density..), fill = "lightgrey", color = "black", bins = 30) +
  geom_density() +
  ggtitle("Distribution of the number of A. fabae larvae across all genotypes")

ggsave(filename = "02.analyses/2023_screen_specifics_251genotypes/02.aphis_fabae_on_faba_collection/01_analysis1_N_higher_than_6_no-zeros/plot01_distribution_larvae_number.pdf") 
ggsave(filename = "02.analyses/2023_screen_specifics_251genotypes/02.aphis_fabae_on_faba_collection/01_analysis1_N_higher_than_6_no-zeros/plot01_distribution_larvae_number.png") 


#############################
# Barplot of fecundity levels
#############################

median_larvae_counts_afabae = 
  afabae_fecundity_n6 %>% 
  group_by(genotype) %>% 
  summarise(median_larvae = median(larvae_total_number)) %>% 
  arrange(median_larvae) 

order_afabae <- median_larvae_counts_afabae %>% pull(genotype) %>% unique()

p1_afabae <- median_larvae_counts_afabae %>% 
  mutate(genotype = factor(genotype, levels = order_afabae)) %>% 
  ggplot(aes(x = genotype, y = median_larvae)) + 
  geom_bar(stat = "identity", fill="grey", alpha = 0.8, color = "black") +
  theme_bw() +
  theme(axis.text.x = element_blank(), 
        axis.title.x = element_text(size = 12), 
        axis.title.y = element_text(size = 12)) + 
  labs(x = "Vicia faba genotype", y = "Total number of larvae\n(from one single N1 larvae)") + 
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
  group_by(genotype) %>% 
  summarise(median_afabae = median(larvae_total_number)) %>% 
  arrange(median_afabae) %>% 
  head(n = 10) %>% 
  mutate(phenotype = "Ten most resistant")

ten_most_sensitive_genotypes_afabae <- 
  afabae_fecundity_n6 %>% 
  group_by(genotype) %>% 
  summarise(median_afabae = median(larvae_total_number)) %>% 
  arrange(median_afabae) %>% 
  tail(n = 10) %>%
  mutate(phenotype = "Ten most sensitive") 

# Number of biological replicate per genotype (to add to the graph)
n_per_genotype <- 
  afabae_fecundity_n6 %>% 
  dplyr::select(genotype, replicate) %>%
  group_by(genotype) %>% 
  tally() %>% 
  dplyr::rename("n_rep" = "n")

most_extreme_genotype_afabae <- bind_rows(ten_most_resistant_genotypes_afabae, ten_most_sensitive_genotypes_afabae) %>% 
  inner_join(n_per_genotype, by = "genotype") %>%
  mutate(n_rep = paste("n=",n_rep, sep = "")) %>% 
  inner_join(igepp_ivip_genetic_group, by = "genotype") %>% 
  dplyr::select(genotype, phenotype, n_rep, clusterName, cluster) %>% # add n replicates + genetic cluster info
  mutate(cluster = paste("group",cluster,sep = ""))
  

most_contrasted_fabae <- 
  afabae_fecundity_n6 %>% 
  inner_join(., most_extreme_genotype_afabae) %>%
  mutate(genotype = factor(genotype, levels = order_afabae)) %>% 
  ggplot(., aes(x = genotype, y = larvae_total_number, fill = phenotype, label = paste(n_rep, sep = " "))) +
  geom_boxplot() +
  geom_point(alpha=0.3) +
  scale_fill_brewer(type = 'qual', palette = 1) +
  ggtitle("Most contrasted genotypes to Aphis fabae\n(N>=6 and zeros removed)") +
  coord_flip() +
  labs(x = "Vicia faba genotype", y = "Total number of Aphis fabae larvae (from one L1)") +
  scale_y_continuous(limits = c(0,90)) +
  geom_text(y = 75, hjust=0) + 
  theme(axis.text = element_text(size = 10));most_contrasted_fabae

ggsave(plot = most_contrasted_fabae, 
       filename = "02.analyses/2023_screen_specifics_251genotypes/02.aphis_fabae_on_faba_collection/01_analysis1_N_higher_than_6_no-zeros/plot03_afabae_most_contrasted_genotypes.pdf", 
       width = 10, height = 8)
ggsave(plot = most_contrasted_fabae, width = 10, height = 8, units = "cm",
       filename = "02.analyses/2023_screen_specifics_251genotypes/02.aphis_fabae_on_faba_collection/01_analysis1_N_higher_than_6_no-zeros/plot03_afabae_most_contrasted_genotypes.png")

#####################
## List of genotypes
#####################

afabae_fecundity_n6 %>% 
  inner_join(., most_extreme_genotype_afabae) %>% 
  dplyr::select(genotype, larvae_total_number, phenotype) %>% 
  group_by(genotype, phenotype) %>% 
  summarise(median_larvae = median(larvae_total_number, na.rm=T)) %>% 
  arrange(median_larvae) %>% 
  knitr::kable()

afabae_fecundity_n6 %>% 
  inner_join(., most_extreme_genotype_afabae) %>% 
  dplyr::select(genotype, larvae_total_number, phenotype) %>% 
  group_by(genotype, phenotype) %>% 
  summarise(median_larvae = median(larvae_total_number, na.rm=T)) %>% 
  arrange(median_larvae) %>% 
  write.csv(file = "02.analyses/2023_screen_specifics_251genotypes/02.aphis_fabae_on_faba_collection/01_analysis1_N_higher_than_6_no-zeros/most_extreme_genotypes_aphis_fabae.csv", 
            quote = F, row.names = F)

