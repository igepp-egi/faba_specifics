rm(list = ls()) # clean workspace
suppressPackageStartupMessages(library("tidyverse"))
suppressPackageStartupMessages(library("patchwork"))
suppressPackageStartupMessages(library("ggrepel"))
suppressPackageStartupMessages(library("ggpubr"))

###################################################
## Execute single analysis of A. pisum and A. fabae
## Analysis 01: without zeros and N>=6
##################################################

#source("02.analyses/2023_screen_specifics_251genotypes/01.acyrtosiphon_pisum_on_faba_collection/01_analysis1_N_higher_than_6_no-zero/analysis01_n_higher_than_6_apisum_.R")
#source("02.analyses/2023_screen_specifics_251genotypes/02.aphis_fabae_on_faba_collection/01_analysis1_N_higher_than_6_no-zeros/analysis01_n_higher_than_6_aphis_fabae.R")

apisum_fecundity_n6 <- read.csv(
  file = "02.analyses/02.Discriminant_Genotypes/2023_screen_specifics_251genotypes/01.acyrtosiphon_pisum_on_faba_collection/01_analysis1_N_higher_than_6_no-zero/apisum_fecundity_on_faba_n_higher_than_6_no_zeros.csv", 
  stringsAsFactors = F)
most_extreme_genotype_apisum <- read.csv(
  file = "02.analyses/02.Discriminant_Genotypes/2023_screen_specifics_251genotypes/01.acyrtosiphon_pisum_on_faba_collection/01_analysis1_N_higher_than_6_no-zero/apisum_most_contrasted_genotypes.csv", 
  stringsAsFactors = F)

afabae_fecundity_n6 <- read.csv(
  file = "02.analyses/02.Discriminant_Genotypes/2023_screen_specifics_251genotypes/02.aphis_fabae_on_faba_collection/01_analysis1_N_higher_than_6_no-zeros/aphis_fabae_fecundity_on_faba_collection_N_higher_than_6_no_zeros.csv", 
  stringsAsFactors = F)
most_extreme_genotype_afabae = read.csv(
  file = "02.analyses/02.Discriminant_Genotypes/2023_screen_specifics_251genotypes/02.aphis_fabae_on_faba_collection/01_analysis1_N_higher_than_6_no-zeros/most_extreme_genotypes_aphis_fabae.csv", 
  stringsAsFactors = F)

#################################
# Paired plots (absolute numbers)
#################################

# 222 genotypes
res_pisum <- 
  apisum_fecundity_n6 %>% 
  group_by(ivip) %>% 
  summarise(median_apisum = median(larvae_total_number)) %>% 
  arrange(median_apisum) 

# 191 genotypes
res_afabae <- 
  afabae_fecundity_n6 %>% 
  group_by(ivip) %>% 
  summarise(median_afabae = median(larvae_total_number)) %>% 
  arrange(median_afabae)

# 190 genotypes communs
res_pisum_and_afabae <- inner_join(res_pisum, res_afabae, by = "ivip") %>% 
  rename("A_pisum" = "median_apisum") %>% 
  rename("A_fabae" = "median_afabae") 

genotypes_res_afabae <- most_extreme_genotype_afabae %>% filter(phenotype == "Ten most resistant") %>% pull(ivip)
genotypes_res_pisum <- most_extreme_genotype_apisum %>% filter(phenotype == "Ten most resistant") %>% pull(ivip)
all_genotypes = c(genotypes_res_afabae, genotypes_res_pisum)

###############
### Scatterplot
###############
points_to_label <- 
  res_pisum_and_afabae %>% 
  filter(ivip %in% all_genotypes) 

additional_labels <- res_pisum_and_afabae %>% filter(A_pisum >= 35) %>% filter(A_fabae >= 40)

points_to_label2 <- bind_rows(points_to_label, additional_labels)

p_scatterplot_continuous <- 
  res_pisum_and_afabae %>% 
  ggplot(., aes(x = A_fabae, y = A_pisum, label = ivip)) +
  geom_point() +
  labs(x = "Number of larvae from one L1 larvae of A. fabae") + 
  labs(y = "Number of larvae from one L1 larvae of A. pisum") +
  geom_label_repel(data = points_to_label2, max.overlaps = 10) +
  scale_fill_brewer(type = "qual", palette = 2) +
  theme_bw() +
  # add regression coefficient equation on the graph
  geom_smooth(method = "lm", se = FALSE, color = "black", formula = y ~ x) +
  # add equation and R sqaured
  stat_cor(method = "pearson", label.x = 20, label.y = 40, 
           aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")), 
           label.x.npc = "right", label.y.npc = "top", size = 4);p_scatterplot_continuous

ggsave(plot = p_scatterplot_continuous, 
       filename = "02.analyses/02.Discriminant_Genotypes/2023_screen_specifics_251genotypes/03.comparison_apisum_vs_afabae/plot01_scatterplot_apisum_vs_afabae.pdf", 
       width = 7, height = 7)
ggsave(plot = p_scatterplot_continuous, width = 10, height = 10, units = "cm",
       filename = "02.analyses/02.Discriminant_Genotypes/2023_screen_specifics_251genotypes/03.comparison_apisum_vs_afabae/plot01_scatterplot_apisum_vs_afabae.png")

###########################################
# Closer look at some interesting genotypes
###########################################


#interesting_genotypes = c("JF2094", "JF2274", "JF2256", "KF2219")
interesting_genotypes = points_to_label2$genotype

afabae_fecundity_n6_subset = afabae_fecundity_n6 %>% dplyr::select(genotype, aphid_species, larvae_total_number) %>% rename("aphid" = "aphid_species")
apisum_fecundity_n6_subset = apisum_fecundity_n6 %>% dplyr::select(genotype, aphid, larvae_total_number) 

afabae_apisum_fecundity_contrasted_genotypes = bind_rows(afabae_fecundity_n6_subset, apisum_fecundity_n6_subset) %>% 
  filter(genotype %in% all_of(interesting_genotypes))

geno_ranked_per_increasing_fabae = afabae_apisum_fecundity_contrasted_genotypes %>% 
  dplyr::filter(aphid == "A.fabae") %>%
  group_by(genotype) %>%
  summarise(larvae_total_number = median(larvae_total_number)) %>%
  arrange(larvae_total_number) %>%
  pull(genotype) %>%
  unique()

afabae_apisum_fecundity_contrasted_genotypes <- 
  afabae_apisum_fecundity_contrasted_genotypes %>% 
  mutate(aphid = factor(x = aphid, levels = c("A.pisum", "A.fabae")))

afabae_apisum_fecundity_contrasted_genotypes <- 
  afabae_apisum_fecundity_contrasted_genotypes %>% 
  mutate(genotype = factor(x = genotype, levels = geno_ranked_per_increasing_fabae))

p_interesting_genotypes = 
  afabae_apisum_fecundity_contrasted_genotypes %>%
  ggplot(aes(x = aphid, y = larvae_total_number, fill = aphid)) + 
  geom_boxplot() +
  geom_point() +
  facet_wrap(~ genotype, scales = "fixed") + 
  scale_fill_brewer(palette = 1, type = "qual") + 
  labs(x = "Aphid species", y = "Number of larvae obtained from one single L1 larvae") +
  theme_bw();p_interesting_genotypes 


ggsave(plot = p_interesting_genotypes, 
       filename = "02.analyses/2023_screen_specifics_251genotypes/03.comparison_apisum_vs_afabae/plot02_interesting_genotypes.pdf")
ggsave(plot = p_interesting_genotypes, 
       filename = "02.analyses/2023_screen_specifics_251genotypes/03.comparison_apisum_vs_afabae/plot02_interesting_genotypes.png")

