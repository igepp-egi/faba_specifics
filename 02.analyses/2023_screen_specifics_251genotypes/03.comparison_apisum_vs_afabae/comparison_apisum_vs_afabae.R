library("ggpubr")
library("tidyverse")
library("patchwork")
library("ggrepel")

###################################################
## Execute single analysis of A. pisum and A. fabae
## Analysis 01: without zeros and N>=6
##################################################

source("02.analyses/2023_screen_specifics_251genotypes/01.acyrtosiphon_pisum_on_faba_collection/01_analysis1_N_higher_than_6_no-zero/analysis01_n_higher_than_6_apisum_.R")
source("02.analyses/2023_screen_specifics_251genotypes/02.aphis_fabae_on_faba_collection/01_analysis1_N_higher_than_6_no-zeros/analysis01_n_higher_than_6_aphis_fabae.R")

#################################
# Paired plots (absolute numbers)
#################################

# 223 genotypes
res_pisum <- 
  apisum_fecundity_n6 %>% 
  group_by(genotype) %>% 
  summarise(median_apisum = median(larvae_total_number)) %>% 
  arrange(median_apisum) 

# 193 genotypes
res_afabae <- 
  afabae_fecundity_n6 %>% 
  group_by(genotype) %>% 
  summarise(median_afabae = median(larvae_total_number)) %>% 
  arrange(median_afabae)

# 181 genotypes communs
res_pisum_and_afabae <- inner_join(res_pisum, res_afabae, by = "genotype") %>% 
  rename("A_pisum" = "median_apisum") %>% 
  rename("A_fabae" = "median_afabae") %>% 
  inner_join(., igepp_ivip_genetic_group) %>% 
  select(genotype, A_fabae, A_pisum, clusterName)

genotypes_res_afabae <- most_extreme_genotype_afabae %>% filter(phenotype == "Ten most resistant") %>% pull(genotype)
genotypes_res_pisum <- most_extreme_genotype_apisum %>% filter(phenotype == "Ten most resistant") %>% pull(genotype)
all_genotypes = c(genotypes_res_afabae, genotypes_res_pisum)

### Paired plot
p_paired_continuous <- ggpaired(res_pisum_and_afabae, 
                                cond1 = "A_pisum", 
                                cond2 = "A_fabae",
                                line.color = "lightgray", 
                                line.size = 0.4,
                                label = "genotype",
                                label.select = c(genotypes_res_afabae, genotypes_res_pisum)) + 
  labs(x = "Aphid species", y = "Number of larvae obtained from one L1 larvae") +
  scale_y_continuous(limits = c(0,60)); p_paired_continuous

ggsave(plot = p_paired_continuous, 
       filename = "02.analyses/2023_screen_specifics_251genotypes/03.comparison_apisum_vs_afabae/plot01_paired_plot_apisum_vs_afabae.pdf", 
       width = 10, height = 8)

### Scatterplot
points_to_label <- 
  res_pisum_and_afabae %>% 
  filter(genotype %in% all_genotypes) 

p_scatterplot_continuous <- 
  res_pisum_and_afabae %>% 
  ggplot(., aes(x = A_fabae, y = A_pisum, label = genotype)) +
  geom_point() +
  labs(x = "Number of larvae from one L1 larvae of A. fabae") + 
  labs(y = "Number of larvae from one L1 larvae of A. pisum") +
  geom_label_repel(data = points_to_label, max.overlaps = 10) +
  scale_fill_brewer(type = "qual", palette = 2) +
  theme_bw();p_scatterplot_continuous

ggsave(plot = p_scatterplot_continuous, 
       filename = "02.analyses/2023_screen_specifics_251genotypes/03.comparison_apisum_vs_afabae/plot02_scatterplot_apisum_vs_afabae.pdf", 
       width = 7, height = 7)
ggsave(plot = p_scatterplot_continuous, width = 10, height = 10, units = "cm",
       filename = "02.analyses/2023_screen_specifics_251genotypes/03.comparison_apisum_vs_afabae/plot02_scatterplot_apisum_vs_afabae.png")

p_scatterplot_continuous + ggrepel::geom_label_repel(max.overlaps = Inf) 


#######
# Ranks
#######

ranks_afabae <-  
  res_afabae %>% 
  mutate(rank_afabae = dense_rank(median_afabae))

ranks_apisum <- 
  res_pisum %>% 
  mutate(rank_apisum = dense_rank(median_apisum))

ranks_both_aphids <- 
  inner_join(ranks_afabae, 
             ranks_apisum, 
             by = "genotype") %>% 
  select(genotype, rank_afabae, rank_apisum)

p_paired_ranks <- ggpaired(ranks_both_aphids, 
                                cond1 = "rank_apisum", 
                                cond2 = "rank_afabae",
                                line.color = "azure3", 
                                line.size = 0.4);p_paired_ranks

 #                          label = "genotype", 
#                           label.select = c(
#                             "JF2118","JF2240","JF2188", # most R to A. fabae
#                             "KF0188", "JF2274", "JF2216"))  +
#  + labs(x = "Aphid species", 
#      y = "Rank based on L1 fecundity (low rank, low fecundity)")

  
p_scatterplot_ranks <- 
  ranks_both_aphids %>% 
  ggplot(., aes(x = rank_afabae, y = rank_apisum, label = genotype)) +
  geom_point() +
  labs(x = "Rank based on A. fabae median fecundity\n(low rank, low fecundity, ties allowed)") +
  labs(y = "Rank based on A. pisum median fecundity\n(low rank, low fecundity, ties allowed)") +
  theme_bw()
p_scatterplot_ranks

# with labels
p_scatterplot_ranks + ggrepel::geom_label_repel(max.overlaps = Inf)
ggsave(filename = "02.analyses/2023_screen_specifics_251genotypes/03.comparison_apisum_vs_afabae/plot03_scatterplot_ranks_all_labeled.pdf")

###########################################
# Closer look at some interesting genotypes
###########################################

interesting_genotypes = c("JF2094", "JF2274", "JF2256", "KF2219")

afabae_fecundity_n6_subset = afabae_fecundity_n6 %>% select(genotype, aphid_species, larvae_total_number) %>% rename("aphid" = "aphid_species")
apisum_fecundity_n6_subset = apisum_fecundity_n6 %>% select(genotype, aphid, larvae_total_number) 
afabae_apisum_fecundity = bind_rows(afabae_fecundity_n6_subset, apisum_fecundity_n6_subset)


p_interesting_genotypes = 
  afabae_apisum_fecundity %>% 
  filter(genotype %in% all_of(interesting_genotypes)) %>% 
  ggplot(aes(x = aphid, y = larvae_total_number, fill = aphid)) + 
  geom_boxplot() +
  geom_point() +
  facet_wrap(~ genotype, scales = "fixed") + 
  scale_fill_brewer(palette = 1, type = "qual") + 
  labs(x = "Aphid species", y = "Number of larvae obtained from one single L1 larvae") +
  theme_bw()
p_interesting_genotypes 


ggsave(plot = p_interesting_genotypes, 
       filename = "02.analyses/2023_screen_specifics_251genotypes/03.comparison_apisum_vs_afabae/plot04_interesting_genotypes.pdf")
ggsave(plot = p_interesting_genotypes, 
       filename = "02.analyses/2023_screen_specifics_251genotypes/03.comparison_apisum_vs_afabae/plot04_interesting_genotypes.png")


ranks_both_aphids %>% 
  filter(genotype %in% all_of(interesting_genotypes)) %>% 
  pivot_longer(- genotype, names_to = "species") %>% 
  mutate(species = gsub(replacement = "", x = species, perl = T, pattern="rank_")) %>% 
  ggplot(aes(x = species, y = value, fill = species)) + 
  geom_boxplot() +
  geom_point() +
  facet_wrap(~ genotype, scales = "fixed") + 
  scale_fill_brewer(palette = 1, type = "qual")



# JF2094
# Resistant to both aphids
res_pisum_and_afabae %>% filter(genotype == "JF2094")

# JF2274
# Resistant to A. pisum and sensitive to A. fabae

# JF2109 

# JF2256 and KF2219
