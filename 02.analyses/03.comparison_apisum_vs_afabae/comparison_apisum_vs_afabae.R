library("ggpubr")
library("tidyverse")
library("patchwork")
library("ggrepel")

###################################################
## Execute single analysis of A. pisum and A. fabae
## Analysis 01: without zeros and N>=6
##################################################

source("02.analyses/01.acyrtosiphon_pisum_on_faba_collection/01_analysis1_N_higher_than_6_no-zero/analysis01_n_higher_than_6_apisum_.R")
source("02.analyses/02.aphis_fabae_on_faba_collection/01_analysis1_N_higher_than_6_no-zeros/analysis01_n_higher_than_6_aphis_fabae.R")

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
  select(genotype, A_fabae, A_pisum, cluster)

genotypes_res_afabae <- m %>% filter(phenotype == "Ten most resistant") %>% pull(genotype)
genotypes_res_pisum <- ten_most_resistant_genotypes_apisum %>% filter(phenotype == "Ten most resistant") %>% pull(genotype)
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
  scale_y_continuous(limits = c(0,60))
p_paired_continuous

ggsave(plot = p_paired_continuous, 
       filename = "02.analyses/03.comparison_apisum_vs_afabae/plot01_paired_plot_apisum_vs_afabae.pdf", 
       width = 10, height = 8)

### Scatterplot
points_to_label <- 
  res_pisum_and_afabae %>% 
  filter(genotype %in% all_genotypes)

p_scatterplot_continuous <- 
  res_pisum_and_afabae %>% 
  ggplot(., aes(x = A_fabae, y = A_pisum, fill = factor(cluster))) +
  geom_point() +
  labs(x = "Number of larvae from one L1 larvae of A. fabae") + 
  labs(y = "Number of larvae from one L1 larvae of A. pisum") +
  geom_label_repel(data = points_to_label, 
                   mapping = aes(x = A_fabae, y = A_pisum, label = genotype)) 
p_scatterplot_continuous

ggsave(plot = p_scatterplot_continuous, 
       filename = "02.analyses/03.comparison_apisum_vs_afabae/plot02_scatterplot_apisum_vs_afabae.pdf", 
       width = 10, height = 8)
ggsave(plot = p_scatterplot_continuous, 
       filename = "02.analyses/03.comparison_apisum_vs_afabae/plot02_scatterplot_apisum_vs_afabae.png")

### All plots
p_continuous <- p_paired_continuous + p_scatterplot_continuous
ggsave(plot = p_continuous, 
       filename = "02.analyses/03.comparison_apisum_vs_afabae/comparison_apisum_and_afabae.pdf", 
       width = 12, height = 8)

#######
# Ranks
#######

stat_res_table_aphis_fabae_ranks <-  
  stat_res_table_aphis_fabae %>% 
  mutate(rank_afabae = dense_rank(median_afabae_larvae))

stat_res_table_pisum_ranks <- 
  stat_res_table_pisum %>% 
  mutate(rank_apisum = dense_rank(median_larvae))

d_ranks <- 
  inner_join(stat_res_table_aphis_fabae_ranks, 
             stat_res_table_pisum_ranks, 
             by = "genotype") %>% 
  select(genotype, rank_afabae, rank_apisum)

p_paired_ranks <- ggpaired(d_ranks, 
                                cond1 = "rank_apisum", 
                                cond2 = "rank_afabae",
                                line.color = "azure3", 
                                line.size = 0.4,
                           label = "genotype", 
                           label.select = c(
                             "JF2118","JF2240","JF2188", # most R to A. fabae
                             "KF0188", "JF2274", "JF2216"))  +
  labs(x = "Aphid species", 
       y = "Rank based on L1 fecundity (low rank, low fecundity)")
p_paired_ranks
  
p_scatterplot_ranks <- 
  d_ranks %>% 
  ggplot(., aes(x = rank_afabae, y = rank_apisum)) +
  geom_point() +
  labs(x = "Rank based on A. fabae median fecundity\n(low rank, low fecundity, ties allowed)") +
  labs(y = "Rank based on A. pisum median fecundity\n(low rank, low fecundity, ties allowed)") 

p_ranks <- p_paired_ranks + p_scatterplot_ranks

###########
# All plots
###########

p_continuous / p_ranks
ggsave(filename = "02.analyses/03.comparison_apisum_vs_afabae/comparison_apisum_and_afabae.pdf", 
       width = 14, height = 10)
