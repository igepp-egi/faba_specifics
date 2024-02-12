library("ggpubr")
library("tidyverse")
library("patchwork")
library("ggrepel")

###################################################
# Execute single analysis of A. pisum and A. fabae
##################################################

source("02.analyses/01.acyrtosiphon_pisum_on_faba_collection/apisum_on_faba_analysis.R")
source("02.analyses/02.aphis_fabae_on_faba_collection/aphis_fabae_analysis.R")

##################
# N rep / genotype
##################
apisum_n_rep_per_genotype + aphis_fabae_n_rep_per_genotype
ggsave(filename = "02.analyses/03.comparison_apisum_vs_afabae/n_rep_per_geno.pdf")

#################################
# Paired plots (absolute numbers)
#################################

d <- inner_join(stat_res_table_aphis_fabae, stat_res_table_pisum, by = "genotype") %>% 
  select(genotype, median_afabae_larvae, median_larvae) %>% 
  rename("A_pisum" = "median_larvae") %>% 
  rename("A_fabae" = "median_afabae_larvae")

p_paired_continuous <- ggpaired(d, 
                                cond1 = "A_pisum", 
                                cond2 = "A_fabae",
                                line.color = "lightgray", 
                                line.size = 0.4, 
                                label = "genotype", 
                                label.select = c(
                                  "JF2118","JF2240","JF2188",  # most R to A. fabae
                                  "KF0188", "JF2274", "JF2216" # Most R to A. pisum
                                  )) +
  labs(x = "Aphid species", y = "Number of larvae obtained from one L1 larvae") +
  scale_y_continuous(limits = c(0,60)) 
p_paired_continuous

scatterplot_continous_to_label <- 
  d %>% 
  filter(genotype %in% c("JF2118","JF2240","JF2188", "KF0188", "JF2274", "JF2216"))

p_scatterplot_continuous <- 
  d %>% 
  ggplot(., aes(x = A_fabae, y = A_pisum)) +
  geom_point() +
  labs(x = "Number of larvae from one L1 larvae of A. fabae") + 
  labs(y = "Number of larvae from one L1 larvae of A. pisum") +
  geom_label_repel(data = scatterplot_continous_to_label, 
                   mapping = aes(x = A_fabae, y = A_pisum, label = genotype))
p_scatterplot_continuous

p_continuous <- p_paired_continuous + p_scatterplot_continuous

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
