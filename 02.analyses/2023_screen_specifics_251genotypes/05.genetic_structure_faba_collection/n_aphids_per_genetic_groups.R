##############################################
## N aphids per genetic group: any difference? 
#############################################

library(tidyverse)
library(patchwork)
suppressPackageStartupMessages(library(multcomp))

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

#######################
### Load fecundity data
#######################
apisum <- read.csv("01.data/2023_screen_specifics_251genotypes/01.acyrthosiphon_pisum_P123amp_on_faba_collection/apisum_fecundity_on_faba_n_higher_than_6_no_zeros.csv", 
                   header = TRUE, sep = ",") %>%
  inner_join(., igepp_ivip_genetic_group, by = "genotype") %>% 
  dplyr::select(genotype, cluster, larvae_total_number) %>% 
  mutate(cluster = as.factor(cluster))

afabae <- read.csv("01.data/2023_screen_specifics_251genotypes/02.aphis_fabae_on_faba_collection/aphis_fabae_fecundity_on_faba_collection_n_higher_6.csv", 
                   header = TRUE, sep = ",") %>%
  filter(larvae_total_number > 0) %>%
  inner_join(., igepp_ivip_genetic_group, by = "genotype") %>%
  dplyr::select(genotype, cluster, larvae_total_number) %>%
  mutate(cluster = as.factor(cluster))

#########################################
### Plot fecundity data per genetic group
#########################################

jitter_width = 0.05
alpha_point = 0.1

p1 = 
  apisum %>%
  ggplot(aes(x = cluster, y = larvae_total_number, fill = cluster)) +
  geom_violin() +
  geom_jitter(width = jitter_width, size = 1, alpha = alpha_point) +
  labs(title = "A. pisum fecundity = f(genetic group)", x = "Genetic group", y = "Number of A. pisum offspring") +
  theme_minimal() + 
  scale_fill_brewer(palette = "Set1") +
  stat_summary(fun = median, geom = "point", shape = 23, size = 3, fill = "white") + 
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 14)) +
  theme(legend.position = "none")

p2 = 
  afabae %>%
  ggplot(aes(x = cluster, y = larvae_total_number, fill = cluster)) +
  geom_violin() +
  geom_jitter(width = jitter_width, size = 1, alpha = alpha_point) +
  labs(title = "A. fabae fecundity = f(genetic group)", x = "Genetic group", y = "Number of A. fabae offspring") +
  theme_minimal() + 
  scale_fill_brewer(palette = "Set2") +
  # add median as a point to violinplot
  stat_summary(fun = median, geom = "point", shape = 23, size = 3, fill = "white") +
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 14)) + 
  # remove legend
  theme(legend.position = "none")
  

#####################################################################################
# Calculate number of values per genetic group for each aphid species and add to plot
#####################################################################################

p1 <- p1 + geom_text(data = apisum %>% group_by(cluster) %>% summarise(n = n()), aes(x = cluster, y = 70, label = paste("n=",n)), vjust = -0.5) 

p2 <- p2 + geom_text(data = afabae %>% group_by(cluster) %>% summarise(n = n()), aes(x = cluster, y = 100, label = paste("n=",n)), vjust = -0.5)

p1 + p2

ggsave("02.analyses/2023_screen_specifics_251genotypes/05.genetic_structure_faba_collection/n_aphids_per_genetic_group.pdf", width = 12, height = 7)


#####################################################################
### Compute a one-way ANOVA and Tukey HSD test for both aphid species
#####################################################################
apisum_aov <- aov(larvae_total_number ~ cluster, data = apisum)
summary(apisum_aov)
tukey_res_apisum = TukeyHSD(apisum_aov) %>% broom::tidy() %>% 
  # add a column with signif if padj < 0.05
  mutate(signif = ifelse(adj.p.value < 0.05, "yes", "no"))
tukey_res_apisum


# same for Aphis fabae
afabae_aov <- aov(larvae_total_number ~ cluster, data = afabae)
summary(afabae_aov)
tukey_res_afabae = TukeyHSD(afabae_aov) %>% broom::tidy() %>% 
  # add a column with signif if padj < 0.05
  mutate(signif = ifelse(adj.p.value < 0.05, "yes", "no"))
tukey_res_afabae
