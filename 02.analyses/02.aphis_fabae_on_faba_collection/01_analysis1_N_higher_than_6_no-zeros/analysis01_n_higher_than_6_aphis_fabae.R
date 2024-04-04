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
afabae_fecundity <- read.csv("01.data/02.aphis_fabae_on_faba_collection/aphis_fabae_fecundity_on_faba_collection_2024-04-02.csv", 
                            stringsAsFactors = F) %>% 
  mutate(type = gsub(pattern = "Equina", x = type, replacement = "Vicia faba var. equina")) %>% 
  mutate(type = gsub(pattern = "Minor", x = type, replacement = "Vicia faba var. minor")) %>% 
  mutate(type = gsub(pattern = "Major", x = type, replacement = "Vicia faba var. major")) %>% 
  mutate(type = gsub(pattern = "Paucijuga", x = type, replacement = "Vicia paucijuga")) %>% 
  filter(larvae_total_number != 0) # Larvae probably died

genotypes_nrep_higher_than_6 = 
  afabae_fecundity %>% 
  select(genotype, larvae_total_number) %>% 
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

ggsave(filename = "02.analyses/02.aphis_fabae_on_faba_collection/01_analysis1_N_higher_than_6_no-zeros/plot01_distribution_larvae_number.pdf") 
ggsave(filename = "02.analyses/02.aphis_fabae_on_faba_collection/01_analysis1_N_higher_than_6_no-zeros/plot01_distribution_larvae_number.png") 

# Model 01: only genotype
fit1 <- 
  afabae_fecundity_n6 %>% 
  mutate(genotype = as.factor(genotype)) %>% 
  with(., aov(formula = larvae_total_number ~ genotype, data = .))
summary(fit1)
#             Df    Sum Sq  Mean Sq     F value   Pr(>F)    
# genotype    192   41430   215.8       1.535     1.31e-05 ***
# Residuals   1512  212579  140.6              
# ---

# Model 02: genotype and session
fit2 <- 
  afabae_fecundity_n6 %>% 
  mutate(genotype = as.factor(genotype)) %>% 
  with(., aov(formula = larvae_total_number ~ genotype + session, data = .))
summary(fit2)
#              Df Sum Sq Mean Sq F value  Pr(>F)    
# genotype     240  34228     143    1.79 3.3e-11 ***
# session        6  26828    4471   56.11 < 2e-16 ***
# Residuals   2028 161603      80                    
# ---

# Model 03: genotype and session and their interaction
fit3 <- 
  afabae_fecundity_n6 %>% 
  mutate(genotype = as.factor(genotype)) %>% 
  with(., aov(formula = larvae_total_number ~ genotype + session + genotype:session, data = .))
summary(fit3)
#                     Df Sum Sq Mean Sq F value   Pr(>F)    
#   genotype          240  34228     143   2.432  < 2e-16 ***
#   session             6  26828    4471  76.244  < 2e-16 ***
#   genotype:session 1129 108881      96   1.644 4.48e-15 ***
#   Residuals         899  52721      59                     
# ---
  
# Posthoc tests
tukey_res_fit1 <- TukeyHSD(fit1) 

# make the letter magic !!
cld = multcompLetters4(fit1, tukey_res_fit1, threshold = 0.01)
cld_df = as.data.frame.list(cld$genotype) %>% 
  as.data.frame() %>% 
  rownames_to_column("genotype") %>% 
  select(genotype, Letters)

# Add letters to data
# Save as csv
stat_res_table_afabae <- afabae_fecundity_n6 %>% 
  select(genotype, replicate, larvae_total_number) %>% 
  group_by(genotype) %>% 
  summarise(median_larvae = median(larvae_total_number, na.rm = T), 
            sd_larvae = sd(larvae_total_number, na.rm = T)) %>% 
  ungroup() %>% 
  inner_join(., cld_df) %>% 
  arrange(median_larvae) %>% 
  mutate(genotype = factor(genotype, levels = genotype))

write.csv(x = stat_res_table_afabae, 
          file = '02.analyses/02.aphis_fabae_on_faba_collection/01_analysis1_N_higher_than_6_no-zeros/table01_stat_afabae_faba_genotypes.csv', 
          row.names = F)

#############
# Geom ribbon
#############
ggplot(stat_res_table_afabae, aes(x = genotype, y = median_larvae)) + 
  geom_point() + 
  geom_segment(aes(x=genotype, xend=genotype, y=0, yend=median_larvae))

order_afabae <- afabae_fecundity_n6 %>% 
  group_by(genotype) %>% 
  summarise(median_afabae = median(larvae_total_number)) %>% 
  arrange(median_afabae) %>% 
  pull(genotype) %>% 
  unique()

ribbon_plot_afabae <- afabae_fecundity_n6 %>% 
  mutate(genotype = factor(genotype, levels = order_afabae)) %>% 
  drop_na() %>% 
  ggplot(aes(x = genotype, y = larvae_total_number)) + 
  stat_lineribbon(.width = c(0.50, 0.90, 0.95), group=1) +
  scale_fill_brewer() + 
  coord_flip() +
  theme(axis.text.y = element_text(size = 4)) +
  labs(x = "Vicia faba genotype", y = "Number of A. fabae larvae obtained from one L1 (Median -/+ CI)")

ribbon_plot_afabae

ggplot(plot = ribbon_plot_afabae, filename="02.analyses/02.aphis_fabae_on_faba_collection/01_analysis1_N_higher_than_6_no-zeros/plot02_afabae_ribbonplot.pdf")
ggplot(plot = ribbon_plot_afabae, filename="02.analyses/02.aphis_fabae_on_faba_collection/01_analysis1_N_higher_than_6_no-zeros/plot02_afabae_ribbonplot.png")

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
  select(genotype, replicate) %>%
  group_by(genotype) %>% 
  tally() %>% 
  dplyr::rename("n_rep" = "n")

most_extreme_genotype_afabae <- bind_rows(ten_most_resistant_genotypes_afabae, ten_most_sensitive_genotypes_afabae) %>% 
  inner_join(n_per_genotype, by = "genotype") %>%
  mutate(n_rep = paste("n=",n_rep, sep = "")) %>% 
  inner_join(igepp_ivip_genetic_group, by = "genotype") %>% 
  select(genotype, phenotype, n_rep, clusterName, cluster) %>% # add n replicates + genetic cluster info
  mutate(cluster = paste("group",cluster,sep = ""))
  

most_contrasted_fabae <- 
  afabae_fecundity_n6 %>% 
  inner_join(., most_extreme_genotype_afabae) %>%
  mutate(genotype = factor(genotype, levels = order_afabae)) %>% 
  ggplot(., aes(x = genotype, y = larvae_total_number, fill = phenotype, label = paste(n_rep, cluster,sep = ", "))) +
  geom_boxplot() +
  geom_point(alpha=0.2) +
  scale_fill_brewer(type = 'qual', palette = 1) +
  ggtitle("Most contrasted genotypes to Aphis fabae\n(N>=6 and zeros removed)") +
  coord_flip() +
  labs(x = "Vicia faba genotype", y = "Total number of Aphis fabae larvae (from one L1)") +
  scale_y_continuous(limits = c(0,100)) +
  geom_text(y = 80, hjust=0) + 
  theme(axis.text = element_text(size = 10))
most_contrasted_fabae

ggsave(plot = most_contrasted_fabae, 
       filename = "02.analyses/02.aphis_fabae_on_faba_collection/01_analysis1_N_higher_than_6_no-zeros/plot03_afabae_most_contrasted_genotypes.pdf", 
       width = 10, height = 8)
ggsave(plot = most_contrasted_fabae, 
       filename = "02.analyses/02.aphis_fabae_on_faba_collection/01_analysis1_N_higher_than_6_no-zeros/plot03_afabae_most_contrasted_genotypes.png")

#####################
## List of genotypes
#####################

afabae_fecundity_n6 %>% 
  inner_join(., most_extreme_genotype_afabae) %>% 
  select(genotype, larvae_total_number, phenotype) %>% 
  group_by(genotype, phenotype) %>% 
  summarise(median_larvae = median(larvae_total_number, na.rm=T)) %>% 
  arrange(median_larvae) %>% 
  knitr::kable()

afabae_fecundity_n6 %>% 
  inner_join(., most_extreme_genotype_afabae) %>% 
  select(genotype, larvae_total_number, phenotype) %>% 
  group_by(genotype, phenotype) %>% 
  summarise(median_larvae = median(larvae_total_number, na.rm=T)) %>% 
  arrange(median_larvae) %>% 
  write.csv(file = "02.analyses/02.aphis_fabae_on_faba_collection/01_analysis1_N_higher_than_6_no-zeros/most_extreme_genotypes_aphis_fabae.csv", quote = F, row.names = F)

