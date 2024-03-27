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

### Fecundity data
apisum_fecundity <- read.csv("01.data/01.acyrthosiphon_pisum_P123amp_on_faba_collection/acyrthospipon_pisum_fecundity_on_faba_collection.csv", 
                            stringsAsFactors = F) %>% 
  mutate(type = gsub(pattern = "Equina", x = type, replacement = "Vicia faba var. equina")) %>% 
  mutate(type = gsub(pattern = "Minor", x = type, replacement = "Vicia faba var. minor")) %>% 
  mutate(type = gsub(pattern = "Major", x = type, replacement = "Vicia faba var. major")) %>% 
  mutate(type = gsub(pattern = "Paucijuga", x = type, replacement = "Vicia paucijuga")) %>% 
  filter(larvae_total_number != 0) # Larvae probably died

##############################
## Faba collection description
##############################

# Type of faba
apisum_fecundity %>% 
  select(genotype, type) %>% 
  distinct() %>% 
  ggplot(., aes(x = type, fill = type)) +
  geom_bar(stat = "count", col = "black") +
  scale_fill_brewer(type = 'qual', palette = 3) +
  geom_text(stat='count', aes(label=..count..), hjust=-1) +
  scale_y_continuous(limits = c(0,125)) +
  guides(fill = "none") +
  labs(x = "Vicia faba classification", y = "Number of occurences in the collection") +
  coord_flip() +
  theme(axis.text = element_text(size = 12))

ggsave(filename = "02.analyses/01.acyrtosiphon_pisum_on_faba_collection/plot01_nb_occurences_per_faba_type.pdf")  

# Number of biological replicate per genotype
n_per_genotype <- 
  apisum_fecundity %>% 
  drop_na() %>% 
  filter(adult_alive != 0) %>%       # Keep only data with adult alive (from first L1 larvae deposited)
  select(genotype, replicate) %>%
  group_by(genotype) %>% 
  tally() %>% 
  dplyr::rename("n_rep" = "n")

apisum_n_rep_per_genotype <- 
  ggplot(n_per_genotype, aes(x = n_rep)) +
  geom_bar() + 
  scale_color_brewer(type = 'seq') +
  labs(x = "Number of biological replicates (plants)", y = "Genotype occurences") +
  ggtitle("A. pisum screen") +
  scale_x_continuous(breaks = seq(0,14,1))
apisum_n_rep_per_genotype

ggsave(plot = apisum_n_rep_per_genotype,
       filename = "02.analyses/01.acyrtosiphon_pisum_on_faba_collection/plot02_nb_rep_per_genotype.pdf")  

###################################################
## Stat models: genotype effect on number of larvae
###################################################

# Continuous normally distributed data of the total number of larvae ?
apisum_fecundity %>% 
  drop_na() %>% 
  filter(adult_alive != 0) %>% 
  ggplot(., aes(x = larvae_total_number)) +
  geom_histogram(aes(y = ..density..), fill = "lightgrey", color = "black", bins = 30) +
  geom_density() +
  ggtitle("Distribution of the number of A. pisum larvae across all genotypes")

ggsave(filename = "02.analyses/01.acyrtosiphon_pisum_on_faba_collection/plot03_distribution_larvae_number.pdf")  

# Model 01: only genotype
fit1 <- 
  apisum_fecundity %>% 
  drop_na() %>% 
  filter(adult_alive != 0) %>% 
  mutate(genotype = as.factor(genotype)) %>% 
  with(., aov(formula = larvae_total_number ~ genotype, data = .))
summary(fit1)
# Df Sum Sq Mean Sq F value   Pr(>F)    
# genotype     240  34228  142.62   1.539 9.97e-07 ***
# Residuals   2034 188431   92.64                     
# ---

# Model 02: genotype and session
fit2 <- 
  apisum_fecundity %>% 
  drop_na() %>% 
  filter(adult_alive != 0) %>% 
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
  apisum_fecundity %>% 
  drop_na() %>% 
  filter(adult_alive != 0) %>% 
  mutate(genotype = as.factor(genotype)) %>% 
  with(., aov(formula = larvae_total_number ~ genotype + session + genotype:session, data = .))
summary(fit3)
#                     Df Sum Sq Mean Sq F value   Pr(>F)    
#   genotype          240  34228     143   2.432  < 2e-16 ***
#   session             6  26828    4471  76.244  < 2e-16 ***
#   genotype:session 1129 108881      96   1.644 4.48e-15 ***
#   Residuals         899  52721      59                     
# ---
  
# AIC criterion: smaller the better
aic_values <- data.frame(fit = paste("fit", seq(1,3,1)),
                   aic = c(
                     AIC(fit1), 
                     AIC(fit2), 
                     AIC(fit3))
                   )
ggplot(aic_values, aes(x = fit, y = aic)) +
  geom_bar(stat = "identity")

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
stat_res_table_pisum <- apisum_fecundity %>% 
  select(genotype, replicate, larvae_total_number) %>% 
  group_by(genotype) %>% 
  summarise(median_larvae = median(larvae_total_number, na.rm = T), 
            sd_larvae = sd(larvae_total_number, na.rm = T)) %>% 
  ungroup() %>% 
  inner_join(., cld_df) %>% 
  arrange(median_larvae) %>% 
  mutate(genotype = factor(genotype, levels = genotype))

write.csv(x = stat_res_table_pisum, 
          file = '02.analyses/01.acyrtosiphon_pisum_on_faba_collection/table01_stat_faba_genotypes.csv', 
          row.names = F)

#############
# Geom ribbon
#############
ggplot(stat_res_table_pisum, aes(x = genotype, y = median_larvae)) + 
  geom_point() + 
  geom_segment(aes(x=genotype, xend=genotype, y=0, yend=median_larvae))

order_apisum <- apisum_fecundity %>% 
  group_by(genotype) %>% 
  summarise(median_apisum = median(larvae_total_number)) %>% 
  arrange(median_apisum) %>% 
  pull(genotype) %>% 
  unique()

apisum_fecundity %>% 
  drop_na() %>% 
  ggplot(aes(x = genotype, y = larvae_total_number)) + 
  stat_lineribbon(.width = c(0.50, 0.90, 0.95), group=1) +
  scale_fill_brewer() + 
  coord_flip() +
  theme(axis.text.y = element_text(size = 4)) +
  labs(x = "Vicia faba genotype", y = "Number of A. pisum larvae obtained from one L1 (Median -/+ CI)")
ggplot(filename="02.analyses/01.acyrtosiphon_pisum_on_faba_collection/plot03_apisum_ribbonplot.pdf")

#################
# 10 most extreme
#################
ten_most_resistant_genotypes_apisum <- 
  apisum_fecundity %>% 
  group_by(genotype) %>% 
  summarise(median_apisum = median(larvae_total_number)) %>% 
  arrange(median_apisum) %>% 
  head(n = 10) %>% 
  mutate(phenotype = "Ten most resistant")

ten_most_sensitive_genotypes_apisum <- 
  apisum_fecundity %>% 
  group_by(genotype) %>% 
  summarise(median_apisum = median(larvae_total_number)) %>% 
  arrange(median_apisum) %>% 
  tail(n = 10) %>%
  mutate(phenotype = "Ten most sensitive") 

most_extreme_genotype_apisum <- bind_rows(ten_most_resistant_genotypes_apisum, ten_most_sensitive_genotypes_apisum) %>% 
  inner_join(n_per_genotype, by = "genotype") %>%
  mutate(n_rep = paste("n=",n_rep, sep = "")) %>% 
  inner_join(igepp_ivip_genetic_group, by = "genotype") %>% 
  select(genotype, phenotype, n_rep, clusterName, cluster) %>% # add n replicates + genetic cluster info
  mutate(cluster = paste("genetic-group",cluster,sep = ""))
  

apisum_fecundity %>% 
  inner_join(., most_extreme_genotype_apisum) %>%
  mutate(genotype = factor(genotype, levels = order_apisum)) %>% 
  ggplot(., aes(x = genotype, y = larvae_total_number, fill = phenotype, label = paste(n_rep, cluster,sep = ", "))) +
  geom_boxplot() +
  geom_point(alpha=0.3) +
  scale_fill_brewer(type = 'qual', palette = 2) +
  ggtitle("Most contrasted genotypes\n (A. pisum screen)") +
  coord_flip() +
  labs(x = "Vicia faba genotype", y = "Total number of A. pisum larvae (from one L1)") +
  scale_y_continuous(limits = c(0,80)) +
  geom_text(y = 75)

ggsave(filename = "02.analyses/01.acyrtosiphon_pisum_on_faba_collection/plot03_apisum_10_most_contrasted_genotypes.pdf")
