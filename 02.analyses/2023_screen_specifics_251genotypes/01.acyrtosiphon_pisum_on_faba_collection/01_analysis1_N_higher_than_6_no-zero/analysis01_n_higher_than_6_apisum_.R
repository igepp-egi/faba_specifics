suppressPackageStartupMessages(library("tidyverse"))
suppressPackageStartupMessages(library("lme4"))
suppressPackageStartupMessages(library("MASS"))
suppressPackageStartupMessages(library("emmeans"))
suppressPackageStartupMessages(library("patchwork"))

##############
## Data import
##############

###########################################################
### Genetic structure info (Jonathan, UMR Agroécologie/ECP)
##########################################################
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
apisum_fecundity <- read.csv("01.data/2023_screen_specifics_251genotypes/01.acyrthosiphon_pisum_P123amp_on_faba_collection/acyrthospipon_pisum_fecundity_on_faba_collection.csv", 
                            stringsAsFactors = F) %>% 
  mutate(type = gsub(pattern = "Equina", x = type, replacement = "Vicia faba var. equina")) %>% 
  mutate(type = gsub(pattern = "Minor", x = type, replacement = "Vicia faba var. minor")) %>% 
  mutate(type = gsub(pattern = "Major", x = type, replacement = "Vicia faba var. major")) %>% 
  mutate(type = gsub(pattern = "Paucijuga", x = type, replacement = "Vicia paucijuga")) %>% 
  filter(larvae_total_number != 0 & adult_alive == 1) %>% # Larvae died but initial L1 larvae developed into an adult
  inner_join(., igepp_ivip_genetic_group, by = "genotype") 

genotypes_nrep_higher_than_6 = 
  apisum_fecundity %>% 
  dplyr::select(ivip, larvae_total_number) %>% 
  group_by(ivip) %>% 
  tally() %>% 
  filter(n >= 6) %>% # keep only if N plants => 6 >> results 223 genotypes out of 248
  pull(ivip)

apisum_fecundity_n6 <- apisum_fecundity %>% filter(ivip %in% genotypes_nrep_higher_than_6)
rm(apisum_fecundity)

write.csv(x = apisum_fecundity_n6, 
          file = "02.analyses/2023_screen_specifics_251genotypes/01.acyrtosiphon_pisum_on_faba_collection/01_analysis1_N_higher_than_6_no-zero/apisum_fecundity_on_faba_n_higher_than_6_no_zeros.csv", 
          row.names = F, quote = F)

################
## Distributions
################

# Continuous normally distributed data of the total number of larvae ?
apisum_fecundity_n6 %>% 
  ggplot(., aes(x = larvae_total_number)) +
  geom_histogram(aes(y = ..density..), fill = "lightgrey", color = "black", bins = 30) +
  geom_density() +
  ggtitle("Distribution of the number of A. pisum larvae across all genotypes")

ggsave(filename = "02.analyses/2023_screen_specifics_251genotypes/01.acyrtosiphon_pisum_on_faba_collection/01_analysis1_N_higher_than_6_no-zero/plot01_distribution_larvae_number.pdf") 
ggsave(filename = "02.analyses/2023_screen_specifics_251genotypes/01.acyrtosiphon_pisum_on_faba_collection/01_analysis1_N_higher_than_6_no-zero/plot01_distribution_larvae_number.png") 

###################################
# Barplot A. pisum fecundity levels 
###################################

median_larvae_counts_apisum = 
  apisum_fecundity_n6 %>% 
  group_by(ivip) %>% 
  summarise(median_larvae = median(larvae_total_number)) %>% 
  arrange(median_larvae) 

order_apisum <- median_larvae_counts_apisum %>% 
  pull(ivip) %>% 
  unique()  

# 223 genotypes with N>6 values
p1 <- median_larvae_counts_apisum %>% 
  mutate(ivip = factor(ivip, levels = order_apisum)) %>% 
  ggplot(aes(x = ivip, y = median_larvae)) + 
  geom_bar(stat = "identity", fill="mediumblue", alpha = 0.3) +
  theme_bw() +
  theme(axis.text.x = element_blank(), 
        axis.title.x = element_text(size = 12), 
        axis.title.y = element_text(size = 12)) + 
  labs(x = "Vicia faba genotype", y = "Total number of larvae\n(from one single N1 larvae)") + 
  ggtitle("A. pisum fecundity on Vicia faba collection");p1

ggsave(plot = p1, width = 12, height = 8, 
       filename="02.analyses/2023_screen_specifics_251genotypes/01.acyrtosiphon_pisum_on_faba_collection/01_analysis1_N_higher_than_6_no-zero/plot02_barplot_apisum.pdf")

ggsave(plot = p1, width = 12, height = 8, units = "cm",
       filename="02.analyses/2023_screen_specifics_251genotypes/01.acyrtosiphon_pisum_on_faba_collection/01_analysis1_N_higher_than_6_no-zero/plot02_barplot_apisum.png")


#################
# 10 most extreme
#################
ten_most_resistant_genotypes_apisum <- 
  apisum_fecundity_n6 %>% 
  group_by(ivip) %>% 
  summarise(median_apisum = median(larvae_total_number)) %>% 
  arrange(median_apisum) %>% 
  head(n = 10) %>% 
  mutate(phenotype = "Ten most resistant")

ten_most_sensitive_genotypes_apisum <- 
  apisum_fecundity_n6 %>% 
  group_by(ivip) %>% 
  summarise(median_apisum = median(larvae_total_number)) %>% 
  arrange(median_apisum) %>% 
  tail(n = 10) %>%
  mutate(phenotype = "Ten most sensitive") 

# Number of biological replicate per genotype (to add to the graph)
n_per_genotype <- 
  apisum_fecundity_n6 %>% 
  dplyr::select(ivip, replicate) %>%
  group_by(ivip) %>% 
  tally() %>% 
  dplyr::rename("n_rep" = "n")

most_extreme_genotype_apisum <- 
  bind_rows(ten_most_resistant_genotypes_apisum, ten_most_sensitive_genotypes_apisum) %>% 
  inner_join(n_per_genotype, by = "ivip") %>%
  mutate(n_rep = paste("n=",n_rep, sep = ""))
  
p_most_contrasted_pisum <- 
  apisum_fecundity_n6 %>% 
  inner_join(., most_extreme_genotype_apisum) %>%
  mutate(ivip = factor(ivip, levels = order_apisum)) %>% 
  ggplot(., aes(x = ivip, 
                y = larvae_total_number, 
                fill = phenotype, 
                label = paste(n_rep))) +
  geom_boxplot() +
  geom_point(alpha=0.3) +
  scale_fill_brewer(type = 'qual', palette = 2) +
  ggtitle("Most contrasted genotypes to A. pisum\n(N>=6 and zeros removed)") +
  coord_flip() +
  labs(x = "Vicia faba genotype", y = "Total number of A. pisum larvae (from one L1)") +
  scale_y_continuous(limits = c(0,80)) +
  geom_text(y = 75);p_most_contrasted_pisum

ggsave(plot = p_most_contrasted_pisum, 
       filename = "02.analyses/2023_screen_specifics_251genotypes/01.acyrtosiphon_pisum_on_faba_collection/01_analysis1_N_higher_than_6_no-zero/plot03_apisum_most_contrasted_genotypes.pdf", 
       width = 10, height = 8)
ggsave(plot = p_most_contrasted_pisum, width = 16, height = 10, units = "cm",
       filename = "02.analyses/2023_screen_specifics_251genotypes/01.acyrtosiphon_pisum_on_faba_collection/01_analysis1_N_higher_than_6_no-zero/plot03_apisum_most_contrasted_genotypes.png")

# Second plot with short code
median_larvae_counts_apisum_by_shortcode = 
  apisum_fecundity_n6 %>% 
  group_by(short_code) %>% 
  summarise(median_larvae = median(larvae_total_number)) %>% 
  arrange(median_larvae) 

order_apisum_short_code <- median_larvae_counts_apisum_by_shortcode %>% 
  pull(short_code) %>% 
  unique()  

p_most_contrasted_pisum_shortcode <- 
  apisum_fecundity_n6 %>% 
  inner_join(., most_extreme_genotype_apisum) %>%
  mutate(short_code = factor(short_code, levels = order_apisum_short_code)) %>% 
  ggplot(., aes(x = short_code, 
                y = larvae_total_number, 
                fill = phenotype, 
                label = paste(n_rep))) +
  geom_boxplot() +
  geom_point(alpha=0.3) +
  scale_fill_brewer(type = 'qual', palette = 2) +
  ggtitle("Most contrasted genotypes to A. pisum\nShort code edition! \n(N>=6 and zeros removed)") +
  coord_flip() +
  labs(x = "Vicia faba genotype", y = "Total number of A. pisum larvae (from one L1)") +
  scale_y_continuous(limits = c(0,80)) +
  geom_text(y = 75);p_most_contrasted_pisum_shortcode

ggsave(plot = p_most_contrasted_pisum_shortcode, 
       filename = "02.analyses/2023_screen_specifics_251genotypes/01.acyrtosiphon_pisum_on_faba_collection/01_analysis1_N_higher_than_6_no-zero/plot04_apisum_most_contrasted_genotypes_shortcode.pdf", 
       width = 10, height = 8)

p_combined = p_most_contrasted_pisum + p_most_contrasted_pisum_shortcode
p_combined 

ggsave(plot = p_combined, 
       filename = "02.analyses/2023_screen_specifics_251genotypes/01.acyrtosiphon_pisum_on_faba_collection/01_analysis1_N_higher_than_6_no-zero/plot05_apisum_most_contrasted_genotypes_ivip_and_shortcode.pdf", 
       width = 16, height = 10)
#####################
## List of genotypes
#####################

apisum_fecundity_n6 %>% 
  inner_join(., most_extreme_genotype_apisum) %>% 
  dplyr::select(ivip, larvae_total_number, phenotype) %>% 
  group_by(ivip, phenotype) %>% 
  summarise(median_larvae = median(larvae_total_number, na.rm=T)) %>% 
  arrange(median_larvae) %>% 
  knitr::kable()

apisum_fecundity_n6 %>% 
  inner_join(., most_extreme_genotype_apisum) %>% 
  dplyr::select(genotype, larvae_total_number, phenotype) %>% 
  group_by(genotype, phenotype) %>% 
  summarise(median_larvae = median(larvae_total_number, na.rm=T)) %>% 
  arrange(median_larvae) %>% 
  write.csv(file = "02.analyses/2023_screen_specifics_251genotypes/01.acyrtosiphon_pisum_on_faba_collection/01_analysis1_N_higher_than_6_no-zero/apisum_most_contrasted_genotypes.csv", 
            quote = F, row.names = F)


########
## Stats
########
# reverse a vector of genotype names
reverse_order_apisum <- rev(order_apisum) %>% intersect(most_extreme_genotype_apisum$ivip)

# Data of fecundity of the 20 most contrasted genotypes
df <- 
  apisum_fecundity_n6 %>% 
  inner_join(., most_extreme_genotype_apisum) %>% 
  dplyr::select(ivip, larvae_total_number, replicate) 
df$ivip = factor(x = df$ivip, levels = reverse_order_apisum)

# Fit a null and genotype GLM negative binomial model
fit0 = glm.nb(formula = larvae_total_number ~ 1, data = df)
fit1 = glm.nb(formula = larvae_total_number ~ 1 + ivip, data = df)
anova(fit0, fit1)
# Likelihood ratio tests of Negative Binomial Models
# Response: larvae_total_number
# Model    theta Resid. df    2 x log-lik.   Test    df LR stat.      Pr(Chi)
# 1            1 5.235807       207       -1640.783                                   
# 2 1 + genotype 9.523030       188       -1550.237 1 vs 2    19 90.54536 2.655054e-11

# Fit a null (only replicate) and a genotype GLMM negative binomial model 
fit2 = lme4::glmer.nb(larvae_total_number ~ 1 + (1|ivip), data = df)
fit3 = lme4::glmer.nb(larvae_total_number ~ 1 + ivip + (1|replicate), data = df)
anova(fit2, fit3)
# Models:
# Data: df
# Models:
#   fit2: larvae_total_number ~ 1 + (1 | ivip)
# fit3: larvae_total_number ~ 1 + ivip + (1 | replicate)
# npar    AIC    BIC  logLik deviance  Chisq Df Pr(>Chisq)    
# fit2    3 1632.7 1642.8 -813.36   1626.7                         
# fit3   22 1616.5 1690.1 -786.26   1572.5 54.192 19  3.088e-05 ***

# Model comparison
AIC(fit0, fit1, fit2, fit3)

# Post-hoc test
pairs(emmeans(fit3, "ivip")) %>% 
  broom::tidy() %>% 
  # filter rows that contains string 'IVIP103069' > most sensitive genotype
  dplyr::filter(., grepl("IVIP103069",contrast))  %>% 
  dplyr::filter(adj.p.value < 0.05) %>% 
  knitr::kable()
 
# |term     |contrast        | null.value|  estimate| std.error|  df| statistic| adj.p.value|
# |:--------|:---------------|----------:|---------:|---------:|---:|---------:|-----------:|
# |genotype |JF2262 - JF2094 |          0| 0.6385175| 0.1624469| Inf|  3.930624|   0.0126658|
# |genotype |JF2262 - KF0160 |          0| 0.5900087| 0.1620285| Inf|  3.641388|   0.0360970|
# |genotype |JF2262 - JF2216 |          0| 0.5877373| 0.1617331| Inf|  3.633995|   0.0370176|
# |genotype |JF2262 - JF2274 |          0| 0.7872507| 0.1905080| Inf|  4.132377|   0.0056994|
# |genotype |JF2262 - KF0188 |          0| 0.8678889| 0.1749321| Inf|  4.961289|   0.0001268|
