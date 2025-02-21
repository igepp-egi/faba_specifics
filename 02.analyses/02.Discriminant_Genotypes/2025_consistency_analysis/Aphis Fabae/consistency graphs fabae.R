#fix me: faire une annalyse bloc et BBCH avec P2xp009 et 10

# Packages ----
suppressPackageStartupMessages(library(corrplot))
suppressPackageStartupMessages(library(Hmisc))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(agricolae))
suppressPackageStartupMessages(library(emmeans))

#################################################################################
#################################################################################


# Import all fecundity measurements from 2022-2023 (N=251 genotypes), P2-xp007, P2-xp009
source("02.analyses/02.Discriminant_Genotypes/2025_consistency_analysis/Aphis Fabae/preparation data fabae.R")


####### ploting median consistency ##########################################
# plot consistency of all genotypes from P2XP007----
# Plotting with ggplot2 
ggplot(df_long2, aes(x = xp, y = median_larvae_number, group = ivip, color = ivip)) +
  geom_line() +              # Add lines connecting the same genotype across years
  geom_point() +             # Add points for each year/genotype
  labs(x = "expe", y = "Larvae number", title = "confirmation of the genotypes A.fabae") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels for better readability

# add the genetic cluster info ----
###################################

df_long_genetic2 <- read.csv("00.info/02.genetic_structure_jonathan_dijon/faba_bean_ivip_and_genetic_cluster.csv", 
                            sep =",", 
                            stringsAsFactors = F) %>% 
  dplyr::select(ivip, cluster_name) %>% 
  distinct() %>%
  inner_join(., df_long2, by = "ivip")

# plot the 8 genotypes selected based on consistency with genetic background

ggplot(df_long_genetic2, aes(x = xp, y = median_larvae_number, group = ivip, color = cluster_name)) +
  geom_line() +              # Add lines connecting the same genotype across years
  geom_point() +             # Add points for each year/genotype
  geom_label(aes(label = ivip), vjust = -0.5, hjust = 0.5, size = 2) + # Add genotype labels
  labs(x = "Experiment", y = "Median larvae number", title = "Genotype consistency (A. fabae fecundity test)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  # Rotate x-axis labels for better readability
  scale_color_manual(values = c("red", "blue","pink","orange", "green", "yellow"))        



# plot consistency of all genotypes from P2XP007 after trimming----
# Plotting with ggplot2 
ggplot(df_long3, aes(x = xp, y = median_larvae_number, group = ivip, color = ivip)) +
  geom_line() +              # Add lines connecting the same genotype across years
  geom_point() +             # Add points for each year/genotype
  labs(x = "expe", y = "Larvae number", title = "confirmation of the genotypes A.fabae") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels for better readability

# add the genetic cluster info ----
###################################

df_long_genetic3 <- read.csv("00.info/02.genetic_structure_jonathan_dijon/faba_bean_ivip_and_genetic_cluster.csv", 
                             sep =",", 
                             stringsAsFactors = F) %>% 
  dplyr::select(ivip, cluster_name) %>% 
  distinct() %>%
  inner_join(., df_long3, by = "ivip")

# plot the 8 genotypes selected based on consistency with genetic background

ggplot(df_long_genetic3, aes(x = xp, y = median_larvae_number, group = ivip, color = cluster_name)) +
  geom_line() +              # Add lines connecting the same genotype across years
  geom_point() +             # Add points for each year/genotype
  geom_label(aes(label = ivip), vjust = -0.5, hjust = 0.5, size = 2) + # Add genotype labels
  labs(x = "Experiment", y = "Median larvae number", title = "Genotype consistency (A. fabae fecundity test)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  # Rotate x-axis labels for better readability
  scale_color_manual(values = c("red", "blue","pink","orange", "green", "yellow"))        








# plot consistency of all genotypes----

# Plotting with ggplot2 
ggplot(df_long, aes(x = xp, y = median_larvae_number, group = ivip, color = ivip)) +
  geom_line() +              # Add lines connecting the same genotype across years
  geom_point() +             # Add points for each year/genotype
  labs(x = "expe", y = "Larvae number", title = "confirmation of the genotypes A.fabae") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels for better readability

###################################
# add the genetic cluster info ----
###################################

df_long_genetic <- read.csv("00.info/02.genetic_structure_jonathan_dijon/faba_bean_ivip_and_genetic_cluster.csv", 
                            sep =",", 
                            stringsAsFactors = F) %>% 
  dplyr::select(ivip, cluster_name) %>% 
  distinct() %>%
  inner_join(., df_long, by = "ivip")

# plot the 8 genotypes selected based on consistency with genetic background

ggplot(df_long_genetic, aes(x = xp, y = median_larvae_number, group = ivip, color = cluster_name)) +
  geom_line() +              # Add lines connecting the same genotype across years
  geom_point() +             # Add points for each year/genotype
  geom_label(aes(label = ivip), vjust = -0.5, hjust = 0.5, size = 2) + # Add genotype labels
  labs(x = "Experiment", y = "Median larvae number", title = "Genotype consistency (A. fabae fecundity test)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  # Rotate x-axis labels for better readability
  scale_color_manual(values = c("red", "blue","pink","orange", "green", "yellow"))        


