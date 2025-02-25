# Packages ----
suppressPackageStartupMessages(library(corrplot))
suppressPackageStartupMessages(library(Hmisc))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(agricolae))

#################################################################################
#################################################################################


# Import all fecundity measurements from 2022-2023 (N=251 genotypes), P2-XP006 (N=16 genotypes), P2-XP008 (N=7 geno.) and P2-XP010 (N=7 geno.)
source("02.analyses/02.Discriminant_Genotypes/2025_consistency_analysis/acyrthosiphon_pisum_P123/prepare_data_from_bigscreen_P2XP006_P2XP008_P2XP010.R")


####### ploting median consistency ##########################################


# plot consistency of all genotypes from P2XP006----
# Plotting with ggplot2 
ggplot(df_long2, aes(x = xp, y = median_larvae_number, group = ivip, color = ivip)) +
  geom_line() +              # Add lines connecting the same genotype across years
  geom_point() +             # Add points for each year/genotype
  labs(x = "expe", y = "Larvae number", title = "confirmation of the genotypes A.pisum") +
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
  labs(x = "Experiment", y = "Median larvae number", title = "Genotype consistency (A. pisum fecundity test)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  # Rotate x-axis labels for better readability
  scale_color_manual(values = c("North" = "red", "South" = "blue", "Paucijuga" = "pink","South_America" = "orange", "Winter" = "green", "East" = "yellow"))              



# plot consistency of all genotypes from P2XP006 after trimming----
# Plotting with ggplot2 
ggplot(df_long3, aes(x = xp, y = median_larvae_number, group = ivip, color = ivip)) +
  geom_line() +              # Add lines connecting the same genotype across years
  geom_point() +             # Add points for each year/genotype
  labs(x = "expe", y = "Larvae number", title = "confirmation of the genotypes A.pisum") +
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
  labs(x = "Experiment", y = "Median larvae number", title = "Genotype consistency (A. pisum fecundity test)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  # Rotate x-axis labels for better readability
  scale_color_manual(values = c("North" = "red", "South" = "blue", "Paucijuga" = "pink","South_America" = "orange", "Winter" = "green", "East" = "yellow"))



# plot consistency of all genotypes and all experiments----

# Plotting with ggplot2 
ggplot(df_long_4_expe, aes(x = xp, y = median_larvae_number, group = ivip, color = ivip)) +
  geom_line() +              # Add lines connecting the same genotype across years
  geom_point() +             # Add points for each year/genotype
  labs(x = "expe", y = "Larvae number", title = "confirmation of the genotypes A.pisum") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels for better readability

###################################
# add the genetic cluster info ----
###################################

df_long_genetic_4_expe <- read.csv("00.info/02.genetic_structure_jonathan_dijon/faba_bean_ivip_and_genetic_cluster.csv", 
                            sep =",", 
                            stringsAsFactors = F) %>% 
  dplyr::select(ivip, cluster_name) %>% 
  distinct() %>%
  inner_join(., df_long_4_expe, by = "ivip")

# plot the 8 genotypes selected based on consistency with genetic background

ggplot(df_long_genetic_4_expe, aes(x = xp, y = median_larvae_number, group = ivip, color = cluster_name)) +
  geom_line() +              # Add lines connecting the same genotype across years
  geom_point() +             # Add points for each year/genotype
  geom_label(aes(label = ivip), vjust = -0.5, hjust = 0.5, size = 2) + # Add genotype labels
  labs(x = "Experiment", y = "Median larvae number", title = "Genotype consistency (A. pisum fecundity test)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  # Rotate x-axis labels for better readability
  scale_color_manual(values = c("North" = "red", "South" = "blue", "Paucijuga" = "pink",                                  "South_America" = "orange", "cluster5" = "green", "cluster6" = "yellow"))              







#Remove some genotypes based on consistency  ----

# Remove rows 
data_median_expe20 <- data_median_XP006 %>% filter( ivip != 106)
data_median_expe20 <- data_median_expe20 %>% filter( ivip != 155)

#extract the 4 highest and lowest
top_4_lowest <- data_median_expe20 %>%
  arrange(.[[2]]) %>%  # Arrange by the second column (ascending order)
  slice_head(n = 4)  # Select the top 4 rows

top_4_highest <- data_median_expe20 %>%
  arrange(desc(.[[2]])) %>%  # Arrange by the second column (descending order)
  slice_head(n = 4)  # Select the top 4 rows

# Combine the two datasets
median_8_20 <- bind_rows(
  top_4_lowest %>% mutate(group = "lowest"),  # Add a group identifier for clarity
  top_4_highest %>% mutate(group = "highest") # Add a group identifier for clarity
)




median_8_250 <- data_median_expe250 %>%
  filter(ivip %in% median_8_20$ivip)

combined_dataxp006 <- merge(median_8_20, median_8_250, by = "ivip", all = TRUE)

df_longxp006 <- combined_dataxp006 %>%
  pivot_longer(cols = starts_with("expe"), 
               names_to = "xp", 
               names_prefix = "expe",
               values_to = "median_larvae_number")



df_long_geneticxp006 <- read.csv("00.info/02.genetic_structure_jonathan_dijon/faba_bean_ivip_and_genetic_cluster.csv", 
                             sep =",", 
                             stringsAsFactors = F) %>% 
  dplyr::select(ivip, cluster_name) %>% 
  distinct() %>%
  inner_join(., df_longxp006, by = "ivip")

ggplot(df_long_geneticxp006, aes(x = xp, y = median_larvae_number, group = ivip, color = cluster_name)) +
  geom_line() +              # Add lines connecting the same genotype across years
  geom_point() +             # Add points for each year/genotype
  geom_label(aes(label = ivip), vjust = -0.5, hjust = 0.5, size = 2) + # Add genotype labels
  labs(x = "Experiment", y = "Median larvae number", title = "Genotype consistency (A. pisum fecundity test)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  # Rotate x-axis labels for better readability
  scale_color_manual(values = c("North" = "red", "South" = "blue", "Paucijuga" = "pink","South_America" = "orange", "Winter" = "green", "East" = "yellow"))       







df_long_geneticxp008 <- read.csv("00.info/02.genetic_structure_jonathan_dijon/faba_bean_ivip_and_genetic_cluster.csv", 
                                 sep =",", 
                                 stringsAsFactors = F) %>% 
  dplyr::select(ivip, cluster_name) %>% 
  distinct() %>%
  inner_join(., df_long_250_xp8, by = "ivip")



ggplot(df_long_geneticxp008, aes(x = xp, y = median_larvae_number, group = ivip, color = cluster_name)) +
  geom_line() +              # Add lines connecting the same genotype across years
  geom_point() +             # Add points for each year/genotype
  geom_label(aes(label = ivip), vjust = -0.5, hjust = 0.5, size = 2) + # Add genotype labels
  labs(x = "Experiment", y = "Median larvae number", title = "Genotype consistency (A. pisum fecundity test)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  # Rotate x-axis labels for better readability
  scale_color_manual(values = c("North" = "red", "South" = "blue", "Paucijuga" = "pink","South_America" = "orange", "Winter" = "green", "East" = "yellow"))       





