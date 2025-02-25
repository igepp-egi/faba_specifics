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
  scale_color_manual(values = c("North" = "red", "South" = "blue", "Paucijuga" = "pink","South_America" = "orange", "Winter" = "green", "East" = "yellow"))        



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
  scale_color_manual(values = c("North" = "red", "South" = "blue", "Paucijuga" = "pink","South_America" = "orange", "Winter" = "green", "East" = "yellow"))        








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
  scale_color_manual(values = c("North" = "red", "South" = "blue", "Paucijuga" = "pink","South_America" = "orange", "Winter" = "green", "East" = "yellow"))        



#Remove some genotypes based on consistency  ----

# Remove rows 
data_median_expe20 <- data_median_xp007 %>% filter( ivip != "IVIP101738")
data_median_expe20 <- data_median_expe20 %>% filter( ivip != "IVIP101240")
data_median_expe20 <- data_median_expe20 %>% filter( ivip != "IVIP101482")



#extract the 4 highest and lowest
top_4_lowest <- data_median_expe20 %>%
  arrange(.[[2]]) %>%  # Arrange by the second column (ascending order)
  slice_head(n = 5)  # Select the top 4 rows

top_4_highest <- data_median_expe20 %>%
  arrange(desc(.[[2]])) %>%  # Arrange by the second column (descending order)
  slice_head(n = 5)  # Select the top 4 rows

# Combine the two datasets
median_8_20 <- bind_rows(
  top_4_lowest %>% mutate(group = "lowest"),  # Add a group identifier for clarity
  top_4_highest %>% mutate(group = "highest") # Add a group identifier for clarity
)




median_8_250 <- data_median_expe250 %>%
  filter(ivip %in% median_8_20$ivip)

combined_dataxp007 <- merge(median_8_20, median_8_250, by = "ivip", all = TRUE)

df_longxp007 <- combined_dataxp007 %>%
  pivot_longer(cols = starts_with("expe"), 
               names_to = "xp", 
               names_prefix = "expe",
               values_to = "median_larvae_number")

# View the reshaped data
print(df_long, n=250)

df_long_geneticxp007 <- read.csv("00.info/02.genetic_structure_jonathan_dijon/faba_bean_ivip_and_genetic_cluster.csv", 
                                 sep =",", 
                                 stringsAsFactors = F) %>% 
  dplyr::select(ivip, cluster_name) %>% 
  distinct() %>%
  inner_join(., df_longxp007, by = "ivip")

ggplot(df_long_geneticxp007, aes(x = xp, y = median_larvae_number, group = ivip, color = cluster_name)) +
  geom_line() +              # Add lines connecting the same genotype across years
  geom_point() +             # Add points for each year/genotype
  geom_label(aes(label = ivip), vjust = -0.5, hjust = 0.5, size = 2) + # Add genotype labels
  labs(x = "Experiment", y = "Median larvae number", title = "Genotype consistency (A. fabae fecundity test)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  # Rotate x-axis labels for better readability
  scale_color_manual(values = c("North" = "red", "South" = "blue", "Paucijuga" = "pink","South_America" = "orange", "Winter" = "green", "East" = "yellow"))       




