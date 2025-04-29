
#data prep screening of the 250 genotypes fabae ----
data250_fabae <- read.csv(
  file = "01.data/2023_screen_specifics_251genotypes/02.aphis_fabae_on_faba_collection/aphis_fabae_fecundity_on_faba_collection_2024-04-02.csv", 
  sep =",", 
  stringsAsFactors = F)

# Remove rows with missing values in 'column_name'
data250_fabae <- data250_fabae %>% filter(!is.na(larvae_total_number))


data250_fabae$ivip<-as.factor(data250_fabae$ivip)
data_filtered_250_fabae <- data250_fabae %>% filter(adult_alive != "0")
summary(data_filtered_250_fabae)

data_filtered_250_inf6_fabae <- data_filtered_250_fabae %>%
  group_by(ivip) %>%           # Group by genotype
  filter(n() >= 6) %>%             # Keep only groups with 6 or more repetitions
  ungroup()                        # Ungroup to return to original data structure

summary(data_filtered_250_inf6_fabae)
data_filtered_250_inf6_fabae$session<-as.factor(data_filtered_250_inf6_fabae$session)

#data prep screening of the 250 genotypes pisum ----
data250_pisum <- read.csv(
  file = "01.data/2023_screen_specifics_251genotypes/01.acyrthosiphon_pisum_P123amp_on_faba_collection/acyrthospipon_pisum_fecundity_on_faba_collection.csv", 
  sep =",", 
  stringsAsFactors = F)

# Remove rows with missing values in 'column_name'
data250_pisum <- data250_pisum %>% filter(!is.na(larvae_total_number))


data250_pisum$ivip<-as.factor(data250_pisum$ivip)
data_filtered_250_pisum <- data250_pisum %>% filter(adult_alive != "0")
summary(data_filtered_250_pisum)

data_filtered_250_inf6_pisum <- data_filtered_250_pisum %>%
  group_by(ivip) %>%           # Group by genotype
  filter(n() >= 6) %>%             # Keep only groups with 6 or more repetitions
  ungroup()                        # Ungroup to return to original data structure

summary(data_filtered_250_inf6_pisum)
data_filtered_250_inf6_pisum$session<-as.factor(data_filtered_250_inf6_pisum$session)



#median extract ----
data_median_expe250_fabae <- data_filtered_250_inf6_fabae %>%
  group_by(ivip) %>%
  summarise(Median_fecondity = median(larvae_total_number))
colnames(data_median_expe250_fabae)[colnames(data_median_expe250_fabae) == "Median_fecondity"] <- "fabae"


# df_long_genetic3 <- read.csv("00.info/02.genetic_structure_jonathan_dijon/faba_bean_ivip_and_genetic_cluster.csv", 
#                              sep =",", 
#                              stringsAsFactors = F) %>% 
#   dplyr::select(ivip, cluster_name) %>% 
#   distinct() %>%
#   inner_join(., df_long3, by = "ivip")

#median extract ----
data_median_expe250_pisum <- data_filtered_250_inf6_pisum %>%
  group_by(ivip) %>%
  summarise(Median_fecondity = median(larvae_total_number))
colnames(data_median_expe250_pisum)[colnames(data_median_expe250_pisum) == "Median_fecondity"] <- "pisum"








#combine data ----

combined_data <- merge(data_median_expe250_pisum, data_median_expe250_fabae, by = "ivip", all = FALSE)


combined_data_genetcluster<- read.csv("00.info/02.genetic_structure_jonathan_dijon/faba_bean_ivip_and_genetic_cluster.csv", 
                             sep =",", 
                             stringsAsFactors = F) %>% 
  dplyr::select(ivip, cluster_name) %>% 
  distinct() %>%
  inner_join(., combined_data, by = "ivip")




# Filter the data for the "North" cluster
north_data <- subset(combined_data_genetcluster, cluster_name == "North")

# Calculate the correlation coefficient
correlation_coefficient <- cor(north_data$pisum, north_data$fabae)
print(paste("Correlation coefficient: ", correlation_coefficient))



ggplot(north_data, aes(x = pisum, y = fabae)) +
  geom_point() +
  geom_smooth(method = "lm", col = "blue", se = FALSE) + # Add regression line
  labs(title = paste("Scatterplot of Pisum vs Fabae for North Cluster\nCorrelation: ", round(correlation_coefficient, 2)),
       x = "Pisum",
       y = "Fabae") +
  theme_minimal()

# Filter the data for the "South" cluster
south_data <- subset(combined_data_genetcluster, cluster_name == "South")

# Calculate the correlation coefficient
correlation_coefficient <- cor(south_data$pisum, south_data$fabae)
print(paste("Correlation coefficient: ", correlation_coefficient))



ggplot(north_data, aes(x = pisum, y = fabae)) +
  geom_point() +
  geom_smooth(method = "lm", col = "blue", se = FALSE) + # Add regression line
  labs(title = paste("Scatterplot of Pisum vs Fabae for North Cluster\nCorrelation: ", round(correlation_coefficient, 2)),
       x = "Pisum",
       y = "Fabae") +
  theme_minimal()



