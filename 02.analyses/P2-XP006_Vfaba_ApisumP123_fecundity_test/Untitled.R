
suppressPackageStartupMessages(library(corrplot))
suppressPackageStartupMessages(library(Hmisc))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))


#data prep 250 ----
data250<-read.csv("../faba_specifics/01.data/2023_screen_specifics_251genotypes/01.acyrthosiphon_pisum_P123amp_on_faba_collection/acyrthospipon_pisum_fecundity_on_faba_collection.csv", sep =",")


# Remove rows with missing values in 'column_name'
data250 <- data250 %>% filter(!is.na(larvae_total_number))


data250$short_code<-as.factor(data250$short_code)
data_filtered_250 <- data250 %>% filter(adult_alive != "0")
summary(data_filtered_250)

data_filtered_250_inf6 <- data_filtered_250 %>%
  group_by(short_code) %>%           # Group by genotype
  filter(n() >= 6) %>%             # Keep only groups with 6 or more repetitions
  ungroup()                        # Ungroup to return to original data structure

summary(data_filtered_250_inf6)
data_filtered_250_inf6$session<-as.factor(data_filtered_250_inf6$session)

# plot 250 median ----

median_hist_data <- data_filtered_250_inf6 %>%
  group_by(short_code) %>%
  summarise(Median_fecondity = median(larvae_total_number))
hist(median_hist_data$Median_fecondity)




median_hist_data$short_code <- reorder(median_hist_data$short_code, median_hist_data$Median_fecondity)


ggplot(median_hist_data, aes(x = short_code, y = Median_fecondity)) +
  geom_bar(stat = "identity", fill = "Skyblue") +  # Use fill color for the bars
  labs(title = "ajusted mean of number of larvae per genotype",
       x = "Genotype",
       y = "Number of larvae") +
  theme_minimal() +
  theme(axis.text.x = element_blank(),  # Hide x-axis labels for readability
        axis.ticks.x = element_blank()) +
  scale_y_continuous(limits = c(0, 100), expand = c(0, 0))  # Set y-axis limits as needed







#### data prep confirmation set ----

data_20<-read.csv("01.data/P2-XP006_faba_vs_pisum_P123_fecundity-test/P2-XP006_larvae_fecundity_results.csv", sep =";")
data_20 <- data_20 %>% filter(!is.na(larvae_total_number))
data_20$aphid_species<-as.factor(data_20$aphid_species)
data_filtered_20 <- data_20 %>% filter(adult_alive != "0")

data_filtered_20$short_code<-as.factor(data_filtered_20$short_code)
summary(data_filtered_20)
data_filtered_20_inf6 <- data_filtered_20 %>%
  group_by(short_code) %>%           # Group by genotype
  filter(n() >= 6) %>%             # Keep only groups with 6 or more repetitions
  ungroup()                        # Ungroup to return to original data structure
summary(data_filtered_20_inf6)


#----

# Using filter with %in% to keep rows in df1 where column_a matches any value in column_b of df2
result <- data_filtered_250_inf6 %>%
  filter(short_code %in% data_filtered_20$short_code)

# median extract 250 ----


data_median_expe250 <- result %>%
  group_by(short_code) %>%
  summarise(Median_fecondity = median(larvae_total_number))
colnames(data_median_expe250)[colnames(data_median_expe250) == "Median_fecondity"] <- "expe1"
#median extract confirmation set ----

data_median_expe20 <- data_filtered_20_inf6 %>%
  group_by(short_code) %>%
  summarise(Median_fecondity = median(larvae_total_number))
colnames(data_median_expe20)[colnames(data_median_expe20) == "Median_fecondity"] <- "expe2"


# Merge the two data frames by Genotype, ensuring all genotypes are kept
combined_data <- merge(data_median_expe250, data_median_expe20, by = "short_code", all = TRUE)

colnames(combined_data) <- trimws(tolower(colnames(combined_data)))  # Trim spaces and make lowercase


# Reshape the data into long format
df_long <- combined_data %>%
  pivot_longer(cols = starts_with("expe"), 
               names_to = "test", 
               names_prefix = "expe",
               values_to = "expe")

# View the reshaped data
print(df_long, n=250)



# Plotting with ggplot2
ggplot(df_long, aes(x = test, y = expe, group = short_code, color = short_code)) +
  geom_line() +              # Add lines connecting the same genotype across years
  geom_point() +             # Add points for each year/genotype
  labs(x = "expe", y = "Larvae number", title = "confirmation of the genotypes A.pissum") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels for better readability



#plot the 8 selected ----
genotype_medians_test <- data_filtered_250_inf6 %>%
  group_by(short_code) %>%
  summarise(Median_fecondity = median(larvae_total_number))
# Remove rows where 'value' is equal to 5 or 'id' is equal to 3
data_median_expe20 <- data_median_expe20 %>% filter( short_code != 106)
data_median_expe20 <- data_median_expe20 %>% filter( short_code != 155)

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

print(combined_data)


median_8_250 <- data_median_expe250 %>%
  filter(short_code %in% median_8_20$short_code)

combined_data2 <- merge(median_8_20, median_8_250, by = "short_code", all = TRUE)

df_long2 <- combined_data2 %>%
  pivot_longer(cols = starts_with("expe"), 
               names_to = "test", 
               names_prefix = "expe",
               values_to = "expe")

# View the reshaped data
print(df_long, n=250)

# Select only the specific column you want to merge from data2
data2 <- data_filtered_20_inf6[, c(1,13)]
# Remove duplicates from data2
data2_unique <- data2[!duplicated(data2$short_code), ] 

# Merge the two datasets
df_long3 <- merge(data2_unique, df_long2, by = "short_code", all.x = TRUE)


ggplot(df_long3, aes(x = test, y = expe, group = short_code, color = Type)) +
  geom_line() +              # Add lines connecting the same genotype across years
  geom_point() +             # Add points for each year/genotype
  geom_label(aes(label = short_code), vjust = -0.5, hjust = 0.5, size = 2) + # Add genotype labels
  labs(x = "expe", y = "Larvae number", title = "Confirmation of the genotypes, A.pissum") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  # Rotate x-axis labels for better readability
  scale_color_manual(values = c("red", "blue","pink","orange", "green", "yellow"))        








# Plotting with ggplot2
ggplot(df_long2, aes(x = test, y = expe, group = short_code, color = short_code)) +
  geom_line() +              # Add lines connecting the same genotype across years
  geom_point() +             # Add points for each year/genotype
  labs(x = "expe", y = "Larvae number", title = "confirmation of the genotypes, A.pissum") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels for better readability


