#plot des 4 test meme si les genotypes ne sont pas le meme

# etape 1: combiner les genotypes 

#attention au genotypes non commun au verifications

# etape 2 plot avec et sans fond genetique

#etape 3 combiner les dataset effet BBCH et bloc

#etape 4 difference significatives



# Packages ----
suppressPackageStartupMessages(library(corrplot))
suppressPackageStartupMessages(library(Hmisc))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))


#data prep screening of the 250 genotypes ----
data250<-read.csv("../01.data/2023_screen_specifics_251genotypes/01.acyrthosiphon_pisum_P123amp_on_faba_collection/acyrthospipon_pisum_fecundity_on_faba_collection.csv", sep =",")


# Remove rows with missing values in 'column_name'
data250 <- data250 %>% filter(!is.na(larvae_total_number))


data250$ivip<-as.factor(data250$ivip)
data_filtered_250 <- data250 %>% filter(adult_alive != "0")
summary(data_filtered_250)

data_filtered_250_inf6 <- data_filtered_250 %>%
  group_by(ivip) %>%           # Group by genotype
  filter(n() >= 6) %>%             # Keep only groups with 6 or more repetitions
  ungroup()                        # Ungroup to return to original data structure

summary(data_filtered_250_inf6)
data_filtered_250_inf6$session<-as.factor(data_filtered_250_inf6$session)

# plot the median of the 250 genotypes ----

median_hist_data <- data_filtered_250_inf6 %>%
  group_by(ivip) %>%
  summarise(Median_fecondity = median(larvae_total_number))
hist(median_hist_data$Median_fecondity)


median_hist_data$ivip <- reorder(median_hist_data$ivip, median_hist_data$Median_fecondity)


ggplot(median_hist_data, aes(x = ivip, y = Median_fecondity)) +
  geom_bar(stat = "identity", fill = "Skyblue") +  # Use fill color for the bars
  labs(title = "ajusted mean of number of larvae per genotype",
       x = "Genotype",
       y = "Number of larvae") +
  theme_minimal() +
  theme(axis.text.x = element_blank(),  # Hide x-axis labels for readability
        axis.ticks.x = element_blank()) +
  scale_y_continuous(limits = c(0, 100), expand = c(0, 0))  # Set y-axis limits as needed







#data prep P2XP006 ----

data_xp006<-read.csv("../01.data/P2-XP006_faba_vs_pisum_P123_fecundity-test/P2-XP006_larvae_fecundity_results.csv", sep =",")
data_xp006 <- data_xp006 %>% filter(!is.na(larvae_total_number))
data_filtered_xp006 <- data_xp006 %>% filter(adult_alive != "0")

data_filtered_xp006$ivip<-as.factor(data_filtered_xp006$ivip)
summary(data_filtered_xp006)
data_filtered_xp006_inf6 <- data_filtered_xp006 %>%
  group_by(ivip) %>%           # Group by genotype
  filter(n() >= 6) %>%             # Keep only groups with 6 or more repetitions
  ungroup()                        # Ungroup to return to original data structure
summary(data_filtered_xp006_inf6)


#trimming data 250 for comparison with P2XP006----

# Using filter with %in% to keep rows in df1 where column_a matches any value in column_b of df2
data250_trimmed1 <- data_filtered_250_inf6 %>%
  filter(ivip %in% data_filtered_xp006$ivip)

# median extract 250 for comparison with P2XP006 ----


data_median_expe250 <- data250_trimmed1 %>%
  group_by(ivip) %>%
  summarise(Median_fecondity = median(larvae_total_number))
colnames(data_median_expe250)[colnames(data_median_expe250) == "Median_fecondity"] <- "expe1"
#median extract P2X006 ----

data_median_expe20 <- data_filtered_xp006_inf6 %>%
  group_by(ivip) %>%
  summarise(Median_fecondity = median(larvae_total_number))
colnames(data_median_expe20)[colnames(data_median_expe20) == "Median_fecondity"] <- "expe2"

# merging of the newly created file for the screening and P2XP006 ----
# Merge the two data frames by Genotype, ensuring all genotypes are kept
combined_data <- merge(data_median_expe250, data_median_expe20, by = "ivip", all = TRUE)

colnames(combined_data) <- trimws(tolower(colnames(combined_data)))  # Trim spaces and make lowercase


# Reshape the data into long format
df_long <- combined_data %>%
  pivot_longer(cols = starts_with("expe"), 
               names_to = "test", 
               names_prefix = "expe",
               values_to = "expe")

# View the reshaped data
print(df_long, n=250)

# plot consistency of all genotypes from P2XP006----


# Plotting with ggplot2 
ggplot(df_long, aes(x = test, y = expe, group = ivip, color = ivip)) +
  geom_line() +              # Add lines connecting the same genotype across years
  geom_point() +             # Add points for each year/genotype
  labs(x = "expe", y = "Larvae number", title = "confirmation of the genotypes A.pisum") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels for better readability


#Remove some genotypes based on consistency  ----

# Remove rows 
data_median_expe20 <- data_median_expe20 %>% filter( ivip != 106)
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

print(combined_data)


median_8_250 <- data_median_expe250 %>%
  filter(ivip %in% median_8_20$ivip)

combined_data2 <- merge(median_8_20, median_8_250, by = "ivip", all = TRUE)

df_long2 <- combined_data2 %>%
  pivot_longer(cols = starts_with("expe"), 
               names_to = "test", 
               names_prefix = "expe",
               values_to = "expe")

# View the reshaped data
print(df_long, n=250)

#add the genetic cluster info ----

datagen<-read.csv("../00.info/02.genetic_structure_jonathan_dijon/faba_bean_ivip_and_genetic_cluster.csv", sep =",")

# Select only the specific column you want to merge from data2
data2 <- datagen[, c(1,8)]
# Remove duplicates from data2
data2_unique <- data2[!duplicated(data2$ivip), ] 

# Merge the two datasets
df_long3 <- merge(df_long2,data2_unique, by = "ivip", all.x = TRUE)

#plot the 8 genotypes selected based on consistency ----
#plot the 8 genotypes selected based on consistency with genetic background

ggplot(df_long3, aes(x = test, y = expe, group = ivip, color = cluster_name)) +
  geom_line() +              # Add lines connecting the same genotype across years
  geom_point() +             # Add points for each year/genotype
  geom_label(aes(label = ivip), vjust = -0.5, hjust = 0.5, size = 2) + # Add genotype labels
  labs(x = "expe", y = "Larvae number", title = "Confirmation of the genotypes, A.pisum") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  # Rotate x-axis labels for better readability
  scale_color_manual(values = c("red", "blue","pink","orange", "green", "yellow"))        


#plot the 8 genotypes selected based on consistency without genetic background
ggplot(df_long2, aes(x = test, y = expe, group = ivip, color = ivip)) +
  geom_line() +              # Add lines connecting the same genotype across years
  geom_point() +             # Add points for each year/genotype
  labs(x = "expe", y = "Larvae number", title = "confirmation of the genotypes, A.pisum") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels for better readability



#data prep P2XP008 ----
  
  data_xp008<-read.csv("../01.data/P2-XP008_faba_apisum_confirmation_2/P2-XP008_fecundity_data_apisum.csv", sep =",")
data_xp008 <- data_xp008 %>% filter(!is.na(larvae_total_number))
data_xp008$ivip <-as.factor(data_xp008$ivip)
data_xp008$ivip <-as.factor(data_xp008$ivip)
data_xp008$bbch_25days <-as.integer(data_xp008$bbch_25days)

data_xp008$genotype <-as.factor(data_xp008$genotype)
data_filtered_xp008 <- data_xp008 %>% filter(adult_alive != "0")


data_filtered_xp008$ivip<-as.factor(data_filtered_xp008$ivip)
summary(data_filtered_xp008)
data_filtered_xp008_inf6 <- data_filtered_xp008 %>%
  group_by(ivip) %>%           # Group by genotype
  filter(n() >= 6) %>%             # Keep only groups with 6 or more repetitions
  ungroup()                        # Ungroup to return to original data structure
summary(data_filtered_xp008_inf6)

#annalyse bloc effect ----
# larvae number explained by other factor --> effect of the bbch when counting and growth bloc
mod <-lm(larvae_total_number~bbch_10days, data = data_filtered_xp008_inf6)
anova(mod)

mod <-lm(larvae_total_number~bbch_25days, data = data_filtered_xp008_inf6)
anova(mod)

mod <-lm(larvae_total_number~genotype, data = data_filtered_xp008_inf6)
anova(mod)

mod <-lm(larvae_total_number~bbch_10days, data = data_filtered_xp008_inf6)
anova(mod)

mod <-lm(larvae_total_number~block, data = data_filtered_xp008_inf6)
anova(mod)

mod <-lm(larvae_total_number~installation_block, data = data_filtered_xp008_inf6)
anova(mod)

# BBCH
mod <-lm(bbch_25days~block, data = data_filtered_xp008_inf6)
anova(mod)

mod <-lm(bbch_25days~genotype, data = data_filtered_xp008_inf6)
anova(mod)

#regression
step_model <- step(lm(larvae_total_number~installation_block+block+bbch_25days, data = data_filtered_xp008_inf6), direction = "both")
summary(step_model)


mod <- lm(larvae_total_number~genotype+block+bbch_25days, data = data_filtered_xp008_inf6)
anova(mod)

mod_anova <- aov(larvae_total_number ~ ivip + block + bbch_25days, 
                 data = data_filtered_xp008_inf6)

tukey_groups <- HSD.test(mod_anova, "ivip", console = TRUE)


#median extract
data_median_xp08 <- data_filtered_xp008_inf6 %>%
  group_by(ivip) %>%
  summarise(Median_fecondity = median(larvae_total_number))
colnames(data_median_expe20)[colnames(data_median_expe20) == "Median_fecondity"] <- "expe2"




#Combining P2XP006 and P2XP008 ----
data_combined_xp6_8<- merge(data_xp008, data_xp006)

#data prep P2XP010 ----

data_xp010<-read.csv("../01.data/P2-XP010_fabae_pisum_third_confirmation/P2XP010.csv", sep =",")
data_xp007 <- data_xp010 %>% filter(!is.na(larvae_total_number))
data_filtered_xp008 <- data_xp008 %>% filter(adult_alive != "0")

data_filtered_xp008$ivip<-as.factor(data_filtered_xp008$ivip)
summary(data_filtered_xp008)
data_filtered_xp008_inf6 <- data_filtered_xp008 %>%
  group_by(ivip) %>%           # Group by genotype
  filter(n() >= 6) %>%             # Keep only groups with 6 or more repetitions
  ungroup()                        # Ungroup to return to original data structure
summary(data_filtered_xp008_inf6)