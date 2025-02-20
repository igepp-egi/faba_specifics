#packages ----
suppressPackageStartupMessages(library(corrplot))
suppressPackageStartupMessages(library(Hmisc))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(agricolae))


#data prep 250 ----
data250<-read.csv("../../../../faba_specifics/01.data/2023_screen_specifics_251genotypes/02.aphis_fabae_on_faba_collection/aphis_fabae_fecundity_on_faba_collection_2024-04-02.csv", sep =",")


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

#### data prep P2XP007 ----

data_xp007<-read.csv("../../../01.data/P2-XP007_fabae_AF001_confirmation/P2-XP007_confirmation_Afabae_10vs10genotypes.csv", sep =",")
data_xp007 <- data_xp007 %>% filter(!is.na(larvae_total_number))
data_filtered_xp007 <- data_xp007 %>% filter(adult_alive != "0")

data_filtered_xp007$ivip<-as.factor(data_filtered_xp007$ivip)
summary(data_filtered_xp007)
data_filtered_xp007_inf6 <- data_filtered_xp007 %>%
  group_by(ivip) %>%           # Group by genotype
  filter(n() >= 6) %>%             # Keep only groups with 6 or more repetitions
  ungroup()                        # Ungroup to return to original data structure
summary(data_filtered_xp007_inf6)

data_filtered_xp007_inf6
#### data prep P2XP009 ----

data_XP009<-read.csv("../../../01.data/P2-XP009_afabae_third_confirmation/P2-XP009_fecundity_data_afabae.csv", sep =",")
data_XP009 <- data_XP009 %>% filter(!is.na(larvae_total_number))
data_filtered_XP009 <- data_XP009 %>% filter(adult_alive != "0")

data_filtered_XP009$ivip<-as.factor(data_filtered_XP009$ivip)
summary(data_filtered_XP009)
data_filtered_XP009_inf6 <- data_filtered_XP009 %>%
  group_by(ivip) %>%           # Group by genotype
  filter(n() >= 6) %>%             # Keep only groups with 6 or more repetitions
  ungroup()                        # Ungroup to return to original data structure
summary(data_filtered_XP009_inf6)



# median extract screening 250 genotypes ----


data_median_expe250 <- data_filtered_250_inf6 %>%
  group_by(ivip) %>%
  summarise(Median_fecondity = median(larvae_total_number))
colnames(data_median_expe250)[colnames(data_median_expe250) == "Median_fecondity"] <- "expe1"
#median extract P2XP007 ----

data_median_XP007 <- data_filtered_xp007_inf6 %>%
  group_by(ivip) %>%
  summarise(Median_fecondity = median(larvae_total_number))
colnames(data_median_XP007)[colnames(data_median_XP007) == "Median_fecondity"] <- "expe2"

#median extract P2XP009 ----

data_median_XP009 <- data_filtered_XP009_inf6 %>%
  group_by(ivip) %>%
  summarise(Median_fecondity = median(larvae_total_number))
colnames(data_median_XP009)[colnames(data_median_XP009) == "Median_fecondity"] <- "expe3"


# merging the newly created file: 250 genotypes screening, P2XP007, and P2XP009 ----
# Merge the data frames by Genotype, ensuring all genotypes are kept

combined_data <- merge(data_median_expe250, data_median_XP007, by = "ivip", all = TRUE)
combined_data <- merge(combined_data, data_median_XP009, by = "ivip", all = TRUE)


combined_data <- combined_data %>%
  rowwise() %>%
  filter(sum(!is.na(c_across(starts_with("expe")))) > 2) %>%
  ungroup()




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
ggplot(df_long, aes(x = test, y = expe, group = ivip, color = ivip)) +
  geom_line() +              # Add lines connecting the same genotype across years
  geom_point() +             # Add points for each year/genotype
  labs(x = "expe", y = "Larvae number", title = "confirmation of the genotypes A.pissum") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels for better readability

# adding the genetic background----

datagen<-read.csv("../../../../faba_specifics/00.info/02.genetic_structure_jonathan_dijon/faba_bean_ivip_and_genetic_cluster.csv", sep =",")

# Select only the specific column you want to merge from data2
data2 <- datagen[, c(1,8)]
# Remove duplicates from data2
data2_unique <- data2[!duplicated(data2$ivip), ] 

# Merge the two datasets
df_long2 <- merge(df_long,data2_unique, by = "ivip", all.x = TRUE)


ggplot(df_long2, aes(x = test, y = expe, group = ivip, color = cluster_name)) +
  geom_line() +              # Add lines connecting the same genotype across years
  geom_point() +             # Add points for each year/genotype
  geom_label(aes(label = ivip), vjust = -0.5, hjust = 0.5, size = 2) + # Add genotype labels
  labs(x = "expe", y = "Larvae number", title = "Confirmation of the genotypes, A.fabae") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  # Rotate x-axis labels for better readability
  scale_color_manual(values = c("red", "blue","pink","orange", "green", "yellow"))        













#################################################################################
#################################################################################

#annalyse P2XP009 ----
data_filtered_XP009_inf6$bloc<-as.factor(data_filtered_XP009_inf6$bloc)
data_filtered_XP009_inf6$bbch_10days<-as.factor(data_filtered_XP009_inf6$bbch_10days)
data_filtered_XP009_inf6$bbch_25days<-as.factor(data_filtered_XP009_inf6$bbch_25days)
data_filtered_XP009_inf6$genotype<-as.factor(data_filtered_XP009_inf6$genotype)


# larvae number explained by other factor --> effect of the bbch when counting and growth bloc
mod <-lm(larvae_total_number~bbch_10days, data = data_filtered_XP009_inf6)
anova(mod)

mod <-lm(larvae_total_number~bbch_25days, data = data_filtered_XP009_inf6)
anova(mod)

mod <-lm(larvae_total_number~genotype, data = data_filtered_XP009_inf6)
anova(mod)

mod <-lm(larvae_total_number~bloc, data = data_filtered_XP009_inf6)
anova(mod)

mod <-lm(larvae_total_number~bloc+bbch_25days+genotype, data = data_filtered_XP009_inf6)
anova(mod)


summary(data_filtered_XP009_inf6)

# BBCH
mod <-lm(bbch_25days~bloc, data = data_filtered_XP009_inf6)
anova(mod)

mod <-lm(bbch_10days~genotype, data = data_filtered_XP009_inf6)
anova(mod)

mod <-lm(bbch_25days~genotype, data = data_filtered_XP009_inf6)
anova(mod)

#regression
step_model <- step(lm(larvae_total_number~ivip+bloc+bbch_25days, data = data_filtered_XP009_inf6), direction = "both")
summary(step_model)


mod <- lm(larvae_total_number~ivip+block+bbch_25days, data = data_filtered_XP009_inf6)
anova(mod)

mod_anova <- aov(larvae_total_number ~ ivip + block + bbch_25days, 
                 data = data_filtered_XP009_inf6)

tukey_groups2 <- HSD.test(mod_anova, "ivip", console = TRUE)
tukey_groups2<-as.data.frame(tukey_groups2$groups)


