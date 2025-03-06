suppressPackageStartupMessages(library("tidyverse"))
suppressPackageStartupMessages(library("lme4"))
suppressPackageStartupMessages(library("Hmisc"))
suppressPackageStartupMessages(library("corrplot"))
suppressPackageStartupMessages(library("emmeans"))
suppressPackageStartupMessages(library("performance"))
suppressPackageStartupMessages(library("MASS"))

#data prep screening of the 250 genotypes ----
data250 <- read.csv(
  file = "01.data/2023_screen_specifics_251genotypes/02.aphis_fabae_on_faba_collection/aphis_fabae_fecundity_on_faba_collection_2024-04-02.csv", 
  sep =",", 
  stringsAsFactors = F)

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

#data prep P2xp007 ----

data_xp007<-read.csv("01.data/P2-XP007_fabae_AF001_confirmation/P2-XP007_confirmation_Afabae_10vs10genotypes.csv", stringsAsFactors = F)

data_xp007 <- data_xp007 %>% filter(!is.na(larvae_total_number))
data_filtered_xp007 <- data_xp007 %>% filter(adult_alive != "0")

data_filtered_xp007$ivip<-as.factor(data_filtered_xp007$ivip)
summary(data_filtered_xp007)
data_filtered_xp007_inf6 <- data_filtered_xp007 %>%
  group_by(ivip) %>%           # Group by genotype
  filter(n() >= 6) %>%             # Keep only groups with 6 or more repetitions
  ungroup()                        # Ungroup to return to original data structure
summary(data_filtered_xp007_inf6)


#data prep P2xp009 ----

data_xp009<-read.csv("01.data/P2-XP009_afabae_third_confirmation/P2-XP009_fecundity_data_afabae.csv", stringsAsFactors = F)

data_xp009 <- data_xp009 %>% filter(!is.na(larvae_total_number))
# data_xp009$ivip <-as.factor(data_xp009$ivip)
# data_xp009$ivip <-as.factor(data_xp009$ivip)
data_xp009$bbch_25days <-as.integer(data_xp009$bbch_25days)

# data_xp009$genotype <-as.factor(data_xp009$genotype)
data_filtered_xp009 <- data_xp009 %>% filter(adult_alive != "0")


data_filtered_xp009$ivip<-as.factor(data_filtered_xp009$ivip)
summary(data_filtered_xp009)
data_filtered_xp009_inf6 <- data_filtered_xp009 %>%
  group_by(ivip) %>%           # Group by genotype
  filter(n() >= 6) %>%             # Keep only groups with 6 or more repetitions
  ungroup()                        # Ungroup to return to original data structure
summary(data_filtered_xp009_inf6)

#as.integer
data_filtered_250_inf6$adult_alive<-as.integer(data_filtered_250_inf6$adult_alive)
data_filtered_xp007_inf6$adult_alive<-as.integer(data_filtered_xp007_inf6$adult_alive)
data_filtered_xp009_inf6$adult_alive<-as.integer(data_filtered_xp009_inf6$adult_alive)

data_filtered_250_inf6$adult_winged<-as.integer(data_filtered_250_inf6$adult_winged)
data_filtered_xp007_inf6$adult_winged<-as.integer(data_filtered_xp007_inf6$adult_winged)
data_filtered_xp009_inf6$adult_winged<-as.integer(data_filtered_xp009_inf6$adult_winged)

data_filtered_250_inf6$block<-as.character(data_filtered_250_inf6$block)
data_filtered_xp009_inf6$block<-as.character(data_filtered_xp009_inf6$block)


# median extract 250 for comparison with P2xp007 ----


data_median_expe250 <- data_filtered_250_inf6 %>%
  group_by(ivip) %>%
  summarise(Median_fecondity = median(larvae_total_number))
colnames(data_median_expe250)[colnames(data_median_expe250) == "Median_fecondity"] <- "expe1"
#median extract P2X006 ----

data_median_xp007 <- data_filtered_xp007_inf6 %>%
  group_by(ivip) %>%
  summarise(Median_fecondity = median(larvae_total_number))
colnames(data_median_xp007)[colnames(data_median_xp007) == "Median_fecondity"] <- "expe2"


#median extract P2X008 ----

data_median_xp009 <- data_filtered_xp009_inf6 %>%
  group_by(ivip) %>%
  summarise(Median_fecondity = median(larvae_total_number))
colnames(data_median_xp009)[colnames(data_median_xp009) == "Median_fecondity"] <- "expe3"


# merging the newly created file: 250 genotypes screening, P2xp007, P2xp009 ----
# Merge the data frames by Genotype, ensuring all genotypes are kept

combined_data <- merge(data_median_expe250, data_median_xp007, by = "ivip", all = TRUE)
combined_data <- merge(combined_data, data_median_xp009, by = "ivip", all = TRUE)



combined_data <- combined_data %>%
  rowwise() %>%
  filter(sum(!is.na(c_across(starts_with("expe")))) > 2) %>%
  ungroup()


colnames(combined_data) <- trimws(tolower(colnames(combined_data)))  # Trim spaces and make lowercase

# Reshape the data into long format
df_long <- combined_data %>%
  pivot_longer(cols = starts_with("expe"), 
               names_to = "xp", 
               names_prefix = "expe",
               values_to = "median_larvae_number")

# View the reshaped data
print(df_long, n=250)



#merging 250screening and P2XP007 ----

combined_data2 <- merge(data_median_expe250, data_median_xp007, by = "ivip", all = TRUE)

combined_data2 <- combined_data2 %>%
  rowwise() %>%
  filter(sum(!is.na(c_across(starts_with("expe")))) > 1) %>%
  ungroup()

colnames(combined_data2) <- trimws(tolower(colnames(combined_data2)))  # Trim spaces and make lowercase
# Reshape the data into long format
df_long2 <- combined_data2 %>%
  pivot_longer(cols = starts_with("expe"), 
               names_to = "xp", 
               names_prefix = "expe",
               values_to = "median_larvae_number")

# View the reshaped data
print(df_long, n=250)

#-----
combined_data3<- combined_data[,c(1,2,3)]


colnames(combined_data3) <- trimws(tolower(colnames(combined_data3)))  # Trim spaces and make lowercase

# Reshape the data into long format
df_long3 <- combined_data3 %>%
  pivot_longer(cols = starts_with("expe"), 
               names_to = "xp", 
               names_prefix = "expe",
               values_to = "median_larvae_number")

# View the reshaped data
print(df_long, n=250)

#import cluster ----
infos_cluster<- read.csv("00.info/02.genetic_structure_jonathan_dijon/faba_bean_ivip_and_genetic_cluster.csv", 
                         sep =",", 
                         stringsAsFactors = F) %>% 
  dplyr::select(ivip, cluster_name) %>% 
  distinct()
