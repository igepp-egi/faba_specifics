suppressPackageStartupMessages(library("tidyverse"))
suppressPackageStartupMessages(library("lme4"))
suppressPackageStartupMessages(library("Hmisc"))
suppressPackageStartupMessages(library("corrplot"))
suppressPackageStartupMessages(library("emmeans"))
suppressPackageStartupMessages(library("performance"))
suppressPackageStartupMessages(library("MASS"))


#data prep screening of the 250 genotypes ----
data250 <- read.csv(
  file = "01.data/2023_screen_specifics_251genotypes/01.acyrthosiphon_pisum_P123amp_on_faba_collection/acyrthospipon_pisum_fecundity_on_faba_collection.csv", 
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

#data prep P2XP006 ----

data_xp006<-read.csv("01.data/P2-XP006_faba_vs_pisum_P123_fecundity-test/P2-XP006_larvae_fecundity_results.csv", stringsAsFactors = F)

#data_xp006<-read.csv("../../../01.data/P2-XP006_faba_vs_pisum_P123_fecundity-test/P2-XP006_larvae_fecundity_results.csv", stringsAsFactors = F)

data_xp006 <- data_xp006 %>% filter(!is.na(larvae_total_number))
data_filtered_xp006 <- data_xp006 %>% filter(adult_alive != "0")

data_filtered_xp006$ivip<-as.factor(data_filtered_xp006$ivip)
summary(data_filtered_xp006)
data_filtered_xp006_inf6 <- data_filtered_xp006 %>%
  group_by(ivip) %>%           # Group by genotype
  filter(n() >= 6) %>%             # Keep only groups with 6 or more repetitions
  ungroup()                        # Ungroup to return to original data structure
summary(data_filtered_xp006_inf6)


#data prep P2XP008 ----

data_xp008<-read.csv("01.data/P2-XP008_faba_apisum_confirmation_2/P2-XP008_fecundity_data_apisum.csv", stringsAsFactors = F)
#data_xp008<-read.csv("../../../01.data/P2-XP008_faba_apisum_confirmation_2/P2-XP008_fecundity_data_apisum.csv", stringsAsFactors = F)

data_xp008 <- data_xp008 %>% filter(!is.na(larvae_total_number))
# data_xp008$ivip <-as.factor(data_xp008$ivip)
# data_xp008$ivip <-as.factor(data_xp008$ivip)
data_xp008$bbch_25days <-as.integer(data_xp008$bbch_25days)

# data_xp008$genotype <-as.factor(data_xp008$genotype)
data_filtered_xp008 <- data_xp008 %>% filter(adult_alive != "0")


data_filtered_xp008$ivip<-as.factor(data_filtered_xp008$ivip)
summary(data_filtered_xp008)
data_filtered_xp008_inf6 <- data_filtered_xp008 %>%
  group_by(ivip) %>%           # Group by genotype
  filter(n() >= 6) %>%             # Keep only groups with 6 or more repetitions
  ungroup()                        # Ungroup to return to original data structure
summary(data_filtered_xp008_inf6)



#data prep P2XP010 ----

data_xp010<-read.csv("01.data/P2-XP010_fabae_pisum_third_confirmation/P2XP010.csv", stringsAsFactors = F)
#data_xp010<-read.csv("../../../01.data/P2-XP010_fabae_pisum_third_confirmation/P2XP010.csv", stringsAsFactors = F)


data_xp010 <- data_xp010 %>% filter(!is.na(larvae_total_number))
data_filtered_xp010 <- data_xp010 %>% filter(adult_alive != "0")
data_filtered_xp010$ivip<-as.factor(data_filtered_xp010$ivip)
summary(data_filtered_xp010)

data_filtered_xp010_inf6 <- data_filtered_xp010 %>%
  group_by(ivip) %>%           # Group by genotype
  filter(n() >= 6) %>%             # Keep only groups with 6 or more repetitions
  ungroup()                        # Ungroup to return to original data structure
summary(data_filtered_xp010_inf6)



# median extract 250 for comparison with P2XP006 ----


data_median_expe250 <- data_filtered_250_inf6 %>%
  group_by(ivip) %>%
  summarise(Median_fecondity = median(larvae_total_number))
colnames(data_median_expe250)[colnames(data_median_expe250) == "Median_fecondity"] <- "expe1"
#median extract P2X006 ----

data_median_XP006 <- data_filtered_xp006_inf6 %>%
  group_by(ivip) %>%
  summarise(Median_fecondity = median(larvae_total_number))
colnames(data_median_XP006)[colnames(data_median_XP006) == "Median_fecondity"] <- "expe2"


#median extract P2X008 ----

data_median_XP008 <- data_filtered_xp008_inf6 %>%
  group_by(ivip) %>%
  summarise(Median_fecondity = median(larvae_total_number))
colnames(data_median_XP008)[colnames(data_median_XP008) == "Median_fecondity"] <- "expe3"

#median extract P2X010 ----
data_median_XP010 <- data_filtered_xp010_inf6 %>%
  group_by(ivip) %>%
  summarise(Median_fecondity = median(larvae_total_number))
colnames(data_median_XP010)[colnames(data_median_XP010) == "Median_fecondity"] <- "expe4"






# merging the newly created file: 250 genotypes screening, P2XP006, P2XP008 and P2XP10 ----
# Merge the data frames by Genotype, ensuring all genotypes are kept

combined_data_250_xp6 <- merge(data_median_expe250, data_median_XP006, by = "ivip", all = TRUE)
combined_data_250_xp6_xp8 <- merge(combined_data_250_xp6, data_median_XP008, by = "ivip", all = TRUE)
combined_data_250_xp6_xp8_xp10 <- merge(combined_data_250_xp6_xp8, data_median_XP010, by = "ivip", all = TRUE)


median_combined_250_xp6_xp8_xp10_filtered <- combined_data_250_xp6_xp8_xp10 %>%
  rowwise() %>%
  filter(sum(!is.na(c_across(starts_with("expe")))) > 2) %>%
  ungroup()




colnames(median_combined_250_xp6_xp8_xp10_filtered) <- trimws(tolower(colnames(median_combined_250_xp6_xp8_xp10_filtered)))  # Trim spaces and make lowercase


# Reshape the data into long format
df_long_median_combined_250_xp6_xp8_xp10_filtered <- median_combined_250_xp6_xp8_xp10_filtered %>%
  pivot_longer(cols = starts_with("expe"), 
               names_to = "xp", 
               names_prefix = "expe",
               values_to = "median_larvae_number")



#merging 250screening and P2XP006 ----

median_combined_250_xp6_filtered <- combined_data_250_xp6 %>%
  rowwise() %>%
  filter(sum(!is.na(c_across(starts_with("expe")))) > 1) %>%
  ungroup()

colnames(median_combined_250_xp6_filtered) <- trimws(tolower(colnames(median_combined_250_xp6_filtered)))  # Trim spaces and make lowercase
# Reshape the data into long format
df_long_250_xp6 <- combined_data_250_xp6 %>%
  pivot_longer(cols = starts_with("expe"), 
               names_to = "xp", 
               names_prefix = "expe",
               values_to = "median_larvae_number")




#merging 250screening and P2XP008 ----

median_combined_250_xp8 <- merge(data_median_expe250, data_median_XP008, by = "ivip", all = TRUE)

median_combined_250_xp8_filtered <- median_combined_250_xp8 %>%
  rowwise() %>%
  filter(sum(!is.na(c_across(starts_with("expe")))) > 1) %>%
  ungroup()

colnames(median_combined_250_xp8_filtered) <- trimws(tolower(colnames(median_combined_250_xp8_filtered)))  # Trim spaces and make lowercase
# Reshape the data into long format
df_long_250_xp8 <- median_combined_250_xp8_filtered %>%
  pivot_longer(cols = starts_with("expe"), 
               names_to = "xp", 
               names_prefix = "expe",
               values_to = "median_larvae_number")

#import cluster ----
infos_cluster<- read.csv("00.info/02.genetic_structure_jonathan_dijon/faba_bean_ivip_and_genetic_cluster.csv", 
                         sep =",", 
                         stringsAsFactors = F) %>% 
  dplyr::select(ivip, cluster_name) %>% 
  distinct()

