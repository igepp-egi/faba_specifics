#fix me: faire une annalyse bloc et BBCH avec P2XP008 et 10

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
####### ploting median consistency ##########################################

#data prep screening of the 250 genotypes ----
data250<-read.csv("../01.data/2023_screen_specifics_251genotypes/01.acyrthosiphon_pisum_P123amp_on_faba_collection/acyrthospipon_pisum_fecundity_on_faba_collection.csv", sep =",")

#data250<-read.csv("../../../01.data/2023_screen_specifics_251genotypes/01.acyrthosiphon_pisum_P123amp_on_faba_collection/acyrthospipon_pisum_fecundity_on_faba_collection.csv", sep =",")

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

data_xp006<-read.csv("../01.data/P2-XP006_faba_vs_pisum_P123_fecundity-test/P2-XP006_larvae_fecundity_results.csv", sep =",")

#data_xp006<-read.csv("../../../01.data/P2-XP006_faba_vs_pisum_P123_fecundity-test/P2-XP006_larvae_fecundity_results.csv", sep =",")

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

data_xp008<-read.csv("../01.data/P2-XP008_faba_apisum_confirmation_2/P2-XP008_fecundity_data_apisum.csv", sep =",")
#data_xp008<-read.csv("../../../01.data/P2-XP008_faba_apisum_confirmation_2/P2-XP008_fecundity_data_apisum.csv", sep =",")

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



#data prep P2XP010 ----

data_xp010<-read.csv("../01.data/P2-XP010_fabae_pisum_third_confirmation/P2XP010.csv", sep =",")
#data_xp010<-read.csv("../../../01.data/P2-XP010_fabae_pisum_third_confirmation/P2XP010.csv", sep =",")


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

combined_data <- merge(data_median_expe250, data_median_XP006, by = "ivip", all = TRUE)
combined_data <- merge(combined_data, data_median_XP008, by = "ivip", all = TRUE)
combined_data <- merge(combined_data, data_median_XP010, by = "ivip", all = TRUE)


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

# plot consistency of all genotypes----


# Plotting with ggplot2 
ggplot(df_long, aes(x = test, y = expe, group = ivip, color = ivip)) +
  geom_line() +              # Add lines connecting the same genotype across years
  geom_point() +             # Add points for each year/genotype
  labs(x = "expe", y = "Larvae number", title = "confirmation of the genotypes A.pisum") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels for better readability



#add the genetic cluster info ----

datagen<-read.csv("../00.info/02.genetic_structure_jonathan_dijon/faba_bean_ivip_and_genetic_cluster.csv", sep =",")

#datagen<-read.csv("../../../00.info/02.genetic_structure_jonathan_dijon/faba_bean_ivip_and_genetic_cluster.csv", sep =",")


# Select only the specific column to merge from data2
data2 <- datagen[, c(1,8)]
# Remove duplicates from data2
data2_unique <- data2[!duplicated(data2$ivip), ] 

# Merge the two datasets
df_long3 <- merge(df_long,data2_unique, by = "ivip", all.x = TRUE)

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







#################################################################################
#################################################################################
#########Including BBCH and block effect ########################################

#dataset formating to unsure compatibility for the merge ----
data_filtered_250_inf6$genotype<-as.factor(data_filtered_250_inf6$genotype)
data_filtered_xp006_inf6$genotype<-as.factor(data_filtered_xp006_inf6$genotype)
data_filtered_xp008_inf6$genotype<-as.factor(data_filtered_xp008_inf6$genotype)
data_filtered_xp010_inf6$genotype<-as.factor(data_filtered_xp010_inf6$genotype)

data_filtered_250_inf6$age_plant<-as.factor(data_filtered_250_inf6$age_plant)
data_filtered_xp006_inf6$age_plant<-as.factor(data_filtered_xp006_inf6$age_plant)
data_filtered_xp008_inf6$age_plant<-as.factor(data_filtered_xp008_inf6$age_plant)
data_filtered_xp010_inf6$age_plant<-as.factor(data_filtered_xp010_inf6$age_plant)

data_filtered_xp008_inf6$block<-as.factor(data_filtered_xp008_inf6$block)
data_filtered_xp010_inf6$block<-as.factor(data_filtered_xp010_inf6$block)

#add a block column in the big screening as a combinaison of tray number and session
data_filtered_250_inf6$block <- paste(data_filtered_250_inf6$tray_id, data_filtered_250_inf6$session, sep = "")
# in this dataset we now have a column session and a column block nested within the session

#P2XP006 is now the only dataset without block

# add a session column 
data_filtered_xp006_inf6$session <-"xp006"
data_filtered_xp008_inf6$session <-"xp008"
data_filtered_xp010_inf6$session <-"xp010"
# the big screening contained 7 session, it is therefore not equilibrated we 
#have each genotypes repeated 2 types per session from the big screen and 
# <15 genotypes per session for all the P2XPOO..

#binding the datasets ----
df_combined <- bind_rows(data_filtered_250_inf6, data_filtered_xp006_inf6, data_filtered_xp008_inf6, data_filtered_xp010_inf6)

#remove useless columns for clarity
df_combined <- df_combined[,c(1,3,7,12,14,15,20)]
# Annalyse with only the genotypes evaluated 3 or 4 time----
data_trimmed1 <- df_combined %>%
  filter(ivip %in% combined_data$ivip)

data_trimmed1 <- merge(data_trimmed1,data2_unique, by = "ivip", all.x = TRUE)
summary(data_trimmed1)
#ivip + block /session
mod<-lm(larvae_total_number ~ ivip+session / block,data=data_trimmed1)
anova(mod)


emm_genotype <- emmeans(mod, ~ ivip)
emm_df <- as.data.frame(emm_genotype)

emmeans_genetcluster1 <- merge(emm_df,data2_unique, by = "ivip", all.x = TRUE)


mod_anova <- aov(larvae_total_number ~ ivip + session/ block 
                 ,data = data_trimmed1)

tukey_groups2 <- HSD.test(mod_anova, "ivip", console = TRUE)
tukey_groups2<-as.data.frame(tukey_groups2$groups)

# IVIP alone

mod<-lm(larvae_total_number ~ ivip,data=data_trimmed1)
anova(mod)

emm_genotype <- emmeans(mod, ~ ivip)
emm_df <- as.data.frame(emm_genotype)


mod_anova <- aov(larvae_total_number ~ ivip, data = data_trimmed1)

tukey_groups2 <- HSD.test(mod_anova, "ivip", console = TRUE)
tukey_groups2<-as.data.frame(tukey_groups2$groups)

data_trimmed1$ivip <- droplevels(data_trimmed1$ivip)

ggplot(data_trimmed1, aes(x = ivip, y = larvae_total_number, fill = ivip)) +
  geom_boxplot() +
  labs(title = "Boxplot of Larvae Total Number by IVIP",
       x = "IVIP", y = "Larvae Total Number") +
  theme_minimal()


data_trimmed2 <- subset(data_trimmed1, ivip != "IVIP100388")
data_trimmed2 <- subset(data_trimmed2, ivip != "IVIP100396")
ggplot(data_trimmed2, aes(x = ivip, y = larvae_total_number, fill = cluster_name)) +
  geom_boxplot() +
  scale_fill_manual(values = c("red", "blue")) +  # Adjust colors as needed
  theme_minimal() +
  labs(title = "Boxplot of Larvae Total Number by IVIP",
       x = "IVIP", y = "Larvae Total Number", fill = "Cluster Name")


#annalyse P2XP008 ----

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
step_model <- step(lm(larvae_total_number~genotype+block+bbch_25days, data = data_filtered_xp008_inf6), direction = "both")
summary(step_model)
summary(data_filtered_xp008_inf6)

mod <- lm(larvae_total_number~genotype+block+bbch_25days, data = data_filtered_xp008_inf6)
anova(mod)

mod_anova <- aov(larvae_total_number ~ ivip + block + bbch_25days, 
                 data = data_filtered_xp008_inf6)

tukey_groups <- HSD.test(mod_anova, "ivip", console = TRUE)

#annalyse P2XP010 ----
# larvae number explained by other factor --> effect of the bbch when counting and growth bloc
mod <-lm(larvae_total_number~bbch_10days, data = data_filtered_xp010_inf6)
anova(mod)

mod <-lm(larvae_total_number~bbch_25days, data = data_filtered_xp010_inf6)
anova(mod)

mod <-lm(larvae_total_number~genotype, data = data_filtered_xp010_inf6)
anova(mod)

mod <-lm(larvae_total_number~block, data = data_filtered_xp010_inf6)
anova(mod)

mod <-lm(larvae_total_number~block+bbch_25days+genotype, data = data_filtered_xp010_inf6)
anova(mod)


# BBCH
mod <-lm(bbch_25days~block, data = data_filtered_xp010_inf6)
anova(mod)

mod <-lm(bbch_10days~genotype, data = data_filtered_xp010_inf6)
anova(mod)

mod <-lm(bbch_25days~genotype, data = data_filtered_xp010_inf6)
anova(mod)

#regression
step_model <- step(lm(larvae_total_number~genotype+block+bbch_25days+bbch_10days, data = data_filtered_xp010_inf6), direction = "both")
summary(step_model)


mod <- lm(larvae_total_number~ivip+block+bbch_25days, data = data_filtered_xp010_inf6)
anova(mod)

mod_anova <- aov(larvae_total_number ~ ivip + block + bbch_25days, 
                 data = data_filtered_xp010_inf6)

tukey_groups2 <- HSD.test(mod_anova, "ivip", console = TRUE)
tukey_groups2<-as.data.frame(tukey_groups2$groups)


