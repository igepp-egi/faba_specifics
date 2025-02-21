#packages ----
suppressPackageStartupMessages(library(corrplot))
suppressPackageStartupMessages(library(Hmisc))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(agricolae))


#data prep 250 ----
data250<-read.csv("../../faba_specifics/01.data/2023_screen_specifics_251genotypes/02.aphis_fabae_on_faba_collection/aphis_fabae_fecundity_on_faba_collection_2024-04-02.csv", sep =",")


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

data_XP007<-read.csv("../01.data/P2-XP007_fabae_AF001_confirmation/P2-XP007_confirmation_Afabae_10vs10genotypes.csv", sep =",")
data_XP007 <- data_XP007 %>% filter(!is.na(larvae_total_number))
data_filtered_XP007 <- data_XP007 %>% filter(adult_alive != "0")

data_filtered_XP007$ivip<-as.factor(data_filtered_XP007$ivip)
summary(data_filtered_XP007)
data_filtered_XP007_inf6 <- data_filtered_XP007 %>%
  group_by(ivip) %>%           # Group by genotype
  filter(n() >= 6) %>%             # Keep only groups with 6 or more repetitions
  ungroup()                        # Ungroup to return to original data structure
summary(data_filtered_XP007_inf6)

data_filtered_XP007_inf6
#### data prep P2XP009 ----

data_XP009<-read.csv("../01.data/P2-XP009_afabae_third_confirmation/P2-XP009_fecundity_data_afabae.csv", sep =",")
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


data_median_eXPe250 <- data_filtered_250_inf6 %>%
  group_by(ivip) %>%
  summarise(Median_fecondity = median(larvae_total_number))
colnames(data_median_eXPe250)[colnames(data_median_eXPe250) == "Median_fecondity"] <- "eXPe1"
#median extract P2XP007 ----

data_median_XP007 <- data_filtered_XP007_inf6 %>%
  group_by(ivip) %>%
  summarise(Median_fecondity = median(larvae_total_number))
colnames(data_median_XP007)[colnames(data_median_XP007) == "Median_fecondity"] <- "eXPe2"

#median extract P2XP009 ----

data_median_XP009 <- data_filtered_XP009_inf6 %>%
  group_by(ivip) %>%
  summarise(Median_fecondity = median(larvae_total_number))
colnames(data_median_XP009)[colnames(data_median_XP009) == "Median_fecondity"] <- "eXPe3"


# merging the newly created file: 250 genotypes screening, P2XP007, and P2XP009 ----
# Merge the data frames by Genotype, ensuring all genotypes are kept

combined_data <- merge(data_median_eXPe250, data_median_XP007, by = "ivip", all = TRUE)
combined_data <- merge(combined_data, data_median_XP009, by = "ivip", all = TRUE)


combined_data <- combined_data %>%
  rowwise() %>%
  filter(sum(!is.na(c_across(starts_with("eXPe")))) > 2) %>%
  ungroup()




colnames(combined_data) <- trimws(tolower(colnames(combined_data)))  # Trim spaces and make lowercase


# Reshape the data into long format
df_long <- combined_data %>%
  pivot_longer(cols = starts_with("eXPe"), 
               names_to = "test", 
               names_prefix = "eXPe",
               values_to = "eXPe")

# View the reshaped data
print(df_long, n=250)







# Plotting with ggplot2
ggplot(df_long, aes(x = test, y = eXPe, group = ivip, color = ivip)) +
  geom_line() +              # Add lines connecting the same genotype across years
  geom_point() +             # Add points for each year/genotype
  labs(x = "eXPe", y = "Larvae number", title = "confirmation of the genotypes A.pissum") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels for better readability

# adding the genetic background----

datagen<-read.csv("../../faba_specifics/00.info/02.genetic_structure_jonathan_dijon/faba_bean_ivip_and_genetic_cluster.csv", sep =",")

# Select only the specific column you want to merge from data2
data2 <- datagen[, c(1,8)]
# Remove duplicates from data2
data2_unique <- data2[!duplicated(data2$ivip), ] 

# Merge the two datasets
df_long2 <- merge(df_long,data2_unique, by = "ivip", all.x = TRUE)


ggplot(df_long2, aes(x = test, y = eXPe, group = ivip, color = cluster_name)) +
  geom_line() +              # Add lines connecting the same genotype across years
  geom_point() +             # Add points for each year/genotype
  geom_label(aes(label = ivip), vjust = -0.5, hjust = 0.5, size = 2) + # Add genotype labels
  labs(x = "eXPe", y = "Larvae number", title = "Confirmation of the genotypes, A.fabae") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  # Rotate x-axis labels for better readability
  scale_color_manual(values = c("red", "blue","pink","orange", "green", "yellow"))        













#################################################################################
#################################################################################

#annalyse P2XP009 ----
data_filtered_XP009_inf6$block<-as.factor(data_filtered_XP009_inf6$block)
data_filtered_XP009_inf6$bbch_10days<-as.integer(data_filtered_XP009_inf6$bbch_10days)
data_filtered_XP009_inf6$bbch_25days<-as.integer(data_filtered_XP009_inf6$bbch_25days)
data_filtered_XP009_inf6$genotype<-as.factor(data_filtered_XP009_inf6$genotype)


# larvae number eXPlained by other factor --> effect of the bbch when counting and growth block
mod <-lm(larvae_total_number~bbch_10days, data = data_filtered_XP009_inf6)
anova(mod)

mod <-lm(larvae_total_number~bbch_25days, data = data_filtered_XP009_inf6)
anova(mod)

mod <-lm(larvae_total_number~genotype, data = data_filtered_XP009_inf6)
anova(mod)

mod <-lm(larvae_total_number~block, data = data_filtered_XP009_inf6)
anova(mod)

mod <-lm(larvae_total_number~block+bbch_25days+genotype, data = data_filtered_XP009_inf6)
anova(mod)


summary(data_filtered_XP009_inf6)

# BBCH
mod <-lm(bbch_25days~block, data = data_filtered_XP009_inf6)
anova(mod)

mod <-lm(bbch_10days~genotype, data = data_filtered_XP009_inf6)
anova(mod)

mod <-lm(bbch_25days~genotype, data = data_filtered_XP009_inf6)
anova(mod)

#regression
step_model <- step(lm(larvae_total_number~ivip+block+bbch_25days, data = data_filtered_XP009_inf6), direction = "both")
summary(step_model)


mod <- lm(larvae_total_number~ivip+block+bbch_25days, data = data_filtered_XP009_inf6)
anova(mod)

mod_anova <- aov(larvae_total_number ~ ivip + block, 
                 data = data_filtered_XP009_inf6)

tukey_groups2 <- HSD.test(mod_anova, "ivip", console = TRUE)
tukey_groups2<-as.data.frame(tukey_groups2$groups)








#dataset formating to unsure compatibility for the merge ----



data_filtered_250_inf6$genotype<-as.factor(data_filtered_250_inf6$genotype)
data_filtered_XP007_inf6$genotype<-as.factor(data_filtered_XP007_inf6$genotype)
data_filtered_XP009_inf6$genotype<-as.factor(data_filtered_XP009_inf6$genotype)

data_filtered_250_inf6$age_plant<-as.factor(data_filtered_250_inf6$age_plant)
data_filtered_XP007_inf6$age_plant<-as.factor(data_filtered_XP007_inf6$age_plant)
data_filtered_XP009_inf6$age_plant<-as.factor(data_filtered_XP009_inf6$age_plant)

data_filtered_XP007_inf6$adult_alive<-as.integer(data_filtered_XP007_inf6$adult_alive)
data_filtered_XP009_inf6$adult_alive<-as.integer(data_filtered_XP009_inf6$adult_alive)
data_filtered_250_inf6$adult_alive<-as.integer(data_filtered_250_inf6$adult_alive)

data_filtered_XP007_inf6$adult_winged<-as.integer(data_filtered_XP007_inf6$adult_winged)
data_filtered_XP009_inf6$adult_winged<-as.integer(data_filtered_XP009_inf6$adult_winged)
data_filtered_250_inf6$adult_winged<-as.integer(data_filtered_250_inf6$adult_winged)

data_filtered_XP009_inf6$block<-as.factor(data_filtered_XP009_inf6$block)

#add a block column in the big screening as a combinaison of tray number and session
data_filtered_250_inf6$block <- paste(data_filtered_250_inf6$tray_id, data_filtered_250_inf6$session, sep = "")
# in this dataset we now have a column session and a column block nested within the session

#P2XP007 is now the only dataset without block

# add a session column 
data_filtered_XP007_inf6$session <-"XP007"
data_filtered_XP009_inf6$session <-"XP009"

# the big screening contained 7 session, it is therefore not equilibrated we 
#have each genotypes repeated 2 types per session from the big screen and 
# <15 genotypes per session for all the P2XPOO..

#binding the datasets ----
df_combined <- bind_rows(data_filtered_250_inf6, data_filtered_XP007_inf6, data_filtered_XP009_inf6)

#remove useless columns for clarity
df_combined <- dplyr::select(df_combined, c(genotype,ivip,block,session,larvae_total_number))

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



n_obs_per_ivip <- data_trimmed1 %>% 
  group_by(ivip) %>% 
  summarise(n = n()) %>% 
  arrange(desc(n))
 
