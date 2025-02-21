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


# Import all fecundity measurements from 2022-2023 (N=251 genotypes), P2-XP006 (N=16 genotypes), P2-XP008 (N=7 geno.) and P2-XP010 (N=7 geno.)
source("02.analyses/02.Discriminant_Genotypes/2025_consistency_analysis/acyrthosiphon_pisum_P123/prepare_data_from_bigscreen_P2XP006_P2XP008_P2XP010.R")


####### ploting median consistency ##########################################
# plot consistency of all genotypes----

# Plotting with ggplot2 
ggplot(df_long, aes(x = xp, y = median_larvae_number, group = ivip, color = ivip)) +
  geom_line() +              # Add lines connecting the same genotype across years
  geom_point() +             # Add points for each year/genotype
  labs(x = "expe", y = "Larvae number", title = "confirmation of the genotypes A.pisum") +
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
  labs(x = "Experiment", y = "Median larvae number", title = "Genotype consistency (A. pisum fecundity test)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  # Rotate x-axis labels for better readability
  scale_color_manual(values = c("red", "blue","pink","orange", "green", "yellow"))        


#################################################################################
######### Including BBCH and block effect ########################################
#################################################################################

# add a XP column identifier
data_filtered_250_inf6$xp_id = "big_screen"
data_filtered_xp006_inf6$xp_id = "P2XP006"
data_filtered_xp008_inf6$xp_id = "P2XP008"
data_filtered_xp010_inf6$xp_id = "P2XP010"


data_filtered_xp006_inf6$session <-"xp006"
data_filtered_xp008_inf6$session <-"xp008"
data_filtered_xp010_inf6$session <-"xp010"
# the big screening contained 7 session, it is therefore not equilibrated we 
#have each genotypes repeated 2 types per session from the big screen and 
# <15 genotypes per session for all the P2XPOO..

# binding the datasets ----
df_combined <- bind_rows(data_filtered_250_inf6, 
                         data_filtered_xp006_inf6, 
                         data_filtered_xp008_inf6, 
                         data_filtered_xp010_inf6)

# Count number of xp for each ivip number
# Keep if the ivip number has been evaluated at least 3 times
ivip_to_keep_n3 <- df_combined %>%
  group_by(ivip) %>%
  summarise(n = n_distinct(xp_id)) %>% 
  filter(n >= 3)

#remove useless columns for clarity
df_combined_subset <- df_combined %>% 
  dplyr::select(genotype, ivip, block, session, xp_id, larvae_total_number) %>% 
  filter(ivip %in% ivip_to_keep_n3$ivip) # Analyse with only the genotypes evaluated 3 or 4 times----

# Import genetic cluster info
genetic_clusters <- read.csv("00.info/02.genetic_structure_jonathan_dijon/faba_bean_ivip_and_genetic_cluster.csv", 
         sep =",", 
         stringsAsFactors = F) %>% 
  dplyr::select(ivip, cluster_name)

df_combined_subset_with_genetic_clusters <- inner_join(df_combined_subset, 
                                                       genetic_clusters, 
                                                       by = "ivip") 

# Boxplots: all 6 genotypes selected through all experiments

df_combined_subset_with_genetic_clusters %>% 
  arrange() %>% 
  ggplot(., 
       aes(x = ivip, y = larvae_total_number, fill = cluster_name)) +
  geom_boxplot() +
  geom_jitter(width = 0.05, alpha = 0.4, aes(shape = xp_id)) +
  scale_fill_manual(values = c("red", "blue")) +  # Adjust colors as needed
  theme_minimal() +
  labs(title = "Larvae Total Number by IVIP (A. pisum P123)",
       x = "IVIP", y = "Larvae Total Number", fill = "Cluster Name") 

# Boxplots: only consistent genotypes selected through all experiments, excluding IVIP100388 and IVIP100396
inconsistent_genotypes <- c("IVIP100388", "IVIP100396")

p_4genotypes <- df_combined_subset_with_genetic_clusters %>% 
  filter(!ivip %in% inconsistent_genotypes) %>%
  ggplot(aes(x = ivip, y = larvae_total_number, fill = cluster_name)) +
  geom_boxplot() +
  geom_jitter(width = 0.05, alpha = 0.4, aes(shape = xp_id)) +
  scale_fill_manual(values = c("red", "blue")) +  # Adjust colors as needed
  theme_minimal() +
  labs(title = "Larvae Total Number per IVIP (A. pisum P123)",
       x = "IVIP", y = "Larvae Total Number", fill = "Cluster Name");p_4genotypes

# Count number of observations per ivip and add to the graph
n_obs_per_ivip <- df_combined_subset_with_genetic_clusters %>% 
  group_by(ivip) %>% 
  summarise(n = n()) %>% 
  arrange(desc(n)) %>% 
  inner_join(., genetic_clusters, by = "ivip") %>% 
  filter(!ivip %in% inconsistent_genotypes)

# Add number of observation to the p_4genotypes graph
p_4genotypes + 
  geom_text(data = n_obs_per_ivip, aes(x = ivip, y = 70, label = paste("N = ",n)), vjust = -0.5, hjust = 0.5, size = 4)

data_trimmed1 <- df_combined %>%
  filter(ivip %in% combined_data$ivip)

data_trimmed1 <- merge(data_trimmed1,infos_cluster, by = "ivip", all.x = TRUE)
summary(data_trimmed1)

#ivip + block /session
mod<-lm(larvae_total_number ~ ivip+session / block,data=data_trimmed1)
anova(mod)


emm_genotype <- emmeans(mod, ~ ivip)
emm_df <- as.data.frame(emm_genotype)

emmeans_genetcluster1 <- merge(emm_df,infos_cluster, by = "ivip", all.x = TRUE)


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


