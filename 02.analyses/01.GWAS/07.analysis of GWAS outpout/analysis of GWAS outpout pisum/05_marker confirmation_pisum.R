library(data.table)
library(tidyverse)
library(dplyr)
library(multcompView)


#extract marker ----
myGD1 <-fread(file = "../../files_GWAS/PEAMUST_vfaba_F1_SNP_imputed_MAF5.012",header = FALSE, stringsAsFactors =  FALSE) #%>% select(-1)
myGD1 <- myGD1[, -1, with = FALSE]  # Removes the first column
myGD2 <-fread(file = "../../files_GWAS/PEAMUST_vfaba_F1_SNP_imputed_MAF5.012.indv",header = FALSE,stringsAsFactors =  FALSE)
# Load datasets
myGM <-fread(file = "../../files_GWAS/PEAMUST_vfaba_F1_SNP_imputed_MAF5.012.pos") %>% rename(
  "chromosome" = "V1",
  "position" = "V2"
)%>% mutate(snp_id=paste(chromosome,position, sep = "_")) %>% relocate(snp_id,.before = chromosome) 
colnames(myGD2)[1]<-"ivip"
colnames(myGD1)<- myGM$snp_id
myGD <-bind_cols(myGD2,myGD1)
extract_marker <- data.frame(myGD$ivip, myGD$chr2_1461736473, myGD$chr2_35683372, myGD$chr2_920208494, myGD$chr1S_1240520268)
colnames(extract_marker)[colnames(extract_marker) == "myGD.chr2_1461736473"] <- "chr2_1461736473"
colnames(extract_marker)[colnames(extract_marker) == "myGD.chr2_35683372"] <- "chr2_35683372"
colnames(extract_marker)[colnames(extract_marker) == "myGD.chr2_920208494"] <- "chr2_920208494"
colnames(extract_marker)[colnames(extract_marker) == "myGD.chr1S_1240520268"] <- "chr1S_1240520268"
colnames(extract_marker)[colnames(extract_marker) == "myGD.ivip"] <- "ivip"
extract_marker$chr2_1461736473<-as.factor(extract_marker$chr2_1461736473)
extract_marker$chr2_35683372<-as.factor(extract_marker$chr2_35683372)
extract_marker$chr2_920208494<-as.factor(extract_marker$chr2_920208494)

#data prep P2xp006 ----
data_xp006<-read.csv("01.data/P2-XP006_faba_vs_pisum_P123_fecundity-test/P2-XP006_larvae_fecundity_results.csv", stringsAsFactors = F)
data_xp006 <- data_xp006 %>% filter(!is.na(larvae_total_number))
data_filtered_xp006 <- data_xp006 %>% filter(adult_alive != "0")

#data prep P2xp008 ----
data_xp008<-read.csv("01.data/P2-XP008_faba_apisum_confirmation_2/P2-XP008_fecundity_data_apisum.csv", stringsAsFactors = F)
data_xp008 <- data_xp008 %>% filter(!is.na(larvae_total_number))
data_filtered_xp008 <- data_xp008 %>% filter(adult_alive != "0")

#data prep P2xp010 ----
data_xp010<-read.csv("01.data/P2-XP010_fabae_pisum_third_confirmation/P2XP010.csv", stringsAsFactors = F)
data_xp010 <- data_xp010 %>% filter(!is.na(larvae_total_number))
data_filtered_xp010 <- data_xp010 %>% filter(adult_alive != "0")

#data prep xp006 , xp008 and xp010 combined ----
# Define the columns 
columns_to_select <- c("ivip", "larvae_total_number")

# Create the subdataset by selecting those columns
subdataset1 <- data_filtered_xp008[, columns_to_select]
subdataset2 <- data_filtered_xp006[, columns_to_select]
subdataset3 <- data_filtered_xp010[, columns_to_select]

data_filtered_combined<-bind_rows(subdataset1,subdataset2,subdataset3)

combined_median <- data_filtered_combined %>%
  group_by(ivip) %>%
  summarise(Median_fecondity = median(larvae_total_number))

combined_median_marker <- merge(combined_median, extract_marker, by = "ivip")


combined_median_marker_no1 <- combined_median_marker %>% filter(chr2_920208494 != "1")
combined_median_marker_no1 <- combined_median_marker_no1 %>% filter(chr2_35683372 != "1")
combined_median_marker_no1 <- combined_median_marker_no1 %>% filter(chr2_1461736473 != "1")
combined_median_marker_no1 <- combined_median_marker_no1 %>% filter(chr1S_1240520268 != "1")


combined_median_marker_no1$chr2_920208494 <- droplevels(combined_median_marker_no1$chr2_920208494)
combined_median_marker_no1$chr2_35683372 <- droplevels(combined_median_marker_no1$chr2_35683372)
combined_median_marker_no1$chr2_1461736473 <- droplevels(combined_median_marker_no1$chr2_1461736473)
combined_median_marker_no1$chr2_1461736473 <- droplevels(combined_median_marker_no1$chr1S_1240520268)

#data prep big screening ----
data<-read.csv("01.data/2023_screen_specifics_251genotypes/01.acyrthosiphon_pisum_P123amp_on_faba_collection/acyrthospipon_pisum_fecundity_on_faba_collection.csv", stringsAsFactors = F)
data <- data %>% filter(!is.na(larvae_total_number))
data_filtered <- data %>% filter(adult_alive != "0")

median <- data_filtered %>%
  group_by(ivip) %>%
  summarise(Median_fecondity = median(larvae_total_number))
median_marker <- merge(median, extract_marker, by = "ivip")


median_marker_no1 <- median_marker %>% filter(chr2_920208494 != "1")
median_marker_no1 <- median_marker_no1 %>% filter(chr2_35683372 != "1")
median_marker_no1 <- median_marker_no1 %>% filter(chr2_1461736473 != "1")
median_marker_no1$chr2_920208494 <- droplevels(median_marker_no1$chr2_920208494)
median_marker_no1$chr2_35683372 <- droplevels(median_marker_no1$chr2_35683372)
median_marker_no1$chr2_1461736473 <- droplevels(median_marker_no1$chr2_1461736473)


#analyse chr2_920208494 ----

# Calculate the number of observations for each level of chr2_920208494
data_count <- combined_median_marker %>%
  group_by(chr2_920208494) %>%
  summarise(n = n())

# Create the boxplot with number of replicates displayed
combined_median_marker %>%
  mutate(chr2_920208494 = factor(chr2_920208494, levels = c(0, 1, 2))) %>%  
  ggplot(aes(x = chr2_920208494, y = Median_fecondity)) +  # Boxplot aesthetics
  geom_boxplot(fill = "skyblue", color = "black", outlier.shape = 16, outlier.size = 2) +  # Customize the boxplot
  geom_text(data = data_count, 
            aes(x = chr2_920208494, 
                y = max(combined_median_marker$Median_fecondity) + 0.5, 
                label = paste("n =", n)), 
            size = 5, vjust = 0) +  # Add number of observations text
  theme_minimal() +
  labs(
    title = "Larvae Number Depending on Marker chr2_920208494_1349197397 (A. pisum)",
    x = "Marker on chr2_920208494",
    y = "Larvae Total Number"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels for better readability

# Calculate the number of observations for each level of chr2_920208494
data_count <- median_marker %>%
  group_by(chr2_920208494) %>%
  summarise(n = n())

median_marker %>% 
  mutate(chr2_920208494 = factor(chr2_920208494, levels = c(0, 1, 2))) %>%  
  ggplot(aes(x = chr2_920208494, y = Median_fecondity)) +  # Boxplot aesthetics
  geom_boxplot(fill = "skyblue", color = "black", outlier.shape = 16, outlier.size = 2) +  # Customize the boxplot
  theme_minimal() +
  labs(
    title = "Larvae Number confirmation depending on marker chr2_920208494 (A. pisum)",
    x = "marker on chr2_920208494",
    y = "Larvae Total Number"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels for better readability


#analyse chr2_1461736473 ----
data_count <- combined_median_marker %>%
  group_by(chr2_1461736473) %>%
  summarise(n = n())

combined_median_marker %>% 
  mutate(chr2_1461736473 = factor(chr2_1461736473, levels = c(0, 1, 2))) %>%  
  ggplot(aes(x = chr2_1461736473, y = Median_fecondity)) +  # Boxplot aesthetics
  geom_boxplot(fill = "skyblue", color = "black", outlier.shape = 16, outlier.size = 2) +  # Customize the boxplot
  theme_minimal() +
  labs(
    title = "Larvae Number confirmation depending on marker chr2_1461736473 (A. pisum)",
    x = "marker on chr2_1461736473",
    y = "Larvae Total Number"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels for better readability

data_count <- median_marker %>%
  group_by(chr2_1461736473) %>%
  summarise(n = n())

median_marker %>% 
  mutate(chr2_1461736473 = factor(chr2_1461736473, levels = c(0, 1, 2))) %>%  
  ggplot(aes(x = chr2_1461736473, y = Median_fecondity)) +  # Boxplot aesthetics
  geom_boxplot(fill = "skyblue", color = "black", outlier.shape = 16, outlier.size = 2) +  # Customize the boxplot
  theme_minimal() +
  labs(
    title = "Larvae Number confirmation depending on marker chr2_1461736473 (A. pisum)",
    x = "marker on chr2_1461736473",
    y = "Larvae Total Number"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels for better readability


#analyse chr2_35683372 ----

data_count <- combined_median_marker %>%
  group_by(chr2_35683372) %>%
  summarise(n = n())

combined_median_marker %>% 
  mutate(chr2_35683372 = factor(chr2_35683372, levels = c(0, 1, 2))) %>%  
  ggplot(aes(x = chr2_35683372, y = Median_fecondity)) +  # Boxplot aesthetics
  geom_boxplot(fill = "skyblue", color = "black", outlier.shape = 16, outlier.size = 2) +  # Customize the boxplot
  theme_minimal() +
  labs(
    title = "Larvae Number confirmation depending on marker chr2_35683372 (A. pisum)",
    x = "marker on chr2_35683372",
    y = "Larvae Total Number"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels for better readability

data_count <- median_marker %>%
  group_by(chr2_35683372) %>%
  summarise(n = n())

median_marker %>% 
  mutate(chr2_35683372 = factor(chr2_35683372, levels = c(0, 1, 2))) %>%  
  ggplot(aes(x = chr2_35683372, y = Median_fecondity)) +  # Boxplot aesthetics
  geom_boxplot(fill = "skyblue", color = "black", outlier.shape = 16, outlier.size = 2) +  # Customize the boxplot
  theme_minimal() +
  labs(
    title = "Larvae Number confirmation depending on marker chr2_35683372(A. pisum)",
    x = "marker on chr2_35683372",
    y = "Larvae Total Number"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels for better readability


#analyse all QTL ----
#big screen ----
# Create a new variable for the interaction of the three markers
median_marker_no1$interaction <- interaction(median_marker_no1$chr2_920208494, 
                                                   median_marker_no1$chr2_35683372, 
                                                   median_marker_no1$chr2_1461736473)

# Calculate the mean larvae_total_number for each interaction level
interaction_means <- median_marker_no1 %>%
  group_by(interaction) %>%
  summarise(mean_larvae = median(Median_fecondity, na.rm = TRUE))

# Reorder the interaction factor based on mean larvae_total_number
ordered_interactions <- interaction_means %>%
  arrange(mean_larvae) %>%
  pull(interaction)

# Reorder the interaction factor in the original dataset
median_marker_no1$interaction <- factor(median_marker_no1$interaction, levels = ordered_interactions)

data_count <- median_marker_no1 %>%
  group_by(interaction) %>%
  summarise(n = n())

ggplot(median_marker_no1, aes(x = interaction, y = Median_fecondity, fill = interaction)) +
  geom_boxplot() +
  geom_text(data = data_count, 
            aes(x = interaction, 
                y = max(median_marker_no1$Median_fecondity) + 0.5, 
                label = paste("n =", n)), 
            size = 5, vjust = 0) +
  labs(x = "Interaction between chr2_920208494, chr2_35683372, and chr2_1461736473", 
       y = "Larvae Total Number", 
       title = "Interaction of chr2_920208494, chr2_35683372, and chr2_1461736473 on Larvae Total Number") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))



# Fit the linear model with the 'interaction' column
model <- lm(Median_fecondity ~ interaction, data = median_marker_no1)

# Perform ANOVA to check if the model is significant
anova_result <- anova(model)
print(anova_result)

# Perform Tukey HSD post-hoc test
tukey_result <- TukeyHSD(aov(model))

# View the Tukey HSD test results (for the interaction factor)
print(tukey_result$interaction)

# Use multcompView to assign group letters based on Tukey HSD results
# Extract the p-values for the interaction factor
group_letters <- multcompLetters(tukey_result$interaction[, "p adj"])

# Print the group labels
print(group_letters)


# confirmation----
# Create a new variable for the interaction of the three markers
combined_median_marker_no1$interaction <- interaction(combined_median_marker_no1$chr2_920208494, 
                                                               combined_median_marker_no1$chr2_35683372, 
                                                               combined_median_marker_no1$chr2_1461736473)

# Calculate the mean larvae_total_number for each interaction level
interaction_means <- combined_median_marker_no1 %>%
  group_by(interaction) %>%
  summarise(mean_larvae = median(Median_fecondity, na.rm = TRUE))

# Reorder the interaction factor based on mean larvae_total_number
ordered_interactions <- interaction_means %>%
  arrange(mean_larvae) %>%
  pull(interaction)

# Reorder the interaction factor in the original dataset
combined_median_marker_no1$interaction <- factor(combined_median_marker_no1$interaction, levels = ordered_interactions)

data_count <- combined_median_marker_no1 %>%
  group_by(interaction) %>%
  summarise(n = n())

ggplot(combined_median_marker_no1, aes(x = interaction, y = Median_fecondity, fill = interaction)) +
  geom_boxplot() +
  geom_text(data = data_count, 
            aes(x = interaction, 
                y = max(combined_median_marker_no1$Median_fecondity) + 0.5, 
                label = paste("n =", n)), 
            size = 5, vjust = 0) +
  labs(x = "Interaction between chr2_920208494, chr2_35683372, and chr2_1461736473", 
       y = "Larvae Total Number", 
       title = "Interaction of chr2_920208494, chr2_35683372, and chr2_1461736473 on Larvae Total Number (Ordered by Mean)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Fit the linear model with the 'interaction' column
model <- lm(Median_fecondity ~ interaction, data = combined_median_marker_no1)

# Perform ANOVA to check if the model is significant
anova(model)


