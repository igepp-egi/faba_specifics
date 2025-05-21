#extract marker ----
myGD1 <-fread(file = "../../files GWAS/PEAMUST_vfaba_F1_SNP_imputed_MAF5.012",header = FALSE, stringsAsFactors =  FALSE) #%>% select(-1)
myGD1 <- myGD1[, -1, with = FALSE]  # Removes the first column
myGD2 <-fread(file = "../../files GWAS/PEAMUST_vfaba_F1_SNP_imputed_MAF5.012.indv",header = FALSE,stringsAsFactors =  FALSE)
# Load datasets
myGM <-fread(file = "../../files GWAS/PEAMUST_vfaba_F1_SNP_imputed_MAF5.012.pos") %>% rename(
  "chromosome" = "V1",
  "position" = "V2"
)%>% mutate(snp_id=paste(chromosome,position, sep = "_")) %>% relocate(snp_id,.before = chromosome) 
colnames(myGD2)[1]<-"ivip"
colnames(myGD1)<- myGM$snp_id
myGD <-bind_cols(myGD2,myGD1)
extract_marker <- data.frame(myGD$ivip, myGD$chr5_1349197397, myGD$chr3_641834492, myGD$chr1L_1430504409)
colnames(extract_marker)[colnames(extract_marker) == "myGD.chr5_1349197397"] <- "chr5"
colnames(extract_marker)[colnames(extract_marker) == "myGD.chr3_641834492"] <- "chr3"
colnames(extract_marker)[colnames(extract_marker) == "myGD.chr1L_1430504409"] <- "chr1L"
colnames(extract_marker)[colnames(extract_marker) == "myGD.ivip"] <- "ivip"
extract_marker$chr1L<-as.factor(extract_marker$chr1L)
extract_marker$chr3<-as.factor(extract_marker$chr3)
extract_marker$chr5<-as.factor(extract_marker$chr5)

#data prep cluster infos ----
geneticcluster<- read.csv("00.info/02.genetic_structure_jonathan_dijon/faba_bean_ivip_and_genetic_cluster.csv", 
                                                           sep =",", 
                                                           stringsAsFactors = F)


#data prep P2xp007 ----
data_xp007<-read.csv("01.data/P2-XP007_fabae_AF001_confirmation/P2-XP007_confirmation_Afabae_10vs10genotypes.csv", stringsAsFactors = F)
data_xp007 <- data_xp007 %>% filter(!is.na(larvae_total_number))
data_filtered_xp007 <- data_xp007 %>% filter(adult_alive != "0")

#data prep P2xp009 ----
data_xp009<-read.csv("01.data/P2-XP009_afabae_third_confirmation/P2-XP009_fecundity_data_afabae.csv", stringsAsFactors = F)
data_xp009 <- data_xp009 %>% filter(!is.na(larvae_total_number))
data_filtered_xp009 <- data_xp009 %>% filter(adult_alive != "0")
#data prep xp007 and XP009 combined ----
# Define the columns 
columns_to_select <- c("ivip", "larvae_total_number")

# Create the subdataset by selecting those columns
subdataset1 <- data_filtered_xp009[, columns_to_select]
subdataset2 <- data_filtered_xp007[, columns_to_select]
data_filtered_xp007_xp009<-bind_rows(subdataset1,subdataset2)

xp007_xp009_median <- data_filtered_xp007_xp009 %>%
  group_by(ivip) %>%
  summarise(Median_fecondity = median(larvae_total_number))

xp007_xp009_median_marker <- merge(xp007_xp009_median, extract_marker, by = "ivip")

median_xp007_xp009_marker_geneticcluster <- geneticcluster %>% dplyr::select(ivip, cluster_name) %>% 
  distinct() %>%
  inner_join(., xp007_xp009_median_marker, by = "ivip")

median_xp007_xp009_marker_south <- subset(median_xp007_xp009_marker_geneticcluster, cluster_name %in% c("South", "South_America"))

median_xp007_xp009_marker_south_no1 <- median_xp007_xp009_marker_south %>% filter(chr5 != "1")
median_xp007_xp009_marker_south_no1 <- median_xp007_xp009_marker_south_no1 %>% filter(chr3 != "1")
median_xp007_xp009_marker_south_no1 <- median_xp007_xp009_marker_south_no1 %>% filter(chr1L != "1")
median_xp007_xp009_marker_south_no1$chr5 <- droplevels(median_xp007_xp009_marker_south_no1$chr5)
median_xp007_xp009_marker_south_no1$chr3 <- droplevels(median_xp007_xp009_marker_south_no1$chr3)






#data prep big screening ----
data<-read.csv("01.data/2023_screen_specifics_251genotypes/02.aphis_fabae_on_faba_collection/aphis_fabae_fecundity_on_faba_collection_2024-04-02.csv", stringsAsFactors = F)
data <- data %>% filter(!is.na(larvae_total_number))
data_filtered <- data %>% filter(adult_alive != "0")

median <- data_filtered %>%
  group_by(ivip) %>%
  summarise(Median_fecondity = median(larvae_total_number))
median_marker <- merge(median, extract_marker, by = "ivip")

median_marker_geneticcluster <- geneticcluster %>% dplyr::select(ivip, cluster_name) %>% 
  distinct() %>%
  inner_join(., median_marker, by = "ivip")

median_marker_south <- subset(median_marker_geneticcluster, cluster_name %in% c("South", "South_America"))

median_marker_south_no1 <- median_marker_south %>% filter(chr5 != "1")
median_marker_south_no1 <- median_marker_south_no1 %>% filter(chr3 != "1")
median_marker_south_no1 <- median_marker_south_no1 %>% filter(chr1L != "1")
median_marker_south_no1$chr5 <- droplevels(median_marker_south_no1$chr5)
median_marker_south_no1$chr3 <- droplevels(median_marker_south_no1$chr3)




#


data_filtered_marker_south <- subset(data_filtered_xp007_xp009_marker_geneticcluster, cluster_name %in% c("South", "South_America"))


#analyse chr5 (South only) ----

# Calculate the number of observations for each level of chr5
data_count <- median_xp007_xp009_marker_south %>%
  group_by(chr5) %>%
  summarise(n = n())

# Create the boxplot with number of replicates displayed
median_xp007_xp009_marker_south %>%
  mutate(chr5 = factor(chr5, levels = c(0, 1, 2))) %>%  
  ggplot(aes(x = chr5, y = Median_fecondity)) +  # Boxplot aesthetics
  geom_boxplot(fill = "skyblue", color = "black", outlier.shape = 16, outlier.size = 2) +  # Customize the boxplot
  geom_text(data = data_count, 
            aes(x = chr5, 
                y = max(median_xp007_xp009_marker_south$Median_fecondity) + 0.5, 
                label = paste("n =", n)), 
            size = 5, vjust = 0) +  # Add number of observations text
  theme_minimal() +
  labs(
    title = "Larvae Number Confirmation Depending on Marker chr5_1349197397 in South Cluster (A. fabae)",
    x = "Marker on chr5",
    y = "Larvae Total Number"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels for better readability

# Calculate the number of observations for each level of chr5
data_count <- median_marker_south %>%
  group_by(chr5) %>%
  summarise(n = n())

median_marker_south %>% 
  mutate(chr5 = factor(chr5, levels = c(0, 1, 2))) %>%  
  ggplot(aes(x = chr5, y = Median_fecondity)) +  # Boxplot aesthetics
  geom_boxplot(fill = "skyblue", color = "black", outlier.shape = 16, outlier.size = 2) +  # Customize the boxplot
  theme_minimal() +
  labs(
    title = "Larvae Number confirmation depending on marker chr5 in south cluster (A. fabae)",
    x = "marker on chr5",
    y = "Larvae Total Number"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels for better readability


#analyse chr1L (South only) ----
data_count <- median_xp007_xp009_marker_south %>%
  group_by(chr1L) %>%
  summarise(n = n())

median_xp007_xp009_marker_south %>% 
  mutate(chr1L = factor(chr1L, levels = c(0, 1, 2))) %>%  
  ggplot(aes(x = chr1L, y = Median_fecondity)) +  # Boxplot aesthetics
  geom_boxplot(fill = "skyblue", color = "black", outlier.shape = 16, outlier.size = 2) +  # Customize the boxplot
  theme_minimal() +
  labs(
    title = "Larvae Number confirmation depending on marker chr1L in south cluster (A. fabae)",
    x = "marker on chr1L",
    y = "Larvae Total Number"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels for better readability

data_count <- median_marker_south %>%
  group_by(chr1L) %>%
  summarise(n = n())

median_marker_south %>% 
  mutate(chr1L = factor(chr1L, levels = c(0, 1, 2))) %>%  
  ggplot(aes(x = chr1L, y = Median_fecondity)) +  # Boxplot aesthetics
  geom_boxplot(fill = "skyblue", color = "black", outlier.shape = 16, outlier.size = 2) +  # Customize the boxplot
  theme_minimal() +
  labs(
    title = "Larvae Number confirmation depending on marker chr1L in south cluster (A. fabae)",
    x = "marker on chr1L",
    y = "Larvae Total Number"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels for better readability


#analyse chr3 ----

data_count <- xp007_xp009_median_marker %>%
  group_by(chr3) %>%
  summarise(n = n())

xp007_xp009_median_marker %>% 
  mutate(chr3 = factor(chr3, levels = c(0, 1, 2))) %>%  
  ggplot(aes(x = chr3, y = Median_fecondity)) +  # Boxplot aesthetics
  geom_boxplot(fill = "skyblue", color = "black", outlier.shape = 16, outlier.size = 2) +  # Customize the boxplot
  theme_minimal() +
  labs(
    title = "Larvae Number confirmation depending on marker chr3 in south cluster (A. fabae)",
    x = "marker on chr3",
    y = "Larvae Total Number"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels for better readability

data_count <- median_marker %>%
  group_by(chr3) %>%
  summarise(n = n())

median_marker %>% 
  mutate(chr3 = factor(chr3, levels = c(0, 1, 2))) %>%  
  ggplot(aes(x = chr3, y = Median_fecondity)) +  # Boxplot aesthetics
  geom_boxplot(fill = "skyblue", color = "black", outlier.shape = 16, outlier.size = 2) +  # Customize the boxplot
  theme_minimal() +
  labs(
    title = "Larvae Number confirmation depending on marker chr3 in south cluster (A. fabae)",
    x = "marker on chr3",
    y = "Larvae Total Number"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels for better readability


#analyse all QTL (South only) ----
#big screen ----
# Create a new variable for the interaction of the three markers
median_marker_south_no1$interaction <- interaction(median_marker_south_no1$chr5, 
                                                   median_marker_south_no1$chr3, 
                                                   median_marker_south_no1$chr1L)

# Calculate the mean larvae_total_number for each interaction level
interaction_means <- median_marker_south_no1 %>%
  group_by(interaction) %>%
  summarise(mean_larvae = median(Median_fecondity, na.rm = TRUE))

# Reorder the interaction factor based on mean larvae_total_number
ordered_interactions <- interaction_means %>%
  arrange(mean_larvae) %>%
  pull(interaction)

# Reorder the interaction factor in the original dataset
median_marker_south_no1$interaction <- factor(median_marker_south_no1$interaction, levels = ordered_interactions)

data_count <- median_marker_south_no1 %>%
  group_by(interaction) %>%
  summarise(n = n())

ggplot(median_marker_south_no1, aes(x = interaction, y = Median_fecondity, fill = interaction)) +
  geom_boxplot() +
  geom_text(data = data_count, 
            aes(x = interaction, 
                y = max(median_marker_south_no1$Median_fecondity) + 0.5, 
                label = paste("n =", n)), 
            size = 5, vjust = 0) +
  labs(x = "Interaction between chr5, chr3, and chr1L", 
       y = "Larvae Total Number", 
       title = "Interaction of chr5, chr3, and chr1L on Larvae Total Number") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))



# Fit the linear model with the 'interaction' column
model <- lm(Median_fecondity ~ interaction, data = median_marker_south_no1)

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
median_xp007_xp009_marker_south_no1$interaction <- interaction(median_xp007_xp009_marker_south_no1$chr5, 
                                                               median_xp007_xp009_marker_south_no1$chr3, 
                                                               median_xp007_xp009_marker_south_no1$chr1L)

# Calculate the mean larvae_total_number for each interaction level
interaction_means <- median_xp007_xp009_marker_south_no1 %>%
  group_by(interaction) %>%
  summarise(mean_larvae = median(Median_fecondity, na.rm = TRUE))

# Reorder the interaction factor based on mean larvae_total_number
ordered_interactions <- interaction_means %>%
  arrange(mean_larvae) %>%
  pull(interaction)

# Reorder the interaction factor in the original dataset
median_xp007_xp009_marker_south_no1$interaction <- factor(median_xp007_xp009_marker_south_no1$interaction, levels = ordered_interactions)

data_count <- median_xp007_xp009_marker_south_no1 %>%
  group_by(interaction) %>%
  summarise(n = n())

ggplot(median_xp007_xp009_marker_south_no1, aes(x = interaction, y = Median_fecondity, fill = interaction)) +
  geom_boxplot() +
  geom_text(data = data_count, 
            aes(x = interaction, 
                y = max(median_xp007_xp009_marker_south_no1$Median_fecondity) + 0.5, 
                label = paste("n =", n)), 
            size = 5, vjust = 0) +
  labs(x = "Interaction between chr5, chr3, and chr1L", 
       y = "Larvae Total Number", 
       title = "Interaction of chr5, chr3, and chr1L on Larvae Total Number (Ordered by Mean)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Fit the linear model with the 'interaction' column
model <- lm(Median_fecondity ~ interaction, data = median_xp007_xp009_marker_south_no1)

# Perform ANOVA to check if the model is significant
anova(model)


