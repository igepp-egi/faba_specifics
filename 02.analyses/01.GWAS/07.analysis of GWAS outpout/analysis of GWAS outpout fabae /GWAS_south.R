#package----
library(data.table)
library(vcfR)   
library(tidyverse)
library(dplyr)

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
extract_marker <- data.frame(myGD$ivip, myGD$chr5_1349197397, myGD$chr1L_1430504409, myGD$chr3_641834492)
colnames(extract_marker)[colnames(extract_marker) == "myGD.chr5_1349197397"] <- "chr5"
colnames(extract_marker)[colnames(extract_marker) == "myGD.chr1L_1430504409"] <- "chr1L"
colnames(extract_marker)[colnames(extract_marker) == "myGD.chr3_641834492"] <- "chr3"

colnames(extract_marker)[colnames(extract_marker) == "myGD.ivip"] <- "ivip"
extract_marker$chr1L<-as.factor(extract_marker$chr1L)
extract_marker$chr5<-as.factor(extract_marker$chr5)
extract_marker$chr3<-as.factor(extract_marker$chr3)


#data prep cluster infos ----
geneticcluster<- read.csv("00.info/02.genetic_structure_jonathan_dijon/faba_bean_ivip_and_genetic_cluster.csv", 
                          sep =",", 
                          stringsAsFactors = F)


#data prep big screening ----
datafabae<-read.csv("02.analyses/01.GWAS/01_phenotipic_analysis/02.ajusted_means/fabae/emmeans_fabae.csv", sep=",", stringsAsFactors = FALSE)

data_marker <- merge(datafabae, extract_marker, by = "ivip")

data_marker_geneticcluster <- geneticcluster %>% dplyr::select(ivip, cluster_name) %>% 
  distinct() %>%
  inner_join(., data_marker, by = "ivip")

data_marker_south <- subset(data_marker_geneticcluster, cluster_name %in% c("South", "South_America"))
data_marker_north <- subset(data_marker_geneticcluster, cluster_name %in% c("North"))
data_marker_east <- subset(data_marker_geneticcluster, cluster_name %in% c("East"))

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






# generic plot----

data_marker_geneticcluster %>% 
  ggplot(aes(x = cluster_name, y = emmean)) + 
  geom_boxplot(fill = "skyblue", color = "black", outlier.shape = NA) +  # Hide outliers
  geom_jitter(width = 0.2, alpha = 0.5, color = "blue") +  # Only your points will show
  theme_minimal() +
  labs(
    title = "Larvae Number depending on the cluster (A. fabae)",
    x = "Genetic cluster",
    y = "Larvae Total Number"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


#distribution chr5 among clusters (big screening)----





# Count chr5 values by cluster
summary_df <- data_marker_geneticcluster %>%
  group_by(cluster_name, chr5) %>%
  summarise(count = n(), .groups = "drop") %>%
  filter(chr5 %in% c(0, 1, 2)) %>%
  mutate(chr5 = as.factor(chr5))

# Ensure all combinations are present
summary_df_complete <- summary_df %>%
  complete(cluster_name, chr5 = factor(0:2), fill = list(count = 0))

# Add "Total" column
total_counts <- summary_df_complete %>%
  group_by(chr5) %>%
  summarise(count = sum(count)) %>%
  mutate(cluster_name = "Total")

# Combine and set cluster_name order (Total last)
final_df <- bind_rows(summary_df_complete, total_counts)
final_df$cluster_name <- factor(final_df$cluster_name, 
                                levels = c(setdiff(unique(final_df$cluster_name), "Total"), "Total"))

# Plot
ggplot(final_df, aes(x = cluster_name, y = count, fill = chr5)) +
  geom_bar(stat = "identity") +
  labs(
    title = "Counts of 0s, 1s, and 2s at chr5 by Cluster",
    x = "Cluster Name",
    y = "Count",
    fill = "marker"
  ) +
  theme_minimal()

data_count <- data_marker_south %>%
  group_by(chr5) %>%
  summarise(n = n())

data_count <- data_marker_geneticcluster %>%
  group_by(chr5) %>%
  summarise(n = n())


summary_df <- data_marker_geneticcluster %>%
  group_by(cluster_name) %>%
  summarise(num_twos = sum(chr5 == 2))




#analyse chr5 (big screening) ----

# Calculate the number of observations for each level of chr5


data_marker_south %>% 
  mutate(chr5 = factor(chr5, levels = c(0, 1, 2))) %>%  
  ggplot(aes(x = chr5, y = emmean)) +  # Boxplot aesthetics
  geom_boxplot(fill = "skyblue", color = "black", outlier.shape = 16, outlier.size = 2) +  # Customize the boxplot
  theme_minimal() +
  labs(
    title = "Larvae Number confirmation depending on marker chr5 in south cluster (A. fabae)",
    x = "marker on chr5",
    y = "Larvae Total Number"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels for better readability




data_marker_north %>% 
  mutate(chr5 = factor(chr5, levels = c(0, 1, 2))) %>%  
  ggplot(aes(x = chr5, y = emmean)) +  # Boxplot aesthetics
  geom_boxplot(fill = "skyblue", color = "black", outlier.shape = 16, outlier.size = 2) +  # Customize the boxplot
  theme_minimal() +
  labs(
    title = "Larvae Number confirmation depending on marker chr5 in north cluster (A. fabae)",
    x = "marker on chr5",
    y = "Larvae Total Number"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels for better readability


data_marker_east %>% 
  mutate(chr5 = factor(chr5, levels = c(0, 1, 2))) %>%  
  ggplot(aes(x = chr5, y = emmean)) +  # Boxplot aesthetics
  geom_boxplot(fill = "skyblue", color = "black", outlier.shape = 16, outlier.size = 2) +  # Customize the boxplot
  theme_minimal() +
  labs(
    title = "Larvae Number confirmation depending on marker chr5 in east cluster (A. fabae)",
    x = "marker on chr5",
    y = "Larvae Total Number"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels for better readability



data_marker_geneticcluster %>% 
  mutate(chr5 = factor(chr5, levels = c(0, 1, 2))) %>%  
  ggplot(aes(x = chr5, y = emmean)) +  # Boxplot aesthetics
  geom_boxplot(fill = "skyblue", color = "black", outlier.shape = 16, outlier.size = 2) +  # Customize the boxplot
  theme_minimal() +
  labs(
    title = "Larvae Number confirmation depending on marker chr5 (A. fabae)",
    x = "marker on chr5",
    y = "Larvae Total Number"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels for better readability

mod<-lm(emmean ~ chr5, data=data_marker_geneticcluster)
anova(mod)
mod<-lm(emmean ~ chr5, data=data_marker_south)
anova(mod)
#verification----

data_count <- median_xp007_xp009_marker_geneticcluster %>%
  group_by(chr5) %>%
  summarise(n = n())




summary_df <- median_xp007_xp009_marker_geneticcluster %>%
  group_by(cluster_name) %>%
  summarise(num_twos = sum(chr5 == 2))












#########################################################################
#########################################################################
#distribution chr1L among clusters (big screening)----





# Count chr1L values by cluster
summary_df <- data_marker_geneticcluster %>%
  group_by(cluster_name, chr1L) %>%
  summarise(count = n(), .groups = "drop") %>%
  filter(chr1L %in% c(0, 1, 2)) %>%
  mutate(chr1L = as.factor(chr1L))

# Ensure all combinations are present
summary_df_complete <- summary_df %>%
  complete(cluster_name, chr1L = factor(0:2), fill = list(count = 0))

# Add "Total" column
total_counts <- summary_df_complete %>%
  group_by(chr1L) %>%
  summarise(count = sum(count)) %>%
  mutate(cluster_name = "Total")

# Combine and set cluster_name order (Total last)
final_df <- bind_rows(summary_df_complete, total_counts)
final_df$cluster_name <- factor(final_df$cluster_name, 
                                levels = c(setdiff(unique(final_df$cluster_name), "Total"), "Total"))

# Plot
ggplot(final_df, aes(x = cluster_name, y = count, fill = chr1L)) +
  geom_bar(stat = "identity") +
  labs(
    title = "Counts of 0s, 1s, and 2s in chr1L by Cluster",
    x = "Cluster Name",
    y = "Count",
    fill = "marker on chr1L"
  ) +
  theme_minimal()

data_count <- data_marker_south %>%
  group_by(chr1L) %>%
  summarise(n = n())

data_count <- data_marker_geneticcluster %>%
  group_by(chr1L) %>%
  summarise(n = n())


summary_df <- data_marker_geneticcluster %>%
  group_by(cluster_name) %>%
  summarise(num_twos = sum(chr1L == 2))




#analyse chr1L (big screening) ----

# Calculate the number of observations for each level of chr1L


data_marker_south %>% 
  mutate(chr1L = factor(chr1L, levels = c(0, 1, 2))) %>%  
  ggplot(aes(x = chr1L, y = emmean)) +  # Boxplot aesthetics
  geom_boxplot(fill = "skyblue", color = "black", outlier.shape = 16, outlier.size = 2) +  # Customize the boxplot
  theme_minimal() +
  labs(
    title = "Larvae Number confirmation depending on marker chr1L in south cluster (A. fabae)",
    x = "marker on chr1L",
    y = "Larvae Total Number"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels for better readability




data_marker_north %>% 
  mutate(chr1L = factor(chr1L, levels = c(0, 1, 2))) %>%  
  ggplot(aes(x = chr1L, y = emmean)) +  # Boxplot aesthetics
  geom_boxplot(fill = "skyblue", color = "black", outlier.shape = 16, outlier.size = 2) +  # Customize the boxplot
  theme_minimal() +
  labs(
    title = "Larvae Number confirmation depending on marker chr1L in north cluster (A. fabae)",
    x = "marker on chr1L",
    y = "Larvae Total Number"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels for better readability



data_marker_geneticcluster %>% 
  mutate(chr1L = factor(chr1L, levels = c(0, 1, 2))) %>%  
  ggplot(aes(x = chr1L, y = emmean)) +  # Boxplot aesthetics
  geom_boxplot(fill = "skyblue", color = "black", outlier.shape = 16, outlier.size = 2) +  # Customize the boxplot
  theme_minimal() +
  labs(
    title = "Larvae Number confirmation depending on marker chr1L (A. fabae)",
    x = "marker on chr1L",
    y = "Larvae Total Number"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels for better readability

mod<-lm(emmean ~ chr1L, data=data_marker_geneticcluster)
anova(mod)
mod<-lm(emmean ~ chr1L, data=data_marker_south)
anova(mod)
#verification----

data_count <- median_xp007_xp009_marker_geneticcluster %>%
  group_by(chr1L) %>%
  summarise(n = n())




summary_df <- median_xp007_xp009_marker_geneticcluster %>%
  group_by(cluster_name) %>%
  summarise(num_twos = sum(chr1L == 2))


########################################################################
#distribution chr3 among clusters (big screening)----





# Count chr3 values by cluster
summary_df <- data_marker_geneticcluster %>%
  group_by(cluster_name, chr3) %>%
  summarise(count = n(), .groups = "drop") %>%
  filter(chr3 %in% c(0, 1, 2)) %>%
  mutate(chr3 = as.factor(chr3))

# Ensure all combinations are present
summary_df_complete <- summary_df %>%
  complete(cluster_name, chr3 = factor(0:2), fill = list(count = 0))

# Add "Total" column
total_counts <- summary_df_complete %>%
  group_by(chr3) %>%
  summarise(count = sum(count)) %>%
  mutate(cluster_name = "Total")

# Combine and set cluster_name order (Total last)
final_df <- bind_rows(summary_df_complete, total_counts)
final_df$cluster_name <- factor(final_df$cluster_name, 
                                levels = c(setdiff(unique(final_df$cluster_name), "Total"), "Total"))

# Plot
ggplot(final_df, aes(x = cluster_name, y = count, fill = chr3)) +
  geom_bar(stat = "identity") +
  labs(
    title = "Counts of 0s, 1s, and 2s in chr3 by Cluster",
    x = "Cluster Name",
    y = "Count",
    fill = "marker on chr3"
  ) +
  theme_minimal()

data_count <- data_marker_south %>%
  group_by(chr3) %>%
  summarise(n = n())

data_count <- data_marker_geneticcluster %>%
  group_by(chr3) %>%
  summarise(n = n())


summary_df <- data_marker_geneticcluster %>%
  group_by(cluster_name) %>%
  summarise(num_twos = sum(chr3 == 2))




#analyse chr3 (big screening) ----

# Calculate the number of observations for each level of chr3


data_marker_south %>% 
  mutate(chr3 = factor(chr3, levels = c(0, 1, 2))) %>%  
  ggplot(aes(x = chr3, y = emmean)) +  # Boxplot aesthetics
  geom_boxplot(fill = "skyblue", color = "black", outlier.shape = 16, outlier.size = 2) +  # Customize the boxplot
  theme_minimal() +
  labs(
    title = "Larvae Number confirmation depending on marker chr3 in south cluster (A. fabae)",
    x = "marker on chr3",
    y = "Larvae Total Number"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels for better readability




data_marker_north %>% 
  mutate(chr3 = factor(chr3, levels = c(0, 1, 2))) %>%  
  ggplot(aes(x = chr3, y = emmean)) +  # Boxplot aesthetics
  geom_boxplot(fill = "skyblue", color = "black", outlier.shape = 16, outlier.size = 2) +  # Customize the boxplot
  theme_minimal() +
  labs(
    title = "Larvae Number confirmation depending on marker chr3 in north cluster (A. fabae)",
    x = "marker on chr3",
    y = "Larvae Total Number"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels for better readability



data_marker_geneticcluster %>% 
  mutate(chr3 = factor(chr3, levels = c(0, 1, 2))) %>%  
  ggplot(aes(x = chr3, y = emmean)) +  # Boxplot aesthetics
  geom_boxplot(fill = "skyblue", color = "black", outlier.shape = 16, outlier.size = 2) +  # Customize the boxplot
  theme_minimal() +
  labs(
    title = "Larvae Number confirmation depending on marker chr3 in south cluster (A. fabae)",
    x = "marker on chr3",
    y = "Larvae Total Number"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels for better readability

mod<-lm(emmean ~ chr3 + chr5 + chr1L, data=data_marker_geneticcluster)
anova(mod)
mod<-lm(emmean ~  chr3 + chr5 + chr1L, data=data_marker_south)
anova(mod)

mod<-lm(Median_fecondity ~ chr3 + chr5 + chr1L, data=xp007_xp009_median_marker)
anova(mod)
mod<-lm(Median_fecondity ~  chr3 + chr5 + chr1L, data=median_xp007_xp009_marker_south)
anova(mod)

#verification----

data_count <- median_xp007_xp009_marker_geneticcluster %>%
  group_by(chr3) %>%
  summarise(n = n())




summary_df <- median_xp007_xp009_marker_geneticcluster %>%
  group_by(cluster_name) %>%
  summarise(num_twos = sum(chr3 == 2))







####
#Part 2
data_filtered<-read.csv("../faba_specifics/02.analyses/01.GWAS/01_phenotipic_analysis/01.repetability_choice_Nrep/01_A.fabae/Final filters/aphis_fabae_data_filtered.csv")

data_filtered_marker <- merge(data_filtered, extract_marker, by = "ivip")

data_marker_geneticcluster <- geneticcluster %>% dplyr::select(ivip, cluster_name) %>% 
  distinct() %>%
  inner_join(., data_filtered_marker, by = "ivip")

data_marker_geneticcluster_south <- subset(data_marker_geneticcluster, cluster_name %in% c("South", "South_America"))



mod<-lm(larvae_total_number ~ chr3 + chr5 + chr1L, data=data_marker_geneticcluster)
anova(mod)
mod<-lm(larvae_total_number ~ ivip + chr3 + chr5 + chr1L, data=data_marker_geneticcluster_south)
anova(mod)

mod<-lm(Median_fecondity ~ chr3 + chr5 + chr1L, data=xp007_xp009_median_marker)
anova(mod)
mod<-lm(Median_fecondity ~  chr3 + chr5 + chr1L, data=median_xp007_xp009_marker_south)
anova(mod)



install.packages("lmerTest")
library(lmerTest)
##genotype as random factor
#all cluster big screening
mod_mixed <- lmer(larvae_total_number ~ chr3 + chr5 + chr1L + (1|ivip), data = data_marker_geneticcluster)
anova(mod_mixed)

#only south big screening
mod_mixed <- lmer(larvae_total_number ~ chr3 + chr5 + chr1L + (1|ivip), data = data_marker_geneticcluster_south)
anova(mod_mixed)

#all cluster confirmation
mod_mixed <- lmer(larvae_total_number ~ chr3 + chr5 + chr1L + (1|ivip), data = data_filtered_xp007_xp009_marker)
anova(mod_mixed)

#only south confirmation
mod_mixed <- lmer(larvae_total_number ~ chr3 + chr5 + chr1L + (1|ivip), data = data_filtered_xp007_xp009_marker_south)
anova(mod_mixed)


##no genotype
#all cluster big screening
mod_mixed <- lm(larvae_total_number ~ chr3 + chr5 + chr1L , data = data_marker_geneticcluster)
anova(mod_mixed)

#only south big screening
mod_mixed <- lm(larvae_total_number ~ chr3 + chr5 + chr1L , data = data_marker_geneticcluster_south)
anova(mod_mixed)

#all cluster confirmation
mod_mixed <- lm(larvae_total_number ~ chr3 + chr5 + chr1L, data = data_filtered_xp007_xp009_marker)
anova(mod_mixed)

#only south confirmation
mod_mixed <- lm(larvae_total_number ~ chr3 + chr5 + chr1L , data = data_filtered_xp007_xp009_marker_south)
anova(mod_mixed)
