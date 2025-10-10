library(dplyr)
library(ggplot2)
library(gridExtra)

original_list <- read.csv(file = "00.info/11.infos_genotypes/morphological_phenotlogical_data_faba_collection_PeaMUST_with_genet_cluster.csv", header = TRUE, stringsAsFactors = FALSE, sep =",") 

fabae_no_filters <- read.csv(file = "01.data/2023_screen_specifics_251genotypes/02.aphis_fabae_on_faba_collection/aphis_fabae_fecundity_on_faba_collection_2024-04-02.csv", header = TRUE, stringsAsFactors = FALSE, sep =",")

fabae_filtered <- read.csv(file = "02.analyses/01.GWAS/01_phenotipic_analysis/01.dataset_filtering/01_A.fabae/final filters/aphis_fabae_data_filtered_6.csv", header = TRUE, stringsAsFactors = FALSE, sep =",")


# Step 1: Get unique genotype names from each dataset
full_list_genotypes <- unique(original_list$CODE.IVIP)
genotype_fabae_no_filters <- unique(fabae_no_filters$ivip)
genotypes_fabae_filtered <- unique(fabae_filtered$ivip)

length(genotypes_fabae_filtered)
length(full_list_genotypes)



# Step 2: Find genotypes not in common
not_in_the_dataframe_fabae <- setdiff(full_list_genotypes, genotype_fabae_no_filters)
lost_in_filtering_fabae <- setdiff(genotype_fabae_no_filters, genotypes_fabae_filtered)
all_missing_fabae <- setdiff(full_list_genotypes, genotypes_fabae_filtered)




df1_subset <- original_list[original_list$CODE.IVIP %in% all_missing_fabae, ]

original_list$dataset <- "full list of genotypes"
df1_subset$dataset <- "missing genotypes"

combined<-rbind(original_list,df1_subset)

ggplot(combined, aes(x = SumTem, fill = dataset)) +
  geom_histogram(alpha = 0.5, position = "identity", bins = 30) +
  labs(title = "Value Distribution", x = "Value", y = "Count") +
  theme_minimal()


cluster_no_filters <- as.data.frame(table(original_list$ClusterName))
colnames(cluster_no_filters) <- c("cluster", "count")

par(mfrow = c(1, 2))

cluster_no_filters<-cluster_no_filters[-c(1),]

# Make pie chart
p1<-ggplot(cluster_no_filters, aes(x = "", y = count, fill = cluster)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar("y") +
  labs(title = "cluster Distribution before filters") +
  theme_void()

cluster_filtered <- as.data.frame(table(df1_subset$ClusterName))
colnames(cluster_filtered) <- c("cluster", "count")
cluster_filtered<-cluster_filtered[-c(1),]
# Make pie chart
p2<-ggplot(cluster_filtered, aes(x = "", y = count, fill = cluster)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar("y") +
  labs(title = "cluster Distribution after filters") +
  theme_void()

grid.arrange(p1, p2, ncol = 2)



df2_subset <- original_list[original_list$CODE.IVIP %in% genotypes_fabae_filtered, ]


df2_subset$dataset <- "genotypes remaining after filter"

combined2<-rbind(original_list,df2_subset)


combined2 <- combined2 %>%
  filter(!is.na(ClusterName) & ClusterName != "")

ggplot(combined2, aes(x = ClusterName, fill = dataset)) +
  geom_bar(position = "dodge") +
  labs(title = "Count of var1 grouped by var2",
       x = "var1",
       y = "Count") +
  theme_minimal()



################################################################################
################################################################################
################################################################################



pisum_no_filters <- read.csv(file = "01.data/2023_screen_specifics_251genotypes/01.acyrthosiphon_pisum_P123amp_on_faba_collection/acyrthospipon_pisum_fecundity_on_faba_collection.csv", header = TRUE, stringsAsFactors = FALSE, sep =",")

pisum_filtered <- read.csv(file = "02.analyses/01.GWAS/01_phenotipic_analysis/01.dataset_filtering/02_A.pisum/final filters/acyrthosiphon_pisum_filtered_6.csv", header = TRUE, stringsAsFactors = FALSE, sep =",")

# Step 1: Get unique genotype names from each dataset
full_list_genotypes <- unique(original_list$CODE.IVIP)
genotype_pisum_no_filters <- unique(pisum_no_filters$ivip)
genotypes_pisum_filtered <- unique(pisum_filtered$ivip)

length(genotypes_pisum_filtered)
length(full_list_genotypes)



# Step 2: Find genotypes not in common
# Genotypes in df1 but not in df2
not_in_the_dataframe_pisum <- setdiff(full_list_genotypes, genotype_pisum_no_filters)
lost_in_filtering_pisum <- setdiff(genotype_pisum_no_filters, genotypes_pisum_filtered)
all_missing_pisum <- setdiff(full_list_genotypes, genotypes_pisum_filtered)




df1_subset <- original_list[original_list$CODE.IVIP %in% all_missing_pisum, ]

original_list$dataset <- "full list of genotypes"
df1_subset$dataset <- "missing genotypes"

combined<-rbind(original_list,df1_subset)

ggplot(combined, aes(x = SumTem, fill = dataset)) +
  geom_histogram(alpha = 0.5, position = "identity", bins = 30) +
  labs(title = "Value Distribution", x = "Value", y = "Count") +
  theme_minimal()


cluster_no_filters <- as.data.frame(table(original_list$ClusterName))
colnames(cluster_no_filters) <- c("cluster", "count")

par(mfrow = c(1, 2))

cluster_no_filters<-cluster_no_filters[-c(1),]

# Make pie chart
p1<-ggplot(cluster_no_filters, aes(x = "", y = count, fill = cluster)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar("y") +
  labs(title = "cluster Distribution before filters") +
  theme_void()

cluster_filtered <- as.data.frame(table(df1_subset$ClusterName))
colnames(cluster_filtered) <- c("cluster", "count")
cluster_filtered<-cluster_filtered[-c(1),]
# Make pie chart
p2<-ggplot(cluster_filtered, aes(x = "", y = count, fill = cluster)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar("y") +
  labs(title = "cluster Distribution after filters") +
  theme_void()

grid.arrange(p1, p2, ncol = 2)



df2_subset <- original_list[original_list$CODE.IVIP %in% genotypes_pisum_filtered, ]


df2_subset$dataset <- "genotypes remaining after filter"

combined2<-rbind(original_list,df2_subset)


combined2 <- combined2 %>%
  filter(!is.na(ClusterName) & ClusterName != "")

ggplot(combined2, aes(x = ClusterName, fill = dataset)) +
  geom_bar(position = "dodge") +
  labs(title = "Count of var1 grouped by var2",
       x = "var1",
       y = "Count") +
  theme_minimal()


################################################################################
################################################################################
################################################################################

germination <- read.csv(file = "../faba_specifics/00.info/11.infos_genotypes/SPECIFICS_Germination.csv",header=TRUE, stringsAsFactors = TRUE, sep = ";")

genotype_list <- read.csv(file = "../faba_specifics/00.info/11.infos_genotypes/SPECIFICS_GénotypesPisumFabae2.csv",header=TRUE, stringsAsFactors = TRUE, sep = ",")

# Example:
merged_data <- merge(genotype_list, germination, by = "genotype")

merged_data <- merged_data %>%
  mutate(missing_fabae = ivip %in% all_missing_fabae)

merged_data <- merged_data %>%
  mutate(missing_pisum = ivip %in% all_missing_pisum)
