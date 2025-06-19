library(dplyr)
library(ggplot2)
library(gridExtra)
original_list <- read.csv(file = "00.info/11.infos_genotypes/morphological_phenotlogical_data_faba_collection_PeaMUST_with_genet_cluster.csv", header = TRUE, stringsAsFactors = FALSE, sep =",") 

fabae1 <- read.csv(file = "01.data/2023_screen_specifics_251genotypes/02.aphis_fabae_on_faba_collection/aphis_fabae_fecundity_on_faba_collection_2024-04-02.csv", header = TRUE, stringsAsFactors = FALSE, sep =",")

fabae2 <- read.csv(file = "02.analyses/01.GWAS/01_phenotipic_analysis/01.dataset_filtering/01_A.fabae/final filters/aphis_fabae_data_filtered_6.csv", header = TRUE, stringsAsFactors = FALSE, sep =",")


# Step 1: Get unique genotype names from each dataset
genotypes1 <- unique(original_list$CODE.IVIP)
genotypes2 <- unique(fabae1$ivip)
genotypes3 <- unique(fabae2$ivip)

length(genotypes3)
length(genotypes1)



# Step 2: Find genotypes not in common
# Genotypes in df1 but not in df2
only_in_df1 <- setdiff(genotypes1, genotypes2)
only_in_df2 <- setdiff(genotypes2, genotypes3)
only_in_df3 <- setdiff(genotypes1, genotypes3)




df1_subset <- original_list[original_list$CODE.IVIP %in% only_in_df3, ]

original_list$dataset <- "Dataset 1"
df1_subset$dataset <- "Dataset 2"

combined<-rbind(original_list,df1_subset)

ggplot(combined, aes(x = SumTem, fill = dataset)) +
  geom_histogram(alpha = 0.5, position = "identity", bins = 30) +
  labs(title = "Value Distribution", x = "Value", y = "Count") +
  theme_minimal()


test <- as.data.frame(table(original_list$ClusterName))
colnames(test) <- c("cluster", "count")

par(mfrow = c(1, 2))

test<-test[-c(1),]

# Make pie chart
p1<-ggplot(test, aes(x = "", y = count, fill = cluster)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar("y") +
  labs(title = "cluster Distribution") +
  theme_void()

test2 <- as.data.frame(table(df1_subset$ClusterName))
colnames(test2) <- c("cluster", "count")
test2<-test2[-c(1),]
# Make pie chart
p2<-ggplot(test2, aes(x = "", y = count, fill = cluster)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar("y") +
  labs(title = "cluster Distribution") +
  theme_void()

grid.arrange(p1, p2, ncol = 2)



library(ggplot2)



df2_subset <- original_list[original_list$CODE.IVIP %in% genotypes3, ]


df2_subset$dataset <- "Dataset 3"

combined2<-rbind(original_list,df2_subset)


combined2 <- combined2 %>%
  filter(!is.na(ClusterName) & ClusterName != "")

ggplot(combined2, aes(x = ClusterName, fill = dataset)) +
  geom_bar(position = "dodge") +
  labs(title = "Count of var1 grouped by var2",
       x = "var1",
       y = "Count") +
  theme_minimal()

######################################################################
######################################################################
######################################################################

original_list <- read.csv(file = "00.info/11.infos_genotypes/morphological_phenotlogical_data_faba_collection_PeaMUST_with_genet_cluster.csv", header = TRUE, stringsAsFactors = FALSE, sep =",") 

pisum1 <- read.csv(file = "01.data/2023_screen_specifics_251genotypes/01.acyrthosiphon_pisum_P123amp_on_faba_collection/acyrthospipon_pisum_fecundity_on_faba_collection.csv", header = TRUE, stringsAsFactors = FALSE, sep =",")

pisum2 <- read.csv(file = "02.analyses/01.GWAS/01_phenotipic_analysis/01.dataset_filtering/02_A.pisum/final filters/acyrthosiphon_pisum_filtered_6.csv", header = TRUE, stringsAsFactors = FALSE, sep =",")


# Step 1: Get unique genotype names from each dataset
genotypes1 <- unique(original_list$CODE.IVIP)
genotypes2 <- unique(pisum1$ivip)
genotypes3 <- unique(pisum2$ivip)


# Step 2: Find genotypes not in common
# Genotypes in df1 but not in df2
only_in_df1_pisum <- setdiff(genotypes1, genotypes2)
only_in_df2_pisum <- setdiff(genotypes2, genotypes3)
only_in_df3_pisum <- setdiff(genotypes1, genotypes3)




df1_subset <- original_list[original_list$CODE.IVIP %in% only_in_df3_pisum, ]

original_list$dataset <- "Dataset 1"
df1_subset$dataset <- "Dataset 2"

combined<-rbind(original_list,df1_subset)

ggplot(combined, aes(x = SumTem, fill = dataset)) +
  geom_histogram(alpha = 0.5, position = "identity", bins = 30) +
  labs(title = "Value Distribution", x = "Value", y = "Count") +
  theme_minimal()


test <- as.data.frame(table(original_list$ClusterName))
colnames(test) <- c("cluster", "count")

par(mfrow = c(1, 2))

test<-test[-c(1),]

# Make pie chart
p1<-ggplot(test, aes(x = "", y = count, fill = cluster)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar("y") +
  labs(title = "cluster Distribution") +
  theme_void()

test2 <- as.data.frame(table(df1_subset$ClusterName))
colnames(test2) <- c("cluster", "count")
test2<-test2[-c(1),]
# Make pie chart
p2<-ggplot(test2, aes(x = "", y = count, fill = cluster)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar("y") +
  labs(title = "cluster Distribution") +
  theme_void()

grid.arrange(p1, p2, ncol = 2)







df2_subset <- original_list[original_list$CODE.IVIP %in% genotypes3, ]


df2_subset$dataset <- "Dataset 3"

combined2<-rbind(original_list,df2_subset)


combined2 <- combined2 %>%
  filter(!is.na(ClusterName) & ClusterName != "")

ggplot(combined2, aes(x = ClusterName, fill = dataset)) +
  geom_bar(position = "dodge") +
  labs(title = "Count of var1 grouped by var2",
       x = "var1",
       y = "Count") +
  theme_minimal()

fabae <- read.csv(file = "01.data/2023_screen_specifics_251genotypes/02.aphis_fabae_on_faba_collection/aphis_fabae_fecundity_on_faba_collection_2024-04-02.csv", header = TRUE, stringsAsFactors = FALSE, sep =",")




common_names <- intersect(only_in_df3_pisum, only_in_df3)
length(only_in_df3_pisum)
length(only_in_df3)
length(common_names)

test <- intersect(common_names, only_in_df1)


