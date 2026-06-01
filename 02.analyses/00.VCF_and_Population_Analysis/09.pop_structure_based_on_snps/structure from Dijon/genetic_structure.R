suppressPackageStartupMessages(library("tidyverse"))

load("00.info/03.SNP_data_exome_capture/PCA_new.Rdata")

# Extract and annotate PCA scores (adegenet DAPC)
ivips <- dimnames(PCA$scores)[[1]] # row names of PCA scores matrix
PCs <- dimnames(PCA$scores)[[2]]  # column names of PCA scores matrix  
pca_scores <- PCA$scores 
rownames(pca_scores) <- ivips
colnames(pca_scores) <- PCs  


##############################
# Calculate variance explained
###############################

# Extract eigenvalues
eigenvalues <- PCA$eig^2

# Calculate the total variance
total_variance <- sum(eigenvalues)

# Calculate the proportion of variance explained by each principal component
variance_explained <- eigenvalues / total_variance
variance_explained_df <- data.frame(PC = PCs, variance_explained = variance_explained) %>% 
  arrange(desc(variance_explained)) %>%
  mutate(variance_explained = round(variance_explained * 100, 2))


ggplot(data = variance_explained_df[1:9, ], 
       aes(x = PC, y = variance_explained)) +
  geom_bar(stat = "identity", color='black', fill = "grey") +
  theme_minimal() +
  labs(title = "Variance explained by each principal component", x = "Principal component", y = "Proportion of variance explained") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

###############################
# Add genetic group information
###############################

genetic_clusters <- read.csv("00.info/02.genetic_structure_jonathan_dijon/faba_bean_ivip_and_genetic_cluster.csv", stringsAsFactors = F) %>% 
  select(ivip, cluster_dapc, cluster_name)

pca_scores_with_clusters <- pca_scores %>%
  as.data.frame() %>%
  rownames_to_column(var = "ivip") %>%
  inner_join(genetic_clusters, by = "ivip")

# PCA plot
ggplot(pca_scores, aes(x = PC1, y = PC2, color = cluster_name)) +
  geom_point(size = 3) +
  theme_minimal() +
  theme(legend.position = "bottom") +
  labs(title = "PCA plot", 
       x = paste("PC1 (", round(PCA$var$contrib[1], 2), "%)", sep = ""),
       y = paste("PC2 (", round(PCA$var$contrib[2], 2), "%)", sep = ""),
       color = "Group") +
  scale_color_manual(values = c("red", "blue", "green", "purple", "orange", "black", "brown", "pink", "grey", "yellow"))
  
  