## import packages----
library(data.table)
library(vcfR)   
library(tidyverse)
library(dplyr)
library(gplots)
library(openxlsx)


##Load datasets----
myGD2 <- fread(file = "../../../../../files_too_heavy_for_git/pop_analysis/vfaba_filtered_imputed_maf.hedin.012.indv",header = FALSE, stringsAsFactors = FALSE)
names(myGD2)[names(myGD2) == "V1"] <- "ivip"
myGD2$ivip <- gsub("hedin2", "Hedin", myGD2$ivip)
myKI <- fread(file = "../../../../../files_too_heavy_for_git/pop_analysis/kinship/kinship_pop_analysis.cXX.txt")
myGD<-fread(file = "../../../../../files_too_heavy_for_git/pop_analysis/vfaba_filtered_imputed_maf.hedin.012")
myGD <- myGD[, -1]

##Prepare the datasets----
# Extract genotype names
geno_names <- as.character(myGD2[[1]])

# Convert to numeric matrix
kinship <- as.matrix(myKI)
mode(kinship) <- "numeric"

storage.mode(kinship) <- "numeric"

# Assign row and column names
rownames(kinship) <- geno_names

colnames(kinship) <- geno_names


##generate the kinship and clustering tree----
# Open PDF device
pdf("Kinship_Heatmap_Ward_D2.pdf", width = 12, height = 12)


hm <- heatmap.2(kinship,
                col = rev(heat.colors(256)),
                scale = "none",
                trace = "none",
                density.info = "none",
                key = TRUE,
                symkey = FALSE,
                cexRow = 0.3,
                cexCol = 0.3,
                margins = c(3,3),
                hclustfun = function(x) hclust(x, method = "ward.D2"))

# Close PDF device
dev.off()

postscript(
  "Kinship_Heatmap_Ward_D2.eps",
  onefile = FALSE,
  paper = "special",
  width = 12,
  height = 12,
  horizontal = FALSE
)

hm <- heatmap.2(
  kinship,
  col = rev(heat.colors(256)),
  scale = "none",
  trace = "none",
  density.info = "none",
  key = TRUE,
  symkey = FALSE,
  cexRow = 0.3,
  cexCol = 0.3,
  margins = c(3, 3),
  hclustfun = function(x) hclust(x, method = "ward.D2")
)

dev.off()

# Extract row dendrogram
row_dend <- hm$rowDendrogram

# Convert dendrogram to hclust
row_hc <- as.hclust(row_dend)

## chose best value for k ----


# Elbow Method for optimal K


# Define range of K values to test
k_values <- 1:9

# Initialize vector to store total within-cluster sum of squares
wss <- numeric(length(k_values))

# Convert kinship to distance matrix (Euclidean)
dist_mat <- dist(kinship, method = "euclidean")

for (i in seq_along(k_values)) {
  k <- k_values[i]
  
  # Cut tree into k clusters
  cluster_assign <- cutree(row_hc, k = k)
  
  # Compute total within-cluster sum of squares
  wss_k <- 0
  
  for (cl in unique(cluster_assign)) {
    cluster_data <- kinship[cluster_assign == cl, , drop = FALSE]
    
    # Cluster centroid
    centroid <- colMeans(cluster_data)
    
    # Sum of squared distances to centroid
    wss_k <- wss_k + sum(rowSums((cluster_data - matrix(centroid,
                                                        nrow = nrow(cluster_data),
                                                        ncol = ncol(cluster_data),
                                                        byrow = TRUE))^2))
  }
  
  wss[i] <- wss_k
}

# Plot elbow curve

plot(k_values, wss,
     type = "b",
     pch = 19,
     frame = FALSE,
     xlab = "Number of clusters (K)",
     ylab = "Total Within-Cluster Sum of Squares",
     main = "Elbow Method for Optimal K")






# silouhette Method for optimal K


library(cluster)

# Convert kinship to distance matrix
dist_mat <- dist(kinship, method = "euclidean")

# Test range of K
k_values <- 2:10
avg_sil <- numeric(length(k_values))

for (i in seq_along(k_values)) {
  k <- k_values[i]
  cluster_assign <- cutree(row_hc, k = k)
  sil <- silhouette(cluster_assign, dist_mat)
  avg_sil[i] <- mean(sil[, 3])  # Average silhouette width
}

# Plot
plot(k_values, avg_sil, type="b", pch=19,
     xlab="Number of clusters (K)",
     ylab="Average silhouette width",
     main="Silhouette Method")


## prepare data to create admixture figure----

# Cut into 3 clusters
groups <- cutree(row_hc, k = 3)
table(groups)

groups<-as.data.frame(groups)
groups$ivip<-rownames(groups)
rownames(groups) <- NULL



groups$groups <- recode(groups$groups,
                 `1` = 2,
                 `2` = 1,
                 `3` = 3)


# Row order (top to bottom)
row_order <- hm$rowInd
rownames_in_plot <- rownames(kinship)[row_order]
rownames_in_plot <-as.data.frame(rownames_in_plot)
colnames(rownames_in_plot)[1] <- "ivip"



# Read admixture Q matrix
tbl <- read.table("../../../../../files_too_heavy_for_git/pop_analysis/admixture/vfaba_filtered_imputed_maf.hedin_pruned_for_admixture.3.Q")

# Assign genotype names from myGD2
tbl$ivip <- myGD2$ivip   # replace $GenotypeName with the correct column name

tbl_ordered <- merge(rownames_in_plot, tbl, by = "ivip", sort = FALSE)
tbl_ordered<-tbl_ordered[,-1]
## Plot admixture ----
pdf("admixture.pdf", width = 12, height = 4)

barplot(t(as.matrix(tbl_ordered)), 
        col=rainbow(3),
        xlab="Individual #", 
        ylab="Ancestry", 
        border=NA)

dev.off()

## Genarate PCOA  ----


# 1. Prepare SNP matrix
geno <- as.matrix(myGD)     # remove ivip column
rownames(geno) <- myGD2$ivip      # set genotype names as rownames

# 2. Compute genetic distance

# Using simple matching (proportion of different SNPs)
dist_mat <- dist(geno, method = "manhattan") / ncol(geno)  # scale to 0-1


# 3. PCoA
# Compute PCoA with top 3 axes
pcoa <- cmdscale(dist_mat, k = 3, eig = TRUE)

# Eigenvalues
eig_vals <- pcoa$eig

# Keep only positive eigenvalues (important for PCoA)
eig_vals_pos <- eig_vals[eig_vals > 0]

# Percent variance explained
var_explained <- eig_vals_pos / sum(eig_vals_pos) * 100

# Show first 3 axes
var_explained[1:3]

# Define colors for each group
palette <- c(
  "3" = "red",
  "2" = "blue",
  "1" = "green"
)

# Map colors to each genotype
color_map <- setNames(palette[as.character(groups$groups)], groups$ivip)
cols <- color_map[rownames(geno)]

pdf("PCoA.pdf", width = 12, height = 12)

# Plot PCoA colored by group
plot(pcoa$points[,1], pcoa$points[,2],
     xlab = "PCoA1",
     ylab = "PCoA2",
     main = "PCoA of Genotypes by Group",
     pch = 19,
     col = cols)

# Add legend
legend("topright",
       legend = names(palette),
       col = palette,
       pch = 19,
       cex = 1)

dev.off()

postscript(
  "PCoA.eps",
  onefile = FALSE,
  paper = "special",
  width = 12,
  height = 12,
  horizontal = FALSE
)

# Plot PCoA colored by group
plot(pcoa$points[,1], pcoa$points[,2],
     xlab = "PCoA1",
     ylab = "PCoA2",
     main = "PCoA of Genotypes by Group",
     pch = 19,
     col = cols)

# Add legend
legend("topright",
       legend = names(palette),
       col = palette,
       pch = 19,
       cex = 1)

dev.off()

## Save the group names in the updated xls file  -----
passport_data<-read.csv("../../../00.info/00.passeport_data/2024-11_passeport_data_and_genetic_cluster_viciafaba_peamust.csv")

comparison_cluster<-merge(passport_data, groups, by = "ivip", all.x = TRUE)
write.xlsx(comparison_cluster, "passport_data_updated.xlsx")
colnames(groups)<-NULL
groups[, c(2, 1)] <- groups[, c(1, 2)]
write.table(groups, "../../../00.info/11.infos_genotypes/metadata.csv", col.names = FALSE, row.names = FALSE, sep = ",")
groups[, c(2, 1)] <- groups[, c(1, 2)]

list_origin <- merge(
  rownames_in_plot,
  comparison_cluster,
  by = "ivip",
  all.x = TRUE,
  sort = FALSE
)

list_origin <- list_origin[, c(1, 6,8,9,10)]
list_origin <- list_origin[nrow(list_origin):1, ]
## generate phenotype figure-----
# Load data
emmeans_fabae <- read.csv("../../../02.analyses/01.GWAS/01_phenotipic_analysis/02.ajusted_means/fabae/emmeans_fabae_4.csv", sep = ",")
emmeans_fabae <- emmeans_fabae[, c(1,2)]
colnames(emmeans_fabae)[2] <- "emmeans_fabae"

emmeans_pisum <- read.csv("../../../02.analyses/01.GWAS/01_phenotipic_analysis/02.ajusted_means/pisum/emmeans_pisum_6.csv", sep = ",")
emmeans_pisum <- emmeans_pisum[, c(1,2)]
colnames(emmeans_pisum)[2] <- "emmeans_pisum"

# Merge data (preserve order)
merged <- merge(rownames_in_plot, emmeans_pisum, by = "ivip", all.x = TRUE, sort = FALSE)
merged <- merge(merged, emmeans_fabae, by = "ivip", all.x = TRUE, sort = FALSE)

# Reorder rows
merged <- merged[match(rownames_in_plot$ivip, merged$ivip), ]
rownames(merged) <- NULL

# Reshape to long format
merged_long <- merged |>
  pivot_longer(cols = c(emmeans_fabae, emmeans_pisum),
               names_to = "trait",
               values_to = "value")



p <- ggplot(merged_long, aes(
  x = 1, 
  y = factor(ivip, levels = rownames_in_plot$ivip), 
  fill = value
)) +
  geom_tile() +
  scale_fill_gradient(
    low = "green", 
    high = "red", 
    na.value = "lightgrey", 
    labels = number_format(accuracy = 1),
    breaks = seq(data_min, data_max, by = 10),
    guide = guide_colorbar(
      title.position = "top",   # move title above
      title.hjust = 0.5,       # center title
      barwidth = 12,           # make bar longer
      barheight = 4            # make bar shorter
    )
  ) +
  facet_wrap(~trait, ncol = 2) +
  labs(x = "", y = "Row (ivip)", fill = "Color Key") +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    legend.title = element_text(size = 18, face = "bold"),
    legend.text = element_text(size = 14),
    legend.position = "bottom",
    legend.direction = "horizontal"
  )

# Save
ggsave("strip_combined.png", plot = p, width = 4, height = 15, units = "in", dpi = 300)

ggsave(
  "strip_combined.eps",
  plot = p,
  device = "eps",
  width = 4,
  height = 15,
  units = "in"
)