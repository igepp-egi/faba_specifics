# in the terminal
# awk 'NR==FNR {a[$1]; next} FNR==1 || $1 in a' chr2_1461736473_small_window.txt vfaba.genotypes.tsv > vfaba.filtered_big.tsv
library(data.table)
library(pheatmap)
library(ggplot2)
library(emmeans)
library(multcomp)

# -----------------------------
# 1. Load filtered TSV
# -----------------------------
df <- fread("vfaba.filtered.tsv")
df <- fread("vfaba.filtered_big.tsv")

# Add reference genome column "hedin" with default value 0|0 for all SNPs
df$Hedin <- "0|0"

# Optional: reorder columns to place hedin first
df <- df[, c("Hedin", setdiff(names(df), "Hedin")), with = FALSE]


# Store snp_id separately
snp_ids <- df$snp_id
df$snp_id <- NULL

# -----------------------------
# 2. Convert genotype strings to numeric
# -----------------------------
geno_to_num <- function(g) {
  if(g %in% c("0|0", "0/0")) return(0)
  else if(g %in% c("1|1", "1/1")) return(2)
  else return(1)  # heterozygous
}
df_num <- as.data.frame(lapply(df, function(col) sapply(col, geno_to_num)))

# -----------------------------
# 3. Transpose: genotypes as rows, SNPs as columns
# -----------------------------
df_num_t <- t(df_num)
rownames(df_num_t) <- colnames(df)

# -----------------------------
# 4. Hierarchical clustering on genotypes
# -----------------------------
d <- dist(df_num_t, method = "euclidean")  # distance between genotypes
hc <- hclust(d, method = "average")        # hierarchical clustering

# -----------------------------
# 5. Cut tree into clusters
# -----------------------------
k <- 5  # number of clusters, adjust as needed
groups <- cutree(hc, k = k)

# -----------------------------
# 6. Create cluster table
# -----------------------------
cluster_table <- data.frame(
  Genotype = names(groups),
  Group = groups
)
fwrite(cluster_table, "vfaba.genotype_clusters.tsv", sep = "\t")

# -----------------------------
# 7. Reorder dataframe by cluster order for TSV
# -----------------------------
ordered_genotypes <- hc$labels[hc$order]
df_ordered <- df_num_t[ordered_genotypes, ]
df_ordered <- cbind(Genotype = rownames(df_ordered), df_ordered)
fwrite(df_ordered, "vfaba.filtered.ordered.tsv", sep = "\t")

# -----------------------------
# 8. Generate heatmap PNG with clusters
# -----------------------------
annotation <- data.frame(Group = factor(groups))
rownames(annotation) <- names(groups)

pheatmap(df_num_t,
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         clustering_method = "average",
         annotation_row = annotation,
         main = "Genotype similarity with groups",
         filename = "vfaba.filtered.heatmap_groups.png",
         width = 10,
         height = 12)


# Create column from row names
annotation$ivip <- rownames(annotation)

# Remove the row names
rownames(annotation) <- NULL

# move ivip to the first column
annotation <- annotation[, c("ivip", setdiff(names(annotation), "ivip"))]




# test association with phenotype:

pheno <-read.csv("../../../../gitrepo/faba_specifics/02.analyses/01.GWAS/01_phenotipic_analysis/01.dataset_filtering/02_A.pisum/final filters/acyrthosiphon_pisum_filtered_6.csv")
summary(pheno)
merged_df <- merge(pheno, annotation, by = "ivip", all.x = TRUE)


# Assuming 'trait1' is your phenotype, 'Group' is haplotype, 'ivip' is genotype
model <- lmer(larvae_total_number ~ Group + (1|ivip), data = merged_df)

# Check ANOVA table
anova(model)
summary(model)

# Estimated marginal means
emm <- emmeans(model, ~ Group)

# Tukey-adjusted pairwise comparisons
cld <- multcomp::cld(emm, Letters = letters, adjust = "tukey")

# View letters
cld


# Basic boxplot using base R
boxplot(larvae_total_number ~ Group, 
        data = merged_df, 
        main = "Boxplot of Larvae Total Number by Group",
        xlab = "Group",
        ylab = "Larvae Total Number",
        col = "lightblue")

