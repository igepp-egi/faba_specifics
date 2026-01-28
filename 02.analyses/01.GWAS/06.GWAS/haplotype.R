# in the terminal
# awk 'NR==FNR {a[$1]; next} FNR==1 || $1 in a' chr1S_1240520268_small_window.txt vfaba.genotypes.tsv > vfaba.filtered_big.tsv

# -----------------------------
# 0. Load required packages
# -----------------------------
library(data.table)
library(pheatmap)
library(ggplot2)
library(emmeans)
library(multcomp)
library(RColorBrewer)
library(lme4)


# -----------------------------
# 1. Load filtered TSV
# -----------------------------
df <- fread("vfaba.filtered_025.tsv")


df$hedin <- "0|0"  # add reference genome column
df <- df[, c("hedin", setdiff(names(df), "hedin")), with = FALSE]

# Store SNP IDs separately
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
rownames(df_num_t) <- colnames(df)               # rownames = genotype IDs
colnames(df_num_t) <- snp_ids

# -----------------------------
# 4. Hierarchical clustering on genotypes
# -----------------------------
d <- dist(df_num_t, method = "euclidean")
hc <- hclust(d, method = "average")

# -----------------------------
# 5. Cut tree into clusters
# -----------------------------
k <- 4   # number of clusters
groups <- cutree(hc, k = k)

# -----------------------------
# 6. Create cluster table
# -----------------------------
cluster_table <- data.frame(
  ivip = names(groups),
  Group = groups, row.names = NULL
)
fwrite(cluster_table, "vfaba.genotype_clusters.tsv", sep = "\t")

# -----------------------------
# 7. Reorder dataframe by cluster order for TSV
# -----------------------------
ordered_genotypes <- hc$labels[hc$order]
df_ordered <- df_num_t[ordered_genotypes, ]
df_ordered <- cbind(Genotype = rownames(df_ordered), df_ordered)
write(df_ordered, "vfaba.filtered.ordered.tsv", sep = "\t")

# -----------------------------
# 8. Load metadata and prepare annotations
# -----------------------------
genet_cluster <- read.table("../../../crosshap/metadata.csv", header = FALSE, sep = ",")
colnames(genet_cluster) <- c("ivip", "cluster")

emmeans <- read.csv("../../../../gitrepo/faba_specifics/02.analyses/01.GWAS/01_phenotipic_analysis/02.ajusted_means/pisum/emmeans_pisum_6.csv")

annotation <- merge(emmeans,genet_cluster, by = "ivip", all.x = TRUE)
annotation <- merge(annotation, cluster_table, by = "ivip", all.x = TRUE)

# -----------------------------
# 9. Prepare annotation for pheatmap
# -----------------------------
# Ensure rownames match df_num_t
rownames(annotation) <- annotation$ivip
annotation$ivip <- NULL

# Ensure types
annotation$cluster <- as.factor(annotation$cluster)
annotation$Group <- as.factor(annotation$Group)
annotation$emmean <- as.numeric(annotation$emmean)

# Colors
cluster_colors <- setNames(brewer.pal(length(levels(annotation$cluster)), "Set2"),
                           levels(annotation$cluster))
emmean_colors <- colorRampPalette(c("white", "black"))(100)

annotation_colors <- list(
  cluster = cluster_colors,
  emmean = emmean_colors
)

# -----------------------------
# 10. Generate heatmap PNG
# -----------------------------
pheatmap(df_num_t,
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         clustering_method = "average",
         cluster_cols = FALSE,
         annotation_row = annotation,
         annotation_colors = annotation_colors,
         main = "Genotype similarity with clusters & emmean",
         show_rownames = FALSE,
         filename = "vfaba.filtered.heatmap_clusters_emmean.png",
         width = 10,
         height = 12)





pheno <-read.csv("../../../../gitrepo/faba_specifics/02.analyses/01.GWAS/01_phenotipic_analysis/01.dataset_filtering/02_A.pisum/final filters/acyrthosiphon_pisum_filtered_6.csv")
summary(pheno)
merged_df <- merge(pheno, cluster_table, by = "ivip", all.x = TRUE)
merged_df$Group <- as.factor(merged_df$Group)

model <- lm(larvae_total_number ~ Group + ivip + session/block, data = merged_df)
anova(model)

# Assuming 'trait1' is your phenotype, 'Group' is haplotype, 'ivip' is genotype
model <- lmer(larvae_total_number ~ Group + (1|ivip) + (1 | session/block), data = merged_df)


# Check ANOVA table
anova(model)
summary(model)

# Estimated marginal means
emm <- emmeans(model, ~ Group)

# Tukey-adjusted pairwise comparisons
cld <- multcomp::cld(emm, Letters = letters, adjust = "tukey")

# View letters
cld







png(
  filename = "boxplot_haplotype_comparison_chr1S_1240520268.png",
  width = 2000,
  height = 1800,
  res = 300
)


# Basic boxplot using base R
boxplot(larvae_total_number ~ Group, 
        data = merged_df, 
        main = "Boxplot of Larvae Total Number by Group",
        xlab = "Group",
        ylab = "Larvae Total Number",
        col = "grey")

dev.off()






myGD1 <-fread(file = "../../../../../Documents/results_with_hedin/vfaba_apisum_n6.filtered_imputed_maf.vcf.hedin.vcf.012",header = FALSE, stringsAsFactors =  FALSE) #%>% select(-1)

myGD1 <- myGD1[, -1, with = FALSE]  # Removes the first column

myGM <-fread(file = "../../../../../Documents/results_with_hedin/vfaba_apisum_n6.filtered_imputed_maf.vcf.hedin.vcf.012.pos") %>% rename(
  "chromosome" = "V1",
  "position" = "V2"
)%>% mutate(snp_id=paste(chromosome,position, sep = "_")) %>% relocate(snp_id,.before = chromosome) 

colnames(myGD1)<- myGM$snp_id

myGD2 <-fread(file = "../../../../results_with_hedin/vfaba_apisum_n6.filtered_imputed_maf.vcf.hedin.vcf.012.indv",header = FALSE, stringsAsFactors =  FALSE) #%>% select(-1)


colnames(myGD2)[1]<-"ivip"


myGD1$ivip <-myGD2$ivip

chr2_1461736473 <- myGD1[,c("ivip","chr2_1461736473")]

# In-place replacement
chr2_1461736473[ivip == "Hedin2", ivip := "hedin"]

snp_haplotype_comparison <- merge(merged_df, chr2_1461736473, by = "ivip", all.x = TRUE)

snp_haplotype_comparison$chr2_1461736473 <- as.factor(snp_haplotype_comparison$chr2_1461736473)

model_snp <- lm(larvae_total_number ~ chr2_1461736473 + ivip + session/block, data = snp_haplotype_comparison)

test <- lm(larvae_total_number ~ chr2_1461736473, data = snp_haplotype_comparison)
anova(model_snp)



model_haplotype <- lm(larvae_total_number ~ Group + ivip + session/block, data = snp_haplotype_comparison)
anova(model_haplotype)
summary(model_haplotype)





