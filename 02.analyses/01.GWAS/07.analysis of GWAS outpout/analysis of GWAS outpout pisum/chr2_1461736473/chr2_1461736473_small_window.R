library(data.table)
library(dplyr)
library(tidyr)

# -----------------------------
# 1. Load TSV
# -----------------------------
df <- fread("vfaba.genotypes.tsv")

# -----------------------------
# 2. Split SNP ID, filter, and rebuild SNP ID
# -----------------------------
df_filtered <- df %>%
  separate(
    col = snp_id,
    into = c("chromosome", "position"),
    sep = "_",
    convert = TRUE
  ) %>%
  filter(
    chromosome == "chr2",
    position >= 1461736473,
    position <= 1461746277
  ) %>%
  mutate(
    snp_id = paste0(chromosome, "_", position)
  )

# -----------------------------
# 3. Add reference genome column first
# -----------------------------
df_filtered$hedin <- "0|0"

df_filtered <- df_filtered %>%
  relocate(hedin, .before = everything())

# -----------------------------
# 4. Store SNP IDs and drop metadata
# -----------------------------
snp_ids <- df_filtered$snp_id

df_filtered <- df_filtered %>%
  select(-chromosome, -position, -snp_id)

# Clean sample names
setnames(
  df_filtered,
  old = names(df_filtered),
  new = sub("_.*$", "", names(df_filtered))
)

# -----------------------------
# 5. Convert genotype strings to numeric
# -----------------------------
geno_to_num <- function(g) {
  if (g %in% c("0|0", "0/0")) return(0)
  else if (g %in% c("1|1", "1/1")) return(2)
  else return(1)
}

df_num <- as.data.frame(
  lapply(df_filtered, function(col) sapply(col, geno_to_num))
)

# -----------------------------
# 6. Transpose: genotypes as rows, SNPs as columns
# -----------------------------
df_num_t <- t(df_num)

rownames(df_num_t) <- colnames(df_filtered)
colnames(df_num_t) <- snp_ids

# -----------------------------
# 4. Hierarchical clustering on genotypes
# -----------------------------
d <- dist(df_num_t, method = "euclidean")
hc <- hclust(d, method = "average")

# -----------------------------
# 5. Cut tree into clusters
# -----------------------------
k <- 3   # number of clusters
groups <- cutree(hc, k = k)

# -----------------------------
# 6. Create cluster table
# -----------------------------
cluster_table <- data.frame(
  ivip = names(groups),
  Group = groups, row.names = NULL
)

# -----------------------------
# 7. Reorder dataframe by cluster order for TSV
# -----------------------------
ordered_genotypes <- hc$labels[hc$order]
df_ordered <- df_num_t[ordered_genotypes, ]
df_ordered <- cbind(Genotype = rownames(df_ordered), df_ordered)

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

summary(annotation)
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
         filename = "vfaba.whole_segment.heatmap_clusters_emmean.png",
         width = 10,
         height = 12)

#####
# Order rownames by emmean (low → high; use decreasing = TRUE if wanted)
row_order <- order(annotation$emmean)

df_num_t_ord <- df_num_t[row_order, ]
annotation_ord <- annotation[row_order, , drop = FALSE]

pheatmap(
  df_num_t_ord,
  cluster_rows = FALSE,               # IMPORTANT
  clustering_distance_cols = "euclidean",
  clustering_method = "average",
  cluster_cols = FALSE,
  annotation_row = annotation_ord,
  annotation_colors = annotation_colors,
  main = "Genotype similarity ordered by emmean",
  show_rownames = FALSE,
  filename = "vfaba.whole_segment.heatmap_emmean_ordered.png",
  width = 10,
  height = 12
)

#####


pheno <-read.csv("../../../../gitrepo/faba_specifics/02.analyses/01.GWAS/01_phenotipic_analysis/01.dataset_filtering/02_A.pisum/final filters/acyrthosiphon_pisum_filtered_6.csv")
summary(pheno)
merged_df <- merge(pheno, cluster_table, by = "ivip", all.x = TRUE)
merged_df$Group <- as.factor(merged_df$Group)

model <- lm(larvae_total_number ~ Group + ivip + session/block, data = merged_df)
anova(model)

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

summary(snp_haplotype_comparison)

model_snp <- lm(larvae_total_number ~ chr2_1461736473 + ivip + session/block, data = snp_haplotype_comparison)
anova(model_snp)
model_haplotype <- lm(larvae_total_number ~ Group + ivip + session/block, data = snp_haplotype_comparison)
anova(model_haplotype)
model_haplotype_seul <- lm(larvae_total_number ~ Group+ session/block, data = snp_haplotype_comparison)

model_snp_seul <- lm(larvae_total_number ~ chr2_1461736473 + session/block, data = snp_haplotype_comparison)

model_nul <- lm(larvae_total_number ~ ivip + session/block, data = snp_haplotype_comparison)

compare_performance(model_haplotype, model_snp, model_nul,model_haplotype_seul, model_snp_seul, rank = FALSE, verbose = TRUE)



model <- lmer(larvae_total_number ~ chr2_1461736473 + (1|ivip) + (1 | session/block), data = snp_haplotype_comparison)

emm_snp <- emmeans(model, ~ chr2_1461736473)
emm
emm_snp


dev.off()

# Basic boxplot using base R
boxplot(larvae_total_number ~ Group, 
        data = snp_haplotype_comparison, 
        main = "Boxplot of Larvae Total Number by Group",
        xlab = "Group",
        ylab = "Larvae Total Number",
        col = "grey",
        outline = FALSE)   # hide default outliers

stripchart(larvae_total_number ~ Group,
           data = snp_haplotype_comparison,
           vertical = TRUE,
           method = "jitter",
           pch = 16,
           col = "black",
           add = TRUE)

