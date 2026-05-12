############################################################

## LIBRARIES

############################################################
library(data.table)
library(dplyr)
library(tidyr)
library(pheatmap)
library(RColorBrewer)
library(vcfR)
library(snpStats)
library(LDheatmap)
library(tibble)

############################################################

## PARAMETERS

############################################################
chrM <- "chr2"
min_group_size <- 10

# small window

left_small  <- 1461736473
right_small <- 1461745735
# big window

left_big  <- 1461732669
right_big <- 1461746839

############################################################

## LOAD VCF

############################################################
vcf <- read.vcfR(
  "../../../../../../../files_too_heavy_for_git/pisum/vfaba_apisum_n6.filtered_imputed_maf.vcf.hedin2.vcf"
)

fix <- getFIX(vcf)

convertToNumeric <- function(x) {
  gdat <- matrix(NA, nrow = nrow(x), ncol = ncol(x))
  for (i in seq_len(nrow(x))) {
    for (j in seq_len(ncol(x))) {
      alleles <- unlist(strsplit(x[i, j], "[|/]"))
      gdat[i, j] <- sum(as.numeric(alleles))
    }
  }
  rownames(gdat) <- rownames(x)
  colnames(gdat) <- colnames(x)
  gdat
}

############################################################

## LD HEATMAP — BIG WINDOW (UNCHANGED)

############################################################
keep_init <- fix[, "CHROM"] == chrM &
  as.numeric(fix[, "POS"]) >= left_big &
  as.numeric(fix[, "POS"]) <= right_big

vcf_init <- vcf[keep_init, ]
gt_init <- t(vcf_init@gt[, -1])
gdat_init <- convertToNumeric(gt_init)

snp_pos_init <- as.numeric(vcf_init@fix[, "POS"])
snp_ids_init <- paste0(chrM, "_", snp_pos_init)

colnames(gdat_init) <- snp_ids_init
rownames(gdat_init) <- rownames(gt_init)

gdat_init <- as(gdat_init, "SnpMatrix")

custom_colors <- colorRampPalette(c("black","darkred","red","orange","yellow","white"))(100)

png("LD_heatmap_big_window.png", 1200, 1200, res = 150)
LDheatmap(gdat_init, color = custom_colors, genetic.distances = snp_pos_init, add.map=TRUE)
dev.off()

############################################################

## LD HEATMAP — SMALL WINDOW (UNCHANGED)

############################################################
keep_init <- fix[, "CHROM"] == chrM &
  as.numeric(fix[, "POS"]) >= left_small &
  as.numeric(fix[, "POS"]) <= right_small

vcf_init <- vcf[keep_init, ]
gt_init <- t(vcf_init@gt[, -1])
gdat_init <- convertToNumeric(gt_init)

snp_pos_init <- as.numeric(vcf_init@fix[, "POS"])
snp_ids_init <- paste0(chrM, "_", snp_pos_init)

colnames(gdat_init) <- snp_ids_init
rownames(gdat_init) <- rownames(gt_init)

gdat_init <- as(gdat_init, "SnpMatrix")

png("LD_heatmap_small_window.png", 1200, 1200, res = 150)
LDheatmap(gdat_init, color = custom_colors, genetic.distances = snp_pos_init, add.map=TRUE)
dev.off()

############################################################

## LOAD TSV

############################################################
df <- fread(
  "../../../../../../../files_too_heavy_for_git/pisum/vfaba_apisum_n6.filtered_imputed_maf.vcf.hedin2.tsv"
)

geno_to_num <- function(g) {
  if (g %in% c("0|0", "0/0")) return(0)
  else if (g %in% c("1|1", "1/1")) return(2)
  else return(1)
}

############################################################

## LOAD METADATA (ADDED)

############################################################
genet_cluster <- read.table(
  "../../../../../00.info/11.infos_genotypes/metadata.csv",
  header = TRUE, sep = ","
)
colnames(genet_cluster) <- c("ivip", "cluster")

emmeans_df <- read.csv(
  "../../../01_phenotipic_analysis/02.ajusted_means/pisum/emmeans_pisum_6.csv"
)
emmeans_df <- emmeans_df[, c("ivip","emmean")]

############################################################

## SMALL WINDOW MATRIX

############################################################
df_small <- df %>%
  separate(snp_id, into = c("chromosome", "position"), sep = "_", convert = TRUE) %>%
  filter(chromosome == chrM,
         position >= left_small,
         position <= right_small) %>%
  mutate(snp_id = paste0(chromosome, "_", position))

snp_ids_small <- df_small$snp_id

df_small <- df_small %>%
  select(-chromosome, -position, -snp_id)

df_small_num <- as.data.frame(lapply(df_small, function(x) sapply(x, geno_to_num)))
df_small_t <- t(df_small_num)
rownames(df_small_t) <- colnames(df_small)
colnames(df_small_t) <- snp_ids_small

############################################################

## HAPLOTYPE GROUPS (UNCHANGED LOGIC)

############################################################
haplotypes <- apply(df_small_t, 1, paste, collapse = "_")

df_groups <- data.frame(
  ivip = rownames(df_small_t),
  haplotype = haplotypes
)

df_groups_summary <- df_groups %>%
  group_by(haplotype) %>%
  summarise(group_size = n(), .groups = "drop") %>%
  filter(group_size >= min_group_size) %>%
  arrange(desc(group_size)) %>%
  mutate(group_id = paste0("Group_", row_number()))

df_groups_filtered <- df_groups %>%
  inner_join(df_groups_summary, by = "haplotype") %>%
  select(ivip, haplotype, group_size, group_id)

############################################################

## SMALL HEATMAP (WITH emmeans + cluster)

############################################################
df_small_selected <- df_small_t[df_groups_filtered$ivip, , drop = FALSE]

annotation1 <- df_groups_filtered %>%
  left_join(genet_cluster, by = "ivip") %>%
  left_join(emmeans_df, by = "ivip") %>%
  column_to_rownames("ivip")

annotation1$group_id <- as.factor(annotation1$group_id)
annotation1$cluster  <- as.factor(annotation1$cluster)
annotation1$emmean   <- as.numeric(annotation1$emmean)

annotation_colors <- list(
  group_id = setNames(
    brewer.pal(length(unique(annotation1$group_id)), "Set2"),
    unique(annotation1$group_id)
  ),
  cluster = setNames(
    brewer.pal(length(unique(annotation1$cluster)), "Set3"),
    unique(annotation1$cluster)
  ),
  emmean = colorRampPalette(c("white","black"))(100)
)

pheatmap(
  df_small_selected,
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "average",
  cluster_cols = FALSE,
  annotation_row = annotation1[, c("group_id","cluster","emmean")],
  annotation_colors = annotation_colors,
  show_rownames = FALSE,
  main = "Small window — groups ≥10 + emmean + cluster",
  filename = "pheatmap_small_groups.png",
  width = 10,
  height = 12
)

############################################################

## BIG WINDOW MATRIX

############################################################
df_big <- df %>%
  separate(snp_id, into = c("chromosome", "position"), sep = "_", convert = TRUE) %>%
  filter(chromosome == chrM,
         position >= left_big,
         position <= right_big) %>% mutate(snp_id = paste0(chromosome, "_", position))

snp_ids_big <- df_big$snp_id

df_big <- df_big %>%
  select(-chromosome, -position, -snp_id)

df_big_num <- as.data.frame(lapply(df_big, function(x) sapply(x, geno_to_num)))
df_big_t <- t(df_big_num)
rownames(df_big_t) <- colnames(df_big)
colnames(df_big_t) <- snp_ids_big

############################################################

## ORDERING (UNCHANGED)

############################################################
grouped_ivip <- df_groups_filtered$ivip
all_ivip <- rownames(df_big_t)
ungrouped_ivip <- setdiff(all_ivip, grouped_ivip)

group_levels <- df_groups_filtered %>%
  distinct(group_id) %>%
  arrange(as.numeric(sub("Group_", "", group_id))) %>%
  pull(group_id)

grouped_order <- c()

for (g in group_levels) {
  ivips_g <- df_groups_filtered %>%
    filter(group_id == g) %>%
    pull(ivip)
  
  if (length(ivips_g) > 2) {
    hc <- hclust(dist(df_big_t[ivips_g, ]))
    ivips_g <- hc$labels[hc$order]
  }
  
  grouped_order <- c(grouped_order, ivips_g)
}

if (length(ungrouped_ivip) > 2) {
  hc <- hclust(dist(df_big_t[ungrouped_ivip, ]))
  ungrouped_order <- hc$labels[hc$order]
} else {
  ungrouped_order <- ungrouped_ivip
}

final_order <- c(grouped_order, ungrouped_order)
df_big_ordered <- df_big_t[final_order, ]

############################################################

## BIG HEATMAP ANNOTATION (FIXED ORDER ISSUE)

############################################################
annotation2 <- data.frame(ivip = rownames(df_big_ordered))

annotation2 <- annotation2 %>%
  left_join(df_groups_filtered[, c("ivip","group_id")], by = "ivip") %>%
  left_join(genet_cluster, by = "ivip") %>%
  left_join(emmeans_df, by = "ivip")

annotation2$group_id[is.na(annotation2$group_id)] <- "Ungrouped"

rownames(annotation2) <- annotation2$ivip
annotation2$ivip <- NULL

annotation2$group_id <- as.factor(annotation2$group_id)
annotation2$cluster  <- as.factor(annotation2$cluster)
annotation2$emmean   <- as.numeric(annotation2$emmean)

############################################################

## FIX COLORS (IMPORTANT FIX)

############################################################
group_levels_big <- levels(annotation2$group_id)

annotation_colors_big <- list(
  group_id = setNames(
    colorRampPalette(brewer.pal(8, "Set2"))(length(group_levels_big)),
    group_levels_big
  ),
  cluster = setNames(
    brewer.pal(length(unique(annotation2$cluster)), "Set3"),
    unique(annotation2$cluster)
  ),
  emmean = colorRampPalette(c("white","black"))(100)
)

############################################################

## FINAL BIG HEATMAP

############################################################
pheatmap(
  df_big_ordered,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  annotation_row = annotation2[, c("group_id","cluster","emmean")],
  annotation_colors = annotation_colors_big,
  show_rownames = FALSE,
  main = "Big window — haplotypes + emmean + cluster",
  filename = "pheatmap_big_all.png",
  width = 10,
  height = 14
)
