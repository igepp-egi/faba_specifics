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

# small window 1461745735 1461745994

left_small  <- 1461736473
right_small <- 1461745994

# big window

left_big  <- 1461732669
right_big <- 1461746839

############################################################

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


df_small <- df %>%
  separate(snp_id, into = c("chromosome", "position"), sep = "_", convert = TRUE) %>%
  filter(chromosome == chrM,
         position >= left_small,
         position <= right_small) %>%
  mutate(snp_id = paste0(chromosome, "_", position))

snp_ids_small <- df_small$snp_id

df_small <- df_small %>%
  dplyr::select(-chromosome, -position, -snp_id)

df_small_num <- as.data.frame(lapply(df_small, function(x) sapply(x, geno_to_num)))
df_small_t <- t(df_small_num)
rownames(df_small_t) <- colnames(df_small)
colnames(df_small_t) <- snp_ids_small

############################################################

## HAPLOTYPE GROUPS (UNCHANGED LOGIC)

############################################################
haplotypes <- apply(df_small_t, 1, paste, collapse = "_")

############################################################
## HAPLOTYPE GROUPING (1 == 2 LOGIC)
############################################################

# --- STEP 1: collapse genotypes treating 1 as 2 ---
df_binary <- df_small_t
df_binary[df_binary == 1] <- 2   # IMPORTANT: merge 1 into 2

haplotypes <- apply(df_binary, 1, paste, collapse = "_")

############################################################
## STEP 1 — ORIGINAL EXACT MATCH GROUPS (haplotype BASELINE)
############################################################

hap_raw <- apply(df_small_t, 1, paste, collapse = "_")

df_groups <- data.frame(
  ivip = rownames(df_small_t),
  haplotype = hap_raw
)

df_groups_summary <- df_groups %>%
  group_by(haplotype) %>%
  summarise(
    group_size = n(),
    .groups = "drop"
  ) %>%
  filter(group_size >= min_group_size) %>%
  mutate(group_id = paste0("haplotype_", row_number()))

df_groups_haplotype <- df_groups %>%
  inner_join(df_groups_summary, by = "haplotype")

############################################################
## STEP 2 — FIND UNGROUPED IVIPs
############################################################

ungrouped_ivip <- setdiff(rownames(df_small_t), df_groups_haplotype$ivip)

df_ungrouped <- df_small_t[ungrouped_ivip, , drop = FALSE]

############################################################
## STEP 3 — RE-GROUP UNGROUPED (1 == 2 LOGIC)
############################################################

df_binary <- df_ungrouped
df_binary[df_binary == 1] <- 2

hap_binary <- apply(df_binary, 1, paste, collapse = "_")

df_groups2 <- data.frame(
  ivip = rownames(df_ungrouped),
  haplotype = hap_binary
)

df_groups2_summary <- df_groups2 %>%
  group_by(haplotype) %>%
  summarise(
    group_size = n(),
    .groups = "drop"
  ) %>%
  filter(group_size >= min_group_size) %>%
  mutate(group_id = paste0("Rescued_", row_number()))

df_groups_rescued <- df_groups2 %>%
  inner_join(df_groups2_summary, by = "haplotype")


############################################################
## STEP 4 — MERGE haplotype + RESCUED GROUPS
############################################################

df_groups_filtered <- bind_rows(
  df_groups_haplotype,
  df_groups_rescued
) %>%
  dplyr::select(ivip, haplotype, group_id)

df_groups_filtered$group_id[df_groups_filtered$group_id == "Rescued_1"] <- "haplotype_6"
############################################################




############################################################
## SMALL HEATMAP (GROUP ONLY — CLEAN)
############################################################
df_small_selected <- df_small_t[df_groups_filtered$ivip, , drop = FALSE]

annotation1 <- df_groups_filtered %>%
  column_to_rownames("ivip")

annotation1$group_id <- as.factor(annotation1$group_id)

############################################################
## DEFINE ONE SHARED COLOR MAP (USED EVERYWHERE)
############################################################
library(RColorBrewer)

group_levels_all <- sort(unique(annotation1$group_id))

group_colors <- setNames(
  colorRampPalette(brewer.pal(8, "Set2"))(length(group_levels_all)),
  group_levels_all
)

annotation_colors <- list(
  group_id = group_colors
)

############################################################
## SMALL HEATMAP
############################################################
pheatmap(
  df_small_selected,
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "average",
  cluster_cols = FALSE,
  annotation_row = annotation1["group_id"],
  annotation_colors = annotation_colors,
  show_rownames = FALSE,
  main = "Small window — groups ≥10",
  filename = "pheatmap_small_groups.png",
  width = 10,
  height = 12
)


############################################################
df_big <- df %>%
  separate(snp_id, into = c("chromosome", "position"), sep = "_", convert = TRUE) %>%
  filter(chromosome == chrM,
         position >= left_big,
         position <= right_big) %>% mutate(snp_id = paste0(chromosome, "_", position))

snp_ids_big <- df_big$snp_id

df_big <- df_big %>%
  dplyr::select(-chromosome, -position, -snp_id)

df_big_num <- as.data.frame(lapply(df_big, function(x) sapply(x, geno_to_num)))
df_big_t <- t(df_big_num)
rownames(df_big_t) <- colnames(df_big)
colnames(df_big_t) <- snp_ids_big


############################################################
## MANUAL GROUP ORDER 
############################################################

group_levels <- c(
  "haplotype_1",
  "haplotype_2",
  "haplotype_3",
  "haplotype_4",
  "haplotype_5",
  "haplotype_6"
)

# safety: keep only groups that actually exist
group_levels <- group_levels[group_levels %in% unique(df_groups_filtered$group_id)]

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

final_order <- unique(c(grouped_order, ungrouped_order))
df_big_ordered <- df_big_t[final_order, ]


############################################################
## BIG HEATMAP ANNOTATION (GROUP ONLY)
############################################################
annotation2 <- data.frame(ivip = rownames(df_big_ordered))

annotation2 <- annotation2 %>%
  left_join(df_groups_filtered[, c("ivip","group_id")], by = "ivip")



annotation2$group_id[is.na(annotation2$group_id)] <- "ungrouped"

rownames(annotation2) <- annotation2$ivip
annotation2$ivip <- NULL

annotation2$group_id <- as.factor(annotation2$group_id)

############################################################
## UPDATE COLORS TO INCLUDE UNGROUPED
############################################################
all_groups_big <- sort(unique(annotation2$group_id))

group_colors_big <- setNames(
  colorRampPalette(brewer.pal(8, "Set2"))(length(all_groups_big)),
  all_groups_big
)

annotation_colors_big <- list(
  group_id = group_colors_big
)

############################################################
## FINAL BIG HEATMAP (GROUP ONLY)
############################################################
pheatmap(
  df_big_ordered,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  annotation_row = annotation2["group_id"],
  annotation_colors = annotation_colors_big,
  show_rownames = FALSE,
  annotation_names_row = TRUE,
  fontsize = 20,
  fontsize_row = 20,
  fontsize_col = 20,
  fontsize_number = 20,
  main = "Big window — haplotypes",
  filename = "pheatmap_big_all.png",
  width = 13,
  height = 14
)

############################################################
## BOXPLOT WITH SAME COLORS
############################################################

pheno <-read.csv("../../../../01.GWAS/01_phenotipic_analysis/01.dataset_filtering/02_A.pisum/final filters/acyrthosiphon_pisum_filtered_6.csv") 


annotation2 <- annotation2 %>%
  rownames_to_column(var = "ivip")
merged_df <- merge(pheno, annotation2, by = "ivip", all.x = TRUE)

# IMPORTANT: ensure same order
merged_df$group_id <- factor(
  merged_df$group_id,
  levels = names(group_colors_big)
)
png( filename = "boxplot_haplotype_comparison_chr2_1461736473.png", width = 3800, height = 1800, res = 300 )

boxplot(larvae_total_number ~ group_id,
        data = merged_df,
        main = "Boxplot of Larvae Total Number by Group",
        xlab = "Group",
        ylab = "Larvae Total Number",
        col = group_colors_big,
        cex.main = 1,    # Title size (1 = default, 2 = 2x larger)
        cex.lab = 1.5,   # Axis label size
        cex.axis = 1.4   # Axis tick label size
)
dev.off()




favorable_haps <- c("haplotype_6")

annotation2$group_id <- ifelse(
  annotation2$group_id %in% favorable_haps,
  "Favorable",
  "Unfavorable"
)

names(annotation2)[2] <- "chr2_1461736473"


write.csv(
  annotation2,
  "haplo_chr2_1461736473.csv",
  row.names = FALSE
)



