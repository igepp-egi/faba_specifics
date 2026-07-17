## LIBRARIES ----

library(data.table)
library(dplyr)
library(tidyr)
library(pheatmap)
library(RColorBrewer)
library(vcfR)
library(snpStats)
library(LDheatmap)
library(tibble)


## LD HEATMAP — BIG WINDOW----
chrM <- "chr2"
min_group_size <- 10

# big window

left_big  <- 1461732669
right_big <- 1461746839

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

## LD HEATMAP — SMALL WINDOW----
# based on LD:
left_small  <- 1461736473
right_small <- 1461745994


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

## LOAD TSV ----

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
## PHASED HAPLOTYPE ANALYSIS ----

df_phase <- as.matrix(df_small)
rownames(df_phase) <- snp_ids_small
df_phase <- t(df_phase)


## SPLIT INTO THE TWO PHASED HAPLOTYPES

hapA <- matrix(
  NA,
  nrow = nrow(df_phase),
  ncol = ncol(df_phase),
  dimnames = dimnames(df_phase)
)

hapB <- hapA

for(i in seq_len(nrow(df_phase))){
  for(j in seq_len(ncol(df_phase))){
    
    alleles <- strsplit(df_phase[i,j], "\\|")[[1]]
    
    hapA[i,j] <- alleles[1]
    hapB[i,j] <- alleles[2]
  }
}


## CREATE CHROMOSOME HAPLOTYPES


hapA_string <- apply(hapA, 1, paste, collapse = "_")
hapB_string <- apply(hapB, 1, paste, collapse = "_")


## FIND MOST COMMON HAPLOTYPE STARTING WITH 1


all_haps <- data.frame(
  haplotype = c(hapA_string, hapB_string)
)

reference_haplotype <- all_haps %>%
  filter(substr(haplotype,1,1)=="1") %>%
  count(haplotype, sort=TRUE) %>%
  slice(1) %>%
  pull(haplotype)

reference_vector <- strsplit(reference_haplotype,"_")[[1]]


## CLASSIFY INDIVIDUALS

df_groups_filtered <- data.frame(
  ivip = rownames(df_phase),
  group_id = "other",
  stringsAsFactors = FALSE
)

for(i in seq_len(nrow(df_phase))){
  
  matchA <- all(hapA[i,] == reference_vector)
  matchB <- all(hapB[i,] == reference_vector)
  
  if(matchA & matchB){
    
    df_groups_filtered$group_id[i] <- "homozygous"
    
  } else if(xor(matchA, matchB)){
    
    df_groups_filtered$group_id[i] <- "heterozygous"
    
  }
  
}



pheno <- read.csv("../../../../01.GWAS/01_phenotipic_analysis/01.dataset_filtering/02_A.pisum/final filters/acyrthosiphon_pisum_filtered_6.csv", sep =",")

merge<-merge(df_groups_filtered,pheno, by ="ivip", all.x = TRUE)
merge$group_id <- as.factor(merge$group_id)


# For pisum
mod <- lm(larvae_total_number ~ group_id + session/block, data = merge)
anova <- anova(mod)
print(anova)

# Run Tukey HSD
tukey_pisum <- TukeyHSD(aov(larvae_total_number ~ group_id, data = merge))
print(tukey_pisum)


emmean<-read.csv("../../../../01.GWAS/01_phenotipic_analysis/02.ajusted_means/pisum/emmeans_pisum_6.csv", sep =",")
merge_emmean<-merge(df_groups_filtered,emmean, by ="ivip", all.x = TRUE)
merge_emmean$groups_id <- as.factor(merge_emmean$group_id)

boxplot(
  emmean ~ group_id,
  data = merge_emmean,
  main = "Boxplot of aphid count by Group",
  xlab = "Group",
  ylab = "aphid count",
  col = "grey"
)


## HAPLOTYPE ANALYSIS ----

## LOAD BIG WINDOW

df_big <- df %>%
  separate(snp_id, into = c("chromosome", "position"),
           sep = "_", convert = TRUE) %>%
  filter(
    chromosome == chrM,
    position >= left_big,
    position <= right_big
  ) %>%
  mutate(snp_id = paste0(chromosome, "_", position))

snp_ids_big <- df_big$snp_id

df_big <- df_big %>%
  dplyr::select(-chromosome, -position, -snp_id)

df_big_num <- as.data.frame(
  lapply(df_big, function(x) sapply(x, geno_to_num))
)

df_big_t <- t(df_big_num)

rownames(df_big_t) <- colnames(df_big)
colnames(df_big_t) <- snp_ids_big



## CREATE HAPLOTYPE GROUPS


## Exact haplotypes
hap_raw <- apply(df_small_t, 1, paste, collapse = "_")

df_groups_joint <- data.frame(
  ivip = rownames(df_small_t),
  haplotype = hap_raw,
  stringsAsFactors = FALSE
)


## FIND REFERENCE (MOST COMMON HAPLOTYPE STARTING WITH 2)


hap_summary <- df_groups_joint %>%
  count(haplotype, sort = TRUE) %>%
  mutate(first_snp = sub("_.*", "", haplotype))

reference_haplotype <- hap_summary %>%
  filter(first_snp == "2") %>%
  slice(1) %>%
  pull(haplotype)

reference_vector <- strsplit(reference_haplotype, "_")[[1]]


## COLLAPSE HOMOZYGOUS + HETEROZYGOUS

for(i in seq_len(nrow(df_groups_joint))){
  
  hap <- strsplit(df_groups_joint$haplotype[i], "_")[[1]]
  
  ## Only haplotypes beginning with 2 can belong
  ## to the reference group
  if(hap[1] == "2"){
    
    hap_test <- hap
    hap_test[hap_test == "1"] <- "2"
    
    if(all(hap_test == reference_vector)){
      df_groups_joint$haplotype[i] <- reference_haplotype
    }
  }
  
}


## CREATE GROUPS


df_groups_joint_summary <- df_groups_joint %>%
  group_by(haplotype) %>%
  summarise(
    group_size = n(),
    .groups = "drop"
  )

## Common haplotypes
df_groups_joint_summary <- df_groups_joint_summary %>%
  filter(group_size >= min_group_size) %>%
  arrange(desc(group_size)) %>%
  mutate(group_id = paste0("haplotype_", row_number()))

df_groups_joint_filtered <- df_groups_joint %>%
  left_join(df_groups_joint_summary, by = "haplotype")

## Pool rare haplotypes
df_groups_joint_filtered$group_id[is.na(df_groups_joint_filtered$group_id)] <- "other"

df_groups_joint_filtered <- df_groups_joint_filtered %>%
  dplyr::select(ivip, haplotype, group_id)

## RENAME HAPLOTYPES

df_groups_joint_filtered$group_id <- recode(
  df_groups_joint_filtered$group_id,
  "haplotype_4" = "haplotype_1",
  "haplotype_6" = "haplotype_4",
  "haplotype_1" = "haplotype_6")


## DEFINE PHEATMAP GROUP ORDER


group_levels <- c(
  "haplotype_1",
  "haplotype_2",
  "haplotype_3",
  "haplotype_4",
  "haplotype_5",
  "haplotype_6",
  "other"
)

# keep only groups present
group_levels <- group_levels[
  group_levels %in% unique(df_groups_joint_filtered$group_id)
]



## ORDER INDIVIDUALS BY HAPLOTYPE GROUP


ordered_ivips <- character(0)

for(g in group_levels){
  
  ivips_g <- df_groups_joint_filtered %>%
    filter(group_id == g) %>%
    pull(ivip)
  
  # cluster individuals inside each haplotype
  if(length(ivips_g) > 2){
    
    hc <- hclust(dist(df_big_t[ivips_g, ]))
    ivips_g <- hc$labels[hc$order]
    
  }
  
  ordered_ivips <- c(ordered_ivips, ivips_g)
  
}

df_big_ordered <- df_big_t[ordered_ivips, ]

## HEATMAP ANNOTATION


annotation2 <- df_groups_joint_filtered %>%
  dplyr::select(ivip, group_id)

rownames(annotation2) <- annotation2$ivip
annotation2$ivip <- NULL

annotation2$group_id <- factor(
  annotation2$group_id,
  levels = group_levels
)


## COLORS


group_colors_big <- setNames(
  colorRampPalette(brewer.pal(8, "Set2"))(length(group_levels)),
  group_levels
)

annotation_colors_big <- list(
  group_id = group_colors_big
)


## BIG HEATMAP


pheatmap(
  df_big_ordered,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  annotation_row = annotation2["group_id"],
  annotation_colors = annotation_colors_big,
  show_rownames = FALSE,
  show_colnames = FALSE,
  annotation_names_row = FALSE,   # remove "group_id"
  annotation_legend = FALSE,      # remove annotation legend
  legend = FALSE,                 # remove heatmap legend
  fontsize = 20,
  fontsize_row = 20,
  fontsize_col = 20,
  fontsize_number = 20,
  filename = "pheatmap_big_all.png",
  width = 13,
  height = 14
)


## BOXPLOT


pheno <- read.csv(
  "../../../../01.GWAS/01_phenotipic_analysis/01.dataset_filtering/02_A.pisum/final filters/acyrthosiphon_pisum_filtered_6.csv"
)

annotation_box <- annotation2 %>%
  tibble::rownames_to_column("ivip")

merged_df <- merge(
  pheno,
  annotation_box,
  by = "ivip",
  all.x = TRUE
)

merged_df$group_id <- factor(
  merged_df$group_id,
  levels = group_levels
)

postscript(
  file = "boxplot_haplotype_comparison.eps",
  width = 12.7,
  height = 6,
  paper = "special",
  horizontal = FALSE)

boxplot(
  larvae_total_number ~ group_id,
  data = merged_df,
  main = "Larvae Total Number by Haplotype Group",
  xlab = "Group",
  ylab = "Larvae Total Number",
  col = group_colors_big,
  cex.main = 1,
  cex.lab = 1.5,
  cex.axis = 1.4
)

dev.off()


## EXPORT HAPLOTYPE ASSIGNMENTS


annotation_export <- annotation_box

names(annotation_export)[2] <- paste0(chrM, "_", left_small)

write.csv(
  annotation_export,
  paste0("haplo_", chrM, "_", left_small, ".csv"),
  row.names = FALSE
)

