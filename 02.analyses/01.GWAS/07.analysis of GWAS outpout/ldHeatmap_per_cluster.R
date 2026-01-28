## Load required libraries
library(vcfR)
library(LDheatmap)


## -----------------------------
## Read VCF
## -----------------------------
vcf <- read.vcfR("vfaba_apisum_n6.filtered_imputed_maf.vcf.hedin2.vcf")

## -----------------------------
## Subset to genomic region
## -----------------------------
fix  <- getFIX(vcf)

keep <- fix[, "CHROM"] == "chr2" &
  as.numeric(fix[, "POS"]) >= 1461636473 &
  as.numeric(fix[, "POS"]) <= 1461836473

snp <- vcf[keep, ]

cluster_info <- read.csv("../../../../gitrepo/faba_specifics/00.info/00.passeport_data/faba_bean_ivip_and_genetic_cluster.csv")

target_cluster <- "South"
cluster_samples <- cluster_info$ivip[cluster_info$cluster_name == target_cluster]


library(vcfR)

gt_matrix <- snp@gt[, -1]  # remove FORMAT column

# Only keep samples present in both VCF and cluster
samples_in_vcf <- colnames(gt_matrix)
cluster_samples_in_vcf <- intersect(cluster_samples, samples_in_vcf)

# Subset genotype matrix
gt_cluster <- gt_matrix[, cluster_samples_in_vcf]

# Transpose: samples = rows, SNPs = columns
snpMat_cluster <- t(gt_cluster)



convertToNumeric <- function(x){
  gdat <- matrix(NA, nrow = nrow(x), ncol = ncol(x))
  for (m in 1:nrow(x)){
    for (n in 1:ncol(x)){
      alleles <- unlist(strsplit(x[m,n], "[|/]"))
      gdat[m,n] <- sum(as.numeric(alleles))
    }
  }
  rownames(gdat) <- rownames(x)
  colnames(gdat) <- colnames(x)
  return(gdat)
}

gdat_cluster <- convertToNumeric(snpMat_cluster)


snp_pos <- as.numeric(snp@fix[,2])  # POS
colnames(gdat_cluster) <- snp@fix[,3]  # SNP IDs


library(snpStats)
library(LDheatmap)

gdat_cluster <- as(gdat_cluster, "SnpMatrix")
n_colors <- 100
low_range <- seq(0, 0.2, length.out = 90)  # 90 colors for 0-0.2
high_range <- seq(0.21, 1, length.out = 10) # 10 colors for 0.2-1
color_breaks <- sort(c(low_range, high_range))

# Define the colors: yellow to red for 0-0.2, red to dark red for 0.2-1
color_list <- c( "black","darkred", "red","orange","yellow","white")

# Generate the custom color ramp
custom_colors <- colorRampPalette(color_list)(n_colors)

LDheatmap(gdat_cluster, snp_pos, color = custom_colors, title=paste("Cluster", target_cluster), add.map=FALSE)
