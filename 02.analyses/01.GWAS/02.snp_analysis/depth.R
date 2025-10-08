# vcftools --vcf vfaba.selected.biallelic.vcf --site-depth --out snp_site_depth
# vcftools --vcf vfaba.selected.biallelic.qc1.vcf --site-depth --out snp_site_depth

snp_depth_before <- read.table("../../QC/scratch_29august2025/filtered/snp_site_depth_before.ldepth", header = TRUE)
snp_depth <- read.table("../../QC/scratch_29august2025/filtered/snp_site_depth.ldepth", header = TRUE)

snp_depth$average <- snp_depth$SUM_DEPTH / 248

hist(snp_depth$average, xlim = c(0, 750), breaks = 100)


snp_depth_before$average <- snp_depth_befor$SUM_DEPTH / 248

hist(snp_depth_before$average, xlim = c(0, 750), breaks = 100)


snp_quality_before <- read.table("../../QC/scratch_29august2025/filtered/snp_qual.tsv", header = FALSE)

hist(snp_quality_before$V3, xlim = c(0, 120000), breaks = 500)


snp_quality_before <- read.table("../../QC/scratch_29august2025/filtered/snp_qual_after.tsv", header = FALSE)

hist(snp_quality_before$V3, xlim = c(0, 120000), breaks = 500)
