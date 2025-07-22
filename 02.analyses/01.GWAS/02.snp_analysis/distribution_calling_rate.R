#run in terminal : vcftools --vcf ../../../Downloads/testgwas/vfaba_hedin_peamust_SNP_GATK.vcf --missing-site --out snp_missing 

# Load the data
lmiss <- read.table("../../files_GWAS/GWAS/snp_missing.lmiss", header = TRUE)

# Calculate call rate
lmiss$CallRate <- 1 - lmiss$F_MISS

# Add SNP index for plotting (1, 2, 3, ...)
lmiss$SNP_Index <- 1:nrow(lmiss)

# Basic line plot
hist(lmiss$CallRate, breaks = 500)  # 500 bins

# Basic line plot
hist(lmiss$CallRate, breaks = 10)  # 30 bins

