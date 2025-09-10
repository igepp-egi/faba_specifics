#run in terminal : vcftools --vcf ../../../Downloads/testgwas/vfaba_hedin_peamust_SNP_GATK.vcf --missing-site --out snp_missing 

# Load the data
lmiss <- read.table("../../QC/scratch_29august2025/filtered/snp_missing.lmiss", header = TRUE)

# Calculate call rate
lmiss$CallRate <- 1 - lmiss$F_MISS

# Add SNP index for plotting (1, 2, 3, ...)
lmiss$SNP_Index <- 1:nrow(lmiss)

# Basic line plot
hist(lmiss$CallRate, breaks = 500)  # 500 bins

# Basic line plot
hist(lmiss$CallRate, breaks = 10)  # 30 bins


# Load the data
lmiss_after <- read.table("../../QC/scratch_29august2025/filtered/snp_missing_after.lmiss", header = TRUE)

# Calculate call rate
lmiss_after$CallRate <- 1 - lmiss_after$F_MISS

# Add SNP index for plotting (1, 2, 3, ...)
lmiss_after$SNP_Index <- 1:nrow(lmiss_after)

# Basic line plot
hist(lmiss_after$CallRate, breaks = 500)  # 500 bins

# Basic line plot
hist(lmiss_after$CallRate, breaks = 10)  # 30 bins






# run in terminal 

het <- read.table("../../files_GWAS/GWAS/heterozygosity_percentages.txt", header = FALSE)

#1
vcftools --vcf vfaba_hedin_peamust_SNP_GATK.vcf --hardy  
#2
GWAS awk 'NR > 1 {                             
split($3, a, "/");
total = a[1] + a[2] + a[3];
het_perc = (a[2] / total) * 100;
print $1, $2, het_perc
}' out.hwe > heterozygosity_percentages.txt

hist(het$V3, breaks = 500)  # 30 bins


