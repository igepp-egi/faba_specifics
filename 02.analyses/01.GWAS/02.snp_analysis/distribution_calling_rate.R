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

#per genotype:
vcftools --vcf ../../../Downloads/testgwas/vfaba_hedin_peamust_SNP_GATK.vcf --missing-indv --out genotype_missing



# run in terminal 

#1
vcftools --vcf vfaba_hedin_peamust_SNP_GATK.vcf --hardy  
#2
GWAS awk 'NR > 1 {                             
split($3, a, "/");
total = a[1] + a[2] + a[3];
het_perc = (a[2] / total) * 100;
print $1, $2, het_perc
}' out.hwe > heterozygosity_percentages.txt

het <- read.table("../../files_GWAS/GWAS/heterozygosity_percentages.txt", header = FALSE)

hist(het$V3, breaks = 500)  # 30 bins


# in terminal
vcftools --vcf vfaba.selected.biallelic.qc1.qc2.maf.miss.vcf --het

awk 'NR==1 { print $0, "H_obs" } NR>1 { h=1-($2/$4); print $0, h }' out.het | \
awk 'NR==1 || $6 > 0.1'

het2 <- read.table("../../QC/scratch/filtered/out.het", header = TRUE)
summary(het2)summary(het2)TRUE
het2$O.HOM.<- as.numeric(het2$O.HOM.)
het2$N_SITES<- as.numeric(het2$N_SITES)
het2$het <- 1-(het2$O.HOM./het2$N_SITES)

hist(het2$het, breaks = 500)
