# Steps

## All commands

vcftools --gzvcf vfaba_hedin_peamust_SNP_GATK.vcf.gz --remove-indels \
         --min-alleles 2 --max-alleles 2 --recode --out vfaba_noindels_bisnps.vcf &&

vcftools --vcf vfaba_noindels_bisnps.vcf --maf 0.01 --max-missing 0.9 \
         --recode --recode-INFO-all \
         --out vfaba_noindels_bisnps_maf_miss 

vcftools --vcf vfaba_noindels_bisnps_maf_miss --het --out heterozygosity_stats 

## Raw file

zcat vfaba_hedin_peamust_SNP_GATK.vcf.gz |grep -v "#" |wc -l

There are 6,425,454 SNPs in total in the original raw file.

## Only keeping biallelic SNPs

vcftools --gzvcf vfaba_hedin_peamust_SNP_GATK.vcf.gz --remove-indels --min-alleles 2 --max-alleles 2 --recode --out vfaba_noindels_bisnps.vcf

cat vfaba_noindels_bisnps.vcf |grep -v "#" |wc -l

There are 2,776,935 biallelic SNPs (no indels left).

## Filter max MAF 1% and missing data max. 10%

Goal: to reproduce **PEAMUST_vfaba_F1b.vcf.gz** (SNPs + INDELS with heterozygosity<15%, MAF<1% and less than 10% missing data)

Exclude sites on the basis of the proportion of missing data
(defined to be between 0 and 1, where 0 allows sites that are completely missing and 1 indicates no missing data allowed).

vcftools --vcf vfaba_noindels_bisnps.vcf --maf 0.01 --max-missing 0.9 --out vfaba_noindels_bisnps_maf_miss --recode --recode-INFO-all


## Filter base on heterozygosity

vcftools --vcf vfaba_noindels_bisnps_maf_miss --het --out heterozygosity_stats 





