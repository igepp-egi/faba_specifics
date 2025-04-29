# Verify numbers of SNPs after each filter

## Summary 

**Goal:** to recreate **PEAMUST_vfaba_F1b.vcf.gz** (SNPs + INDELS with heterozygosity<15%, MAF<1% and less than 10% missing data)

The co
Filters: 
- filter 0: 6,425,454 SNPs in the raw `vfaba_hedin_peamust_SNP_GATK.vcf.gz` VCF.
- filter 1: 5,819,008 biallelic SNPs.
- filter 2: 4,007,444 biallelic SNPs with a MAF higher or equal to 1%.
- filter 3: 1,819,778 biallelic SNPs with a MAF higher than 1% and less than 10% missing data.  
- filter 4: 1,498,907 biallelic SNPs with a MAF higher than 1% and less than 10% missing data and less than 15% heterozygosity.

## Complete script check_numbers.sh

```bash
MY_VCF="vfaba_hedin_peamust_SNP_GATK.vcf.gz"

bcftools view --max-alleles 2 --with-header -Oz --output vcf_biallelic.vcf.gz $MY_VCF
bcftools view -i 'MAF[0]> 0.01' --threads 5 -Oz -o vcf_biallelic_maf001.vcf.gz vcf_biallelic.vcf.gz
bcftools view -i 'F_MISSING < 0.1' -Oz -o vcf_biallelic_maf001_miss.vcf.gz --threads 5 vcf_biallelic_maf001.vcf.gz
bcftools query -f '%CHROM\t%POS\t[%GT\t]\n' vcf_biallelic_maf001_miss.vcf.gz > genotypes.txt

awk -F '\t' '{
    split($3, genotypes, /[|\/]/);
    count_hom = 0;
    count_het = 0;
    for (i in genotypes) {
        if (genotypes[i] != "." && genotypes[i] == genotypes[i+1]) {
            count_hom++;
        } else if (genotypes[i] != "." && genotypes[i] != genotypes[i+1]) {
            count_het++;
        }
    }
    total = count_hom + count_het;
    if (total > 0) {
        heterozygosity = count_het / total;
        print $1 "\t" $2 "\t" heterozygosity;
    }
}' genotypes.txt > heterozygosity.txt

awk '$3 < 0.15 {print $1 "\t" $2}' heterozygosity.txt > sites_to_keep_less15pc_het.txt

bcftools view --targets-file sites_to_keep_less15pc_het.txt -Oz -o vcf_biallelic_maf001_miss_maxhet.vcf.gz --threads 5 vcf_biallelic_maf001_miss.vcf.gz
initial_snp_counts=$(zcat vfaba_hedin_peamust_SNP_GATK.vcf.gz           | grep -v "#" |wc -l)
after_filter1_snp_counts=$(zcat vcf_biallelic.vcf.gz                    | grep -v "#" |wc -l)
after_filter2_snp_counts=$(zcat vcf_biallelic_maf001.vcf.gz             | grep -v "#" |wc -l)
after_filter3_snp_counts=$(zcat vcf_biallelic_maf001_miss.vcf.gz        | grep -v "#" |wc -l)
after_filter4_snp_counts=$(zcat vcf_biallelic_maf001_miss_maxhet.vcf.gz | grep -v "#" |wc -l)

echo "The raw $MY_VCF has $initial_snp_counts lines."
echo "There are $after_filter1_snp_counts biallelic SNPs."
echo "There are $after_filter2_snp_counts biallelic SNPs with a MAF higher or equal to 1%."
echo "There are $after_filter3_snp_counts biallelic SNPs with a MAF higher than 1% and less than 10% missing data."
echo "There are $after_filter4_snp_counts biallelic SNPs with a MAF higher than 1% and less than 10% missing data and less than 15% heterozygosity."

```
## Filter 0: raw file

zcat vfaba_hedin_peamust_SNP_GATK.vcf.gz |grep -v "#" |wc -l

There are 6,425,454 SNPs in total in the original raw file.

## Filter 1: keep only biallelic SNPs (remove indels and multiallelic SNPs) 

```bash
bcftools view --max-alleles 2 --with-header -Oz --output vcf_biallelic.vcf.gz $MY_VCF
```

## Filter 2: max MAF 1% 

Exclude sites on the basis of the proportion of missing data. 
(defined to be between 0 and 1, where 0 allows sites that are completely missing and 1 indicates no missing data allowed).

```bash
bcftools view -i 'MAF[0]> 0.01' --threads 5 -Oz -o vcf_biallelic_maf001.vcf.gz vcf_biallelic.vcf.gz

zcat vcf_biallelic_maf001.vcf.gz |grep -v "#" |wc -l
```

## Filter 3: max missing data max. 10%

```bash
bcftools view -i 'F_MISSING < 0.1' -Oz -o vcf_biallelic_maf001_miss.vcf.gz --threads 5 vcf_biallelic_maf001.vcf.gz

zcat vcf_biallelic_maf001_miss.vcf.gz |grep -v "#" |wc -l
```

## Filter 4: filter based on heterozygosity

First calculate heterozygosity: 
```bash
bcftools query -f '%CHROM\t%POS\t[%GT\t]\n' vcf_biallelic_maf001_miss.vcf.gz > genotypes.txt

awk -F '\t' '{
    split($3, genotypes, /[|\/]/);
    count_hom = 0;
    count_het = 0;
    for (i in genotypes) {
        if (genotypes[i] != "." && genotypes[i] == genotypes[i+1]) {
            count_hom++;
        } else if (genotypes[i] != "." && genotypes[i] != genotypes[i+1]) {
            count_het++;
        }
    }
    total = count_hom + count_het;
    if (total > 0) {
        heterozygosity = count_het / total;
        print $1 "\t" $2 "\t" heterozygosity;
    }
}' genotypes.txt > heterozygosity.txt
```

Select variants with less than 15% heterozygosity. 
```bash
bcftools view --targets-file sites_to_keep_less15pc_het.txt -Oz -o vcf_biallelic_maf001_miss_maxhet.vcf.gz --threads 5 vcf_biallelic_maf001_miss.vcf.gz
```

## Complete script 

```bash
initial_snp_counts=$(zcat vfaba_hedin_peamust_SNP_GATK.vcf.gz           | grep -v "#" |wc -l)
after_filter1_snp_counts=$(zcat vcf_biallelic.vcf.gz                    | grep -v "#" |wc -l)
after_filter2_snp_counts=$(zcat vcf_biallelic_maf001.vcf.gz             | grep -v "#" |wc -l)
after_filter3_snp_counts=$(zcat vcf_biallelic_maf001_miss.vcf.gz        | grep -v "#" |wc -l)
after_filter4_snp_counts=$(zcat vcf_biallelic_maf001_miss_maxhet.vcf.gz | grep -v "#" |wc -l)

echo "The raw $MY_VCF has $initial_snp_counts lines."
echo "There are $after_filter1_snp_counts biallelic SNPs."
echo "There are $after_filter2_snp_counts biallelic SNPs with a MAF higher or equal to 1%."
echo "There are $after_filter3_snp_counts biallelic SNPs with a MAF higher than 1% and less than 10% missing data."
echo "There are $after_filter4_snp_counts biallelic SNPs with a MAF higher than 1% and less than 10% missing data and less than 15% heterozygosity."
```

