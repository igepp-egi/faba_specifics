# README

This README has 4 sections:   
1. Check SNP numbers  
2. Information from Jonathan on his favorite steps  
3. VCF file description  
4. Complete header of PeaMUST Vicia faba VCF file  

<!-- MarkdownTOC autolink="true" levels="1,2" -->

- [Check SNP numbers](#check-snp-numbers)
- [Information from Jonathan on his favorite steps on Pisum sativum](#information-from-jonathan-on-his-favorite-steps-on-pisum-sativum)

<!-- /MarkdownTOC -->

# Check SNP numbers 

To verify the numbers in [the different VCF files provided](./snp_different_filters.png). 
A script called `check_numbers.sh` has been created (see below for its content): 

The related conda environment has to be created first and actived using: 
```bash
conda create --name vcf -c conda-forge -c bioconda bcftools && conda activate vcf
```
Followed by execution of the script using: `bash check_numbers.sh`.

See [check_numbers.md](./check_numbers.md) for details. 

# Information from Jonathan on his favorite steps on Pisum sativum

First, we filtrate our whole matrix obtained with GATK to keep only 240 genotypes from AMS. Due to those filters, we will maybe see some invariants sites or with only missing data. Use a really low MAF/high missing data to filtrate them now if needed. If not will be done later. 

```bash
grep -f AMS_theo.txt DonneesPasseport_AMS_GRaSP_AK090523.csv  | cut -f1 IPIP_to_keep.lst
bcftools view -S IPIP_to_keep.lst -o AMS_v2_raw.vcf.gz -O  z cameor_v2_r3_peamust_genotyped_SNP_hardfilter.vcf.gz
``` 

I use a really bad custom script to compute a few stats and do some figures. We use one of the file produced to filter site with more than 15% of Htz. This is arbitrary and based on what was done by **R.Ollivier**. I think it works well with GWAS but it could be impossible to detect SNPs in balancing selection in our largest subpopulation. 

```bash
Rscript --vanilla /GEAPSI/peamust/FABA_2022/Data/vcf_stats.R AMS_v2_raw.vcf.gz False AMS_v2_raw
awk 'NR>=2 {if( $5 < 240 && ($4)/($2 + $3 + $4 ) <= 0.15){print $0}}' AMS_v2_raw_bySite.csv  |  cut -f1 > SNP_to_keep.csv
bcftools view --include ID==@SNP_to_keep.csv AMS_v2_raw.vcf.gz  -O z -o   AMS_v2_preF_SNP.vcf.gz
```

Then we remove SNPs with a MAF < 0.01 and with more than 10% of missing data. Since the precedent work, we saw some difficulties to work with complex traits and better results with a more stringent missing data values. 

```bash
bcftools view -O u -m2 -M2 AMS_v2_preF_SNP.vcf.gz  | bcftools filter -e 'F_MISSING > 0.1 || MAF <= 0.01' -O z -o AMS_v2_F1_SNP.vcf.gz 
bcftools sort -o AMS_v2_F1_sorted_SNP.vcf.gz -O z AMS_v2_F1_SNP.vcf.gz 
```

F1 can be imputed to do GWAS or filtered for LD to work on structure. 

Imputation is always a little tricky, the only best practices actually is to use BEAGLE 5 and not BEAGLE 4. If you want to test some parameters, there is a nice paper in G3 ( https://doi.org/10.1534/g3.119.400798). As we don't have a reference panel, DR2 can't be computed, so testing could be difficult. We then filter for a MAF at 0.05. 

```bash
 beagle -Xmx80908m  gt=AMS_v2_F1_sorted_SNP.vcf.gz nthreads=20 out=AMS_v2_F1_SNP_imputed impute=true 
bcftools filter -e 'MAF <= 0.05' -O z -o AMS_v2_F1_SNP_imputed_MAF5.vcf.gz AMS_v2_F1_SNP_imputed.vcf.gz
```


For kinship, I recommand gemma :

```bash
plink --vcf AMS_v2_F1_SNP_imputed_MAF5.vcf.gz  --make-bed --out AMS_v2_F1_SNP_imputed_MAF5 --allow-extra-chr
sed -i 's/-9/1/' AMS_v2_F1_SNP_imputed_MAF5.fam 
gemma -bfile AMS_v2_F1_SNP_imputed_MAF5 -gk 1 -o AMS_v2_F1_SNP_imputed_MAF5
```
For Structure, data shouldn't be really different as what is already published in Ollivier et al. or Aubert et al. If you would like to work on it, you need to prune for LD :

```bash
plink --vcf-idspace-to "_" --double-id --allow-extra-chr --vcf AMS_v2_F1_sorted_SNP.vcf.gz --indep-pairwise 50 10 0.2 --out AMS_v2_F2_SNP
bcftools view --include ID==@AMS_v2_F2_SNP.prune.in AMS_v2_F1_sorted_SNP.vcf.gz -O z -o AMS_v2_F2_SNP.vcf.gz
```
My favorite tools to use for structure are DAPC (package adegenet), fastStructure and sometimes admixture.


