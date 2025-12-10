# VCF processing 

Tables in markdown format for before/after filtering the VCF file

# Make a common table with 3 columns (metric, before and after) in markdown to compare the number of variants before and after filtering

| Metric                      | Before | After |
|-----------------------------|--------|-------|
| Total Variants               | 100000 | 70000 |
| High Quality Variants        | 80000  | 70000 |
| Variants with Sufficient Depth| 75000  | 70000 |
| Variants with MAF > 0.01    | 70000  | 70000 |
| Variants in Coding Regions   | 20000  | 15000 |
| Transition/Transversion Ratio| 2.0    | 2.1   |
| Heterozygosity Rate          | 0.3    | 0.28  |
| Missing Data Rate            | 0.05   | 0.02  |
| Hardy-Weinberg Equilibrium   | 0.95   | 0.98  |
| Linkage Disequilibrium (r^2) | 0.8    | 0.85  |


## Comparison F1 GWAS ready file (J. Kreplak, personal communication) with ours

Comparing 
the F1 GWAS ready file (J. Kreplak, personal communication) with our filtered VCF file to see how many variants are common and how many are specific to each file.
```bash
bcftools isec  \
         --targets chr1S,chr1L,chr2,chr3,chr4,chr5,chr6 \
         -p comparisons_raw_filtered/ \
         --threads 4 \
         vcf_files/PEAMUST_vfaba_F1_SNP_imputed_MAF5.vcf.gz \
         results/filtered/vfaba_filtered_imputed_maf.vcf.gz
```

Note: F1 GWAS ready file has 247 individuals, our filtered VCF has 248 individuals.

| Metric                        | F1 GWAS Ready | Our Filtered VCF |
|-------------------------------|---------------|------------------|
| Total Variants  (only SNPs)   | 1,004,719     | 793,624          |
| Variants specific to 1 VCF    | 211,249       | 1,059            |
| Common variants               | 792,565       | 792,565          |


