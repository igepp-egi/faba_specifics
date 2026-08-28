# Structure analysis

## Step 1: filtering raw VCF file for structure analysis

Individuals: only individuals used in GWAS analyses, phenotyped for A. pisum and/or A. fabae resistance as well as genotyped, are retained. See the [`all_ivip.tsv`](./all_ivip.tsv).  

Minor Allele Frequency was set to 1%, genotypes were imputed for SNPs with less than 5% missing values.   

All settings available in `config.yaml`. The [Snakemake pipeline snakevcf v0.4.1](https://github.com/igepp-egi/snakemake_vcf/releases/tag/v0.4.1) was used.  

## Step 2: Structure analysis using STRUCTURE or ADMIXTURE

[FIXME] Description