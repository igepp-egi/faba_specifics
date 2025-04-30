# SNP analysis 

## Analysis 1: check numbers

Goal: to verify numbers after filtering going from the raw VCF file to the GWAS ready file.   

## Analysis 2: descriptive stats

Total Number of Variants:  
    The total count of variants (SNPs and indels) identified in the dataset.  
 
Variant Types:  
    The distribution of different types of variants, such as SNPs, insertions, deletions, and complex variants.  

Variant Quality Scores:  
    The distribution of quality scores (e.g., QUAL field in VCF) for the variants, which indicates the confidence in the variant calls.  

Allele Frequency Spectrum:  
    The distribution of allele frequencies across the samples, which can help identify common and rare variants.  

Depth of Coverage:  
    The average read depth at variant sites, which affects the confidence in variant calls.  

Genotype Quality:
    The distribution of genotype quality scores (e.g., GQ field in VCF) for the samples, indicating the confidence in the genotype calls.

Heterozygosity Rate:
    The proportion of heterozygous genotypes, which can vary depending on the population and the type of sample (e.g., inbred lines vs. outbred populations).

Mendelian Errors:
    The number of Mendelian inconsistencies in family-based studies, which can indicate errors in variant calling or sample relationships.

Hardy-Weinberg Equilibrium (HWE):
    Deviations from Hardy-Weinberg equilibrium, which can indicate genotyping errors or population stratification.

Missingness Rate:
    The proportion of missing genotype calls, which can affect downstream analyses.

Functional Annotation:
    The distribution of variants across different functional categories (e.g., exonic, intronic, intergenic) and their predicted effects (e.g., synonymous, nonsynonymous, stop-gain).
 
## Analysis 3: calculate LD decay 

[FIXME]

## Analysis 4: genetic diversity between and per genetic cluster

Intra-cluster nucleotide diversity.   
Comparison of Fst between clusters.   