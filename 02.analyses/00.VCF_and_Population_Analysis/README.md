# VCF processing and faba bean population-related analyses 

<!-- MarkdownTOC autolink="true" levels="1,2,3" -->

- [00. Summary and objectives](#00-summary-and-objectives)
- [1. VCF processing folder](#1-vcf-processing-folder)
    - [1.1 Original VCF file description](#11-original-vcf-file-description)
    - [1.2 INFO section](#12-info-section)
    - [1.3 FORMAT section](#13-format-section)
    - [1.4 Complete header of PeaMUST Vicia faba VCF file](#14-complete-header-of-peamust-vicia-faba-vcf-file)
    - [1.5 QC metrics upon VCF processing](#15-qc-metrics-upon-vcf-processing)
- [2. Population metrics](#2-population-metrics)
    - [2.1 LD decay and distribution across the genome](#21-ld-decay-and-distribution-across-the-genome)

<!-- /MarkdownTOC -->

# 00. Summary and objectives

This folder contains (i) VCF processing and (ii) population-related analyses. 
The two main goals is to:
1. produce a VCF file that can be used in GWAS analyses (see [01.GWAS/](../01.GWAS/)) folder.    
2. extract relevant metrics and create plots related to the _Vicia faba_ PEAMUST panel (see Aubert et al. 2023 _PCI Genomics_. 

# 1. VCF processing folder

**The goal of this analysis is to:**  
1. Verify the numbers of SNPs (bi- and multiallelic variants) and InDels in the 'raw' VCF file produced from GATK best practices (J. Kreplak, pers. comm.).  
2. Filter this 'raw' VCF file using a custom [Snakemake pipeline](https://github.com/igepp-egi/snakemake_vcf) to produce a VCF file filtered for SNP quality, MAF, genotype missingness and that is imputed using Beagle v5 to recover missing SNPs.  
3. Produce QC plots before/after filtering to check that no biases have been introduced.  

## 1.1 Original VCF file description 

The 'raw' `vfaba_hedin_peamust_SNP_GATK.vcf.gz`_Vicia faba_ VCF file contains:   
- 248 genotypes.       
- 6,425,454 SNPs (bi- and multiallelic SNPs).    
- md5sum: `645c8d1fa683938b5147ea6be281db70`  

See [VCF v4.2 format](./VCFv4.2_format_specifications.pdf) for explanation of FILTER, FORMAT, etc. fields.

Meta-information lines:   

From the VCF specs:   
> File meta-information is included after the ## string and must be key=value pairs. It is strongly encouraged that
information lines describing the INFO, FILTER and FORMAT entries used in the body of the VCF file be included
in the meta-information section. Although they are optional, if these lines are present then they must be completely
well-formed.

## 1.2 INFO section 

- AC=2;  
- AF=0.333;  
- AN=6;  
- DP=7;  
- ExcessHet=0;  
- FS=0;     Phred-scaled p-value using Fisher's exact test to detect strand bias
- InbreedingCoeff=0.3366;  
- MLEAC=56;  
- MLEAF=1;  
- MQ=41;  
- QD=25.36;  
- SOR=0.693:  

| Field           | Description                                                             | Example                |
|-----------------|-------------------------------------------------------------------------|------------------------|
| AC              | Allele count in genotypes, for each ALT allele                          | AC=2                   |
| AF              | Allele Frequency, for each ALT allele                                   | AF=0.333               | 
| AN              | Total number of alleles in called genotypes                             | AN=6                   |
| DP              | Approximate read depth                                                  | DP=7                   |
| ExcessHet       | Phred-scaled p-value for exact test of excess heterozygosity            | ExcessHet=0            |
| FS              | Phred-scaled p-value using Fisher's exact test to detect strand bias    | FS=0                   |
| InbreedingCoeff | Inbreeding coefficient (1)                                              | InbreedingCoeff=0.3366 |
| MLEAC           | Maximum likelihood expectation (MLE) for the allele counts (2)          | MLEAC=56               |
| MLEAF           | Maximum likelihood expectation (MLE) for the allele frequency (2)       | MLEAF=1                |
| MQ              | RMS Mapping Quality                                                     | MQ=41                  |
| QD              | Variant Confidence/Quality by Depth                                     | QD=25.36               |
| SOR             | Symmetric Odds Ratio of 2x2 contingency table to detect strand bias     | SOR=0.693              |

(1) Inbreeding coefficient as estimated from the genotype likelihoods per-sample when compared against the Hardy-Weinberg expectation.   
(2) Maximum likelihood expectation (MLE) for the allele counts or frequency (not necessarily the same as the AC), for each ALT allele, in the same order as listed.   

## 1.3 FORMAT section

GT:AD:DP:GQ:PL  

- **GT: Genotype**. Genotype encoded as allele values separated by /: 0 is for the reference allele (REF field), 1 for the first allele listed in ALT field, 2 for the second allele, etc. Diploid values are 0/0, 0/1, 1/2 etc. '/' indicate that the genotype is unphased.    
- **AD: Allelic depths for the REF and ALT alleles in the order listed**.  
- **DP: Approximate read depth**. Approximate read depth (reads with MQ=255 or with bad mates are filtered).    
- **GQ: Genotype Quality**.   
- **PL: Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification**.   

## 1.4 Complete header of PeaMUST Vicia faba VCF file

See [VCF v4.2 format](./VCFv4.2_format_specifications.pdf) for explanation of FILTER, FORMAT, etc. fields.

```bash
fileformat=VCFv4.2
FILTER=<ID=PASS,Description="All filters passed">
ALT=<ID=NON_REF,Description="Represents any possible alternative allele at this location">
FILTER=<ID=FS60,Description="FS > 60.0">
FILTER=<ID=LowQual,Description="Low quality">
FILTER=<ID=MQ40,Description="MQ < 40.0">
FILTER=<ID=MQRankSum-12.5,Description="MQRankSum < -12.5">
FILTER=<ID=QD2,Description="QD < 2.0">
FILTER=<ID=QUAL30,Description="QUAL < 30.0">
FILTER=<ID=ReadPosRankSum-8,Description="ReadPosRankSum < -8.0">
FILTER=<ID=SOR3,Description="SOR > 3.0">

FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">
FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth (reads with MQ=255 or with bad mates are filtered)">
FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
FORMAT=<ID=MIN_DP,Number=1,Type=Integer,Description="Minimum DP observed within the GVCF block">
FORMAT=<ID=PGT,Number=1,Type=String,Description="Physical phasing haplotype information, describing how the alternate alleles are phased in relation to one another">
FORMAT=<ID=PID,Number=1,Type=String,Description="Physical phasing ID information, where each unique ID within a given sample (but not across samples) connects records within a phasing group">
FORMAT=<ID=PL,Number=G,Type=Integer,Description="Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification">
FORMAT=<ID=PS,Number=1,Type=Integer,Description="Phasing set (typically the position of the first variant in the set)">
FORMAT=<ID=RGQ,Number=1,Type=Integer,Description="Unconditional reference genotype confidence, encoded as a phred quality -10*log10 p(genotype call is wrong)">
FORMAT=<ID=SB,Number=4,Type=Integer,Description="Per-sample component statistics which comprise the Fisher's Exact Test to detect strand bias.">

GATKCommandLine=<ID=CombineGVCFs,CommandLine="CombineGVCFs --output vfaba_hedin_peamust_combined.vcf --variant UNZIP-VCF/DVF0034-id-21.output.vcf --variant UNZIP-VCF/DVF0081-id-01.output.vcf --variant UNZIP-VCF/DVF0120-id-01.output.vcf --variant UNZIP-VCF/DVF0121-id-02.output.vcf --variant UNZIP-VCF/DVF0123-id-03.output.vcf --variant UNZIP-VCF/DVF0124-id-04.output.vcf --variant UNZIP-VCF/DVF0127-id-05.output.vcf --variant UNZIP-VCF/DVF0159-id-06.output.vcf --variant UNZIP-VCF/DVF0197-id-07.output.vcf --variant UNZIP-VCF/DVF0199-id-08.output.vcf --variant UNZIP-VCF/DVF0200-id-09.output.vcf --variant UNZIP-VCF/DVF0204-id-10.output.vcf --variant UNZIP-VCF/DVF0205-id-11.output.vcf --variant UNZIP-VCF/DVF0217-id-22.output.vcf --variant UNZIP-VCF/DVF0225-id-12.output.vcf --variant UNZIP-VCF/DVF0227-id-23.output.vcf --variant UNZIP-VCF/DVF0234-id-13.output.vcf --variant UNZIP-VCF/DVF0236-id-14.output.vcf --variant UNZIP-VCF/DVF0238-id-24.output.vcf --variant UNZIP-VCF/DVF0239-id-02.output.vcf --variant UNZIP-VCF/DVF0240-id-15.output.vcf --variant UNZIP-VCF/DVF0244-id-16.output.vcf --variant UNZIP-VCF/DVF0249-id-17.output.vcf --variant UNZIP-VCF/DVF0252-id-18.output.vcf --variant UNZIP-VCF/DVF0267-id-19.output.vcf --variant UNZIP-VCF/DVF0276-id-20.output.vcf --variant UNZIP-VCF/DVF0277-id-21.output.vcf --variant UNZIP-VCF/DVF0278-id-22.output.vcf --variant UNZIP-VCF/DVF0279-id-31.output.vcf --variant UNZIP-VCF/DVF0280-id-23.output.vcf --variant UNZIP-VCF/DVF0281-id-24.output.vcf --variant UNZIP-VCF/DVF0285-id-26.output.vcf --variant UNZIP-VCF/DVF0287-id-27.output.vcf --variant UNZIP-VCF/DVF0288-id-28.output.vcf --variant UNZIP-VCF/DVF0289-id-31.output.vcf --variant UNZIP-VCF/DVF0310-id-26.output.vcf --variant UNZIP-VCF/DVF0313-id-27.output.vcf --variant UNZIP-VCF/DVF0319-id-28.output.vcf --variant UNZIP-VCF/DVF0320-id-29.output.vcf --variant UNZIP-VCF/DVF0325-id-30.output.vcf --variant UNZIP-VCF/DVF0328-id-31.output.vcf --variant UNZIP-VCF/DVF0334-id-32.output.vcf --variant UNZIP-VCF/DVF0336-id-33.output.vcf --variant UNZIP-VCF/DVF0337-id-34.output.vcf --variant UNZIP-VCF/DVF0340-id-35.output.vcf --variant UNZIP-VCF/DVF0351-id-36.output.vcf --variant UNZIP-VCF/DVF0355-id-37.output.vcf --variant UNZIP-VCF/DVF0364-id-38.output.vcf --variant UNZIP-VCF/DVF0375-id-39.output.vcf --variant UNZIP-VCF/DVF0376-id-40.output.vcf --variant UNZIP-VCF/DVF0377-id-41.output.vcf --variant UNZIP-VCF/DVF0381-id-42.output.vcf --variant UNZIP-VCF/DVF0382-id-43.output.vcf --variant UNZIP-VCF/DVF0387-id-44.output.vcf --variant UNZIP-VCF/DVF0388-id-07.output.vcf --variant UNZIP-VCF/DVF0390-id-29.output.vcf --variant UNZIP-VCF/DVF0391-id-45.output.vcf --variant UNZIP-VCF/DVF0393-id-46.output.vcf --variant UNZIP-VCF/DVF0396-id-03.output.vcf --variant UNZIP-VCF/DVF0397-id-47.output.vcf --variant UNZIP-VCF/DVF0437-id-48.output.vcf --variant UNZIP-VCF/DVF0455-id-01.output.vcf --variant UNZIP-VCF/DVF0612-id-02.output.vcf --variant UNZIP-VCF/DVF0613-id-03.output.vcf --variant UNZIP-VCF/DVF0615-id-04.output.vcf --variant UNZIP-VCF/DVF0633-id-05.output.vcf --variant UNZIP-VCF/DVF0703-id-06.output.vcf --variant UNZIP-VCF/DVF0783-id-07.output.vcf --variant UNZIP-VCF/DVF0796-id-08.output.vcf --variant UNZIP-VCF/DVF0839-id-09.output.vcf --variant UNZIP-VCF/DVF0860-id-30.output.vcf --variant UNZIP-VCF/DVF0864-id-10.output.vcf --variant UNZIP-VCF/DVF0869-id-11.output.vcf --variant UNZIP-VCF/DVF0870-id-12.output.vcf --variant UNZIP-VCF/DVF0872-id-08.output.vcf --variant UNZIP-VCF/DVF0874-id-13.output.vcf --variant UNZIP-VCF/DVF0878-id-14.output.vcf --variant UNZIP-VCF/DVF0879-id-15.output.vcf --variant UNZIP-VCF/DVF0887-id-16.output.vcf --variant UNZIP-VCF/DVF0888-id-31.output.vcf --variant UNZIP-VCF/DVF0891-id-17.output.vcf --variant UNZIP-VCF/DVF0899-id-18.output.vcf --variant UNZIP-VCF/DVF0901-id-32.output.vcf --variant UNZIP-VCF/DVF0922-id-19.output.vcf --variant UNZIP-VCF/DVF0928-id-20.output.vcf --variant UNZIP-VCF/DVF0947-id-21.output.vcf --variant UNZIP-VCF/DVF0950-id-22.output.vcf --variant UNZIP-VCF/DVF0953-id-23.output.vcf --variant UNZIP-VCF/DVF1049-id-24.output.vcf --variant UNZIP-VCF/DVF1050-id-31.output.vcf --variant UNZIP-VCF/DVF1056-id-26.output.vcf --variant UNZIP-VCF/DVF1059-id-27.output.vcf --variant UNZIP-VCF/DVF1062-id-28.output.vcf --variant UNZIP-VCF/DVF1160-id-29.output.vcf --variant UNZIP-VCF/DVF1161-id-30.output.vcf --variant UNZIP-VCF/DVF1162-id-31.output.vcf --variant UNZIP-VCF/DVF1179-id-32.output.vcf --variant UNZIP-VCF/DVF1204-id-33.output.vcf --variant UNZIP-VCF/DVF1205-id-04.output.vcf --variant UNZIP-VCF/DVF1216-id-33.output.vcf --variant UNZIP-VCF/DVF1219-id-34.output.vcf --variant UNZIP-VCF/DVF1235-id-35.output.vcf --variant UNZIP-VCF/DVF1236-id-36.output.vcf --variant UNZIP-VCF/DVF1240-id-37.output.vcf --variant UNZIP-VCF/DVF1251-id-38.output.vcf --variant UNZIP-VCF/DVF1262-id-39.output.vcf --variant UNZIP-VCF/DVF1274-id-40.output.vcf --variant UNZIP-VCF/DVF1289-id-41.output.vcf --variant UNZIP-VCF/DVF1293-id-34.output.vcf --variant UNZIP-VCF/DVF1309-id-42.output.vcf --variant UNZIP-VCF/DVF1320-id-35.output.vcf --variant UNZIP-VCF/DVF1327-id-43.output.vcf --variant UNZIP-VCF/DVF1371-id-44.output.vcf --variant UNZIP-VCF/DVF1414-id-45.output.vcf --variant UNZIP-VCF/DVF1451-id-46.output.vcf --variant UNZIP-VCF/DVF1454-id-47.output.vcf --variant UNZIP-VCF/DVF1460-id-48.output.vcf --variant UNZIP-VCF/DVF1461-id-01.output.vcf --variant UNZIP-VCF/DVF1462-id-02.output.vcf --variant UNZIP-VCF/DVF1463-id-03.output.vcf --variant UNZIP-VCF/DVF1473-id-04.output.vcf --variant UNZIP-VCF/DVF1482-id-05.output.vcf --variant UNZIP-VCF/DVF1487-id-06.output.vcf --variant UNZIP-VCF/DVF1496-id-07.output.vcf --variant UNZIP-VCF/DVF1497-id-08.output.vcf --variant UNZIP-VCF/DVF1498-id-09.output.vcf --variant UNZIP-VCF/DVF1508-id-10.output.vcf --variant UNZIP-VCF/DVF1515-id-36.output.vcf --variant UNZIP-VCF/DVF1517-id-11.output.vcf --variant UNZIP-VCF/DVF1536-id-12.output.vcf --variant UNZIP-VCF/DVF1541-id-13.output.vcf --variant UNZIP-VCF/DVF1550-id-14.output.vcf --variant UNZIP-VCF/DVF1561-id-15.output.vcf --variant UNZIP-VCF/DVF1582-id-16.output.vcf --variant UNZIP-VCF/DVF1598-id-17.output.vcf --variant UNZIP-VCF/DVF1599-id-18.output.vcf --variant UNZIP-VCF/DVF1631-id-19.output.vcf --variant UNZIP-VCF/DVF1648-id-20.output.vcf --variant UNZIP-VCF/DVF1649-id-21.output.vcf --variant UNZIP-VCF/DVF1653-id-37.output.vcf --variant UNZIP-VCF/DVF1659-id-22.output.vcf --variant UNZIP-VCF/DVF1676-id-23.output.vcf --variant UNZIP-VCF/DVF1692-id-24.output.vcf --variant UNZIP-VCF/DVF1697-id-31.output.vcf --variant UNZIP-VCF/DVF1698-id-26.output.vcf --variant UNZIP-VCF/DVF1699-id-27.output.vcf --variant UNZIP-VCF/DVF1701-id-28.output.vcf --variant UNZIP-VCF/DVF1704-id-29.output.vcf --variant UNZIP-VCF/DVF1709-id-05.output.vcf --variant UNZIP-VCF/DVF1723-id-30.output.vcf --variant UNZIP-VCF/DVF1731-id-31.output.vcf --variant UNZIP-VCF/DVF1736-id-32.output.vcf --variant UNZIP-VCF/DVF1737-id-33.output.vcf --variant UNZIP-VCF/DVF1738-id-34.output.vcf --variant UNZIP-VCF/DVF1747-id-35.output.vcf --variant UNZIP-VCF/DVF1748-id-36.output.vcf --variant UNZIP-VCF/DVF1757-id-37.output.vcf --variant UNZIP-VCF/DVF1777-id-38.output.vcf --variant UNZIP-VCF/DVF1783-id-39.output.vcf --variant UNZIP-VCF/DVF1791-id-40.output.vcf --variant UNZIP-VCF/DVF1799-id-41.output.vcf --variant UNZIP-VCF/DVF1812-id-42.output.vcf --variant UNZIP-VCF/DVF1838-id-43.output.vcf --variant UNZIP-VCF/DVF1847-id-44.output.vcf --variant UNZIP-VCF/DVF1859-id-45.output.vcf --variant UNZIP-VCF/DVF1880-id-46.output.vcf --variant UNZIP-VCF/DVF1908-id-47.output.vcf --variant UNZIP-VCF/DVF1913-id-48.output.vcf --variant UNZIP-VCF/DVF1957-id-01.output.vcf --variant UNZIP-VCF/DVF1978-id-02.output.vcf --variant UNZIP-VCF/DVF1986-id-38.output.vcf --variant UNZIP-VCF/DVF2062-id-03.output.vcf --variant UNZIP-VCF/DVF2063-id-04.output.vcf --variant UNZIP-VCF/DVF2070-id-39.output.vcf --variant UNZIP-VCF/DVF2073-id-05.output.vcf --variant UNZIP-VCF/DVF2077-id-06.output.vcf --variant UNZIP-VCF/DVF2080-id-07.output.vcf --variant UNZIP-VCF/DVF2081-id-08.output.vcf --variant UNZIP-VCF/DVF2088-id-09.output.vcf --variant UNZIP-VCF/DVF2121-id-10.output.vcf --variant UNZIP-VCF/DVF2124-id-11.output.vcf --variant UNZIP-VCF/DVF2161-id-12.output.vcf --variant UNZIP-VCF/DVF2180-id-13.output.vcf --variant UNZIP-VCF/DVF2184-id-14.output.vcf --variant UNZIP-VCF/DVF2188-id-15.output.vcf --variant UNZIP-VCF/DVF2191-id-16.output.vcf --variant UNZIP-VCF/DVF2192-id-17.output.vcf --variant UNZIP-VCF/DVF2198-id-18.output.vcf --variant UNZIP-VCF/DVF2200-id-19.output.vcf --variant UNZIP-VCF/DVF2220-id-20.output.vcf --variant UNZIP-VCF/DVF2244-id-06.output.vcf --variant UNZIP-VCF/DVF2249-id-21.output.vcf --variant UNZIP-VCF/DVF2261-id-40.output.vcf --variant UNZIP-VCF/DVF2275-id-22.output.vcf --variant UNZIP-VCF/DVF2277-id-23.output.vcf --variant UNZIP-VCF/DVF2300-id-24.output.vcf --variant UNZIP-VCF/DVF2302-id-31.output.vcf --variant UNZIP-VCF/DVF2305-id-26.output.vcf --variant UNZIP-VCF/DVF2307-id-27.output.vcf --variant UNZIP-VCF/DVF2327-id-41.output.vcf --variant UNZIP-VCF/DVF2335-id-28.output.vcf --variant UNZIP-VCF/DVF2337-id-29.output.vcf --variant UNZIP-VCF/DVF2345-id-30.output.vcf --variant UNZIP-VCF/DVF2348-id-31.output.vcf --variant UNZIP-VCF/DVF2355-id-32.output.vcf --variant UNZIP-VCF/DVF2356-id-33.output.vcf --variant UNZIP-VCF/DVF2361-id-34.output.vcf --variant UNZIP-VCF/DVF2365-id-35.output.vcf --variant UNZIP-VCF/DVF2366-id-36.output.vcf --variant UNZIP-VCF/DVF2372-id-42.output.vcf --variant UNZIP-VCF/DVF2378-id-37.output.vcf --variant UNZIP-VCF/DVF2389-id-38.output.vcf --variant UNZIP-VCF/DVF2394-id-39.output.vcf --variant UNZIP-VCF/DVF2397-id-40.output.vcf --variant UNZIP-VCF/DVF2405-id-41.output.vcf --variant UNZIP-VCF/DVF2442-id-42.output.vcf --variant UNZIP-VCF/DVF2446-id-43.output.vcf --variant UNZIP-VCF/DVF2454-id-43.output.vcf --variant UNZIP-VCF/DVF2485-id-44.output.vcf --variant UNZIP-VCF/DVF2492-id-44.output.vcf --variant UNZIP-VCF/DVF2516-id-45.output.vcf --variant UNZIP-VCF/DVF2529-id-46.output.vcf --variant UNZIP-VCF/DVF2547-id-45.output.vcf --variant UNZIP-VCF/DVF2578-id-47.output.vcf --variant UNZIP-VCF/DVF2649-id-48.output.vcf --variant UNZIP-VCF/DVF2658-id-01.output.vcf --variant UNZIP-VCF/DVF2702-id-02.output.vcf --variant UNZIP-VCF/DVF2732-id-03.output.vcf --variant UNZIP-VCF/DVF2737-id-04.output.vcf --variant UNZIP-VCF/DVF3046-id-05.output.vcf --variant UNZIP-VCF/DVF3069-id-06.output.vcf --variant UNZIP-VCF/DVF3078-id-07.output.vcf --variant UNZIP-VCF/DVF3087-id-08.output.vcf --variant UNZIP-VCF/DVF3156-id-09.output.vcf --variant UNZIP-VCF/DVF3212-id-10.output.vcf --variant UNZIP-VCF/DVF3228-id-46.output.vcf --variant UNZIP-VCF/DVF3263-id-11.output.vcf --variant UNZIP-VCF/DVF3268-id-12.output.vcf --variant UNZIP-VCF/DVF3274-id-13.output.vcf --variant UNZIP-VCF/DVF3278-id-14.output.vcf --variant UNZIP-VCF/DVF3279-id-15.output.vcf --variant UNZIP-VCF/DVF3284-id-47.output.vcf --variant UNZIP-VCF/DVF3287-id-16.output.vcf --variant UNZIP-VCF/DVF3288-id-48.output.vcf --variant UNZIP-VCF/DVF3295-id-17.output.vcf --variant UNZIP-VCF/DVF3304-id-18.output.vcf --variant UNZIP-VCF/DVF3312-id-19.output.vcf --variant UNZIP-VCF/DVF3314-id-20.output.vcf --reference Data/vfaba_hedin_v1.fasta --convert-to-base-pair-resolution false --break-bands-at-multiples-of 0 --input-is-somatic false --drop-somatic-filtering-annotations false --call-genotypes false --ignore-variants-starting-outside-interval false --combine-variants-distance 0 --max-distance 2147483647 --ref-padding 1 --interval-set-rule UNION --interval-padding 0 --interval-exclusion-padding 0 --interval-merging-rule ALL --read-validation-stringency SILENT --seconds-between-progress-updates 10.0 --disable-sequence-dictionary-validation false --create-output-bam-index true --create-output-bam-md5 false --create-output-variant-index true --create-output-variant-md5 false --max-variants-per-shard 0 --lenient false --add-output-sam-program-record true --add-output-vcf-command-line true --cloud-prefetch-buffer 40 --cloud-index-prefetch-buffer -1 --disable-bam-index-caching false --sites-only-vcf-output false --help false --version false --showHidden false --verbosity INFO --QUIET false --use-jdk-deflater false --use-jdk-inflater false --gcs-max-retries 20 --gcs-project-for-requester-pays  --disable-tool-default-read-filters false --disable-tool-default-annotations false --enable-all-annotations false --allow-old-rms-mapping-quality-annotation-data false",Version="4.2.6.1",Date="July 27, 2022 9:20:09 AM GMT">

GATKCommandLine=<ID=SelectVariants,CommandLine="SelectVariants --output vfaba_hedin_peamust_genotyped_SNP.vcf --select-type-to-include SNP --variant vfaba_hedin_peamust_genotyped.vcf --invertSelect false --exclude-non-variants false --exclude-filtered false --preserve-alleles false --remove-unused-alternates false --restrict-alleles-to ALL --keep-original-ac false --keep-original-dp false --mendelian-violation false --invert-mendelian-violation false --mendelian-violation-qual-threshold 0.0 --select-random-fraction 0.0 --remove-fraction-genotypes 0.0 --fully-decode false --max-indel-size 2147483647 --min-indel-size 0 --max-filtered-genotypes 2147483647 --min-filtered-genotypes 0 --max-fraction-filtered-genotypes 1.0 --min-fraction-filtered-genotypes 0.0 --max-nocall-number 2147483647 --max-nocall-fraction 1.0 --set-filtered-gt-to-nocall false --allow-nonoverlapping-command-line-samples false --suppress-reference-path false --genomicsdb-max-alternate-alleles 50 --call-genotypes false --genomicsdb-use-bcf-codec false --genomicsdb-shared-posixfs-optimizations false --genomicsdb-use-gcs-hdfs-connector false --interval-set-rule UNION --interval-padding 0 --interval-exclusion-padding 0 --interval-merging-rule ALL --read-validation-stringency SILENT --seconds-between-progress-updates 10.0 --disable-sequence-dictionary-validation false --create-output-bam-index true --create-output-bam-md5 false --create-output-variant-index true --create-output-variant-md5 false --max-variants-per-shard 0 --lenient false --add-output-sam-program-record true --add-output-vcf-command-line true --cloud-prefetch-buffer 40 --cloud-index-prefetch-buffer -1 --disable-bam-index-caching false --sites-only-vcf-output false --help false --version false --showHidden false --verbosity INFO --QUIET false --use-jdk-deflater false --use-jdk-inflater false --gcs-max-retries 20 --gcs-project-for-requester-pays  --disable-tool-default-read-filters false",Version="4.2.6.1",Date="August 9, 2022 8:52:23 AM GMT">

GATKCommandLine=<ID=VariantFiltration,CommandLine="VariantFiltration --output vfaba_hedin_peamust_genotyped_SNP_filtered.vcf --filter-expression QD < 2.0 --filter-expression QUAL < 30.0 --filter-expression SOR > 3.0 --filter-expression FS > 60.0 --filter-expression MQ < 40.0 --filter-expression MQRankSum < -12.5 --filter-expression ReadPosRankSum < -8.0 --filter-name QD2 --filter-name QUAL30 --filter-name SOR3 --filter-name FS60 --filter-name MQ40 --filter-name MQRankSum-12.5 --filter-name ReadPosRankSum-8 --variant vfaba_hedin_peamust_genotyped_SNP.vcf --cluster-size 3 --cluster-window-size 0 --mask-extension 0 --mask-name Mask --filter-not-in-mask false --missing-values-evaluate-as-failing false --invalidate-previous-filters false --invert-filter-expression false --invert-genotype-filter-expression false --set-filtered-genotype-to-no-call false --apply-allele-specific-filters false --interval-set-rule UNION --interval-padding 0 --interval-exclusion-padding 0 --interval-merging-rule ALL --read-validation-stringency SILENT --seconds-between-progress-updates 10.0 --disable-sequence-dictionary-validation false --create-output-bam-index true --create-output-bam-md5 false --create-output-variant-index true --create-output-variant-md5 false --max-variants-per-shard 0 --lenient false --add-output-sam-program-record true --add-output-vcf-command-line true --cloud-prefetch-buffer 40 --cloud-index-prefetch-buffer -1 --disable-bam-index-caching false --sites-only-vcf-output false --help false --version false --showHidden false --verbosity INFO --QUIET false --use-jdk-deflater false --use-jdk-inflater false --gcs-max-retries 20 --gcs-project-for-requester-pays  --disable-tool-default-read-filters false",Version="4.2.6.1",Date="9 août 2022 à 13:59:06 CEST">

INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count in genotypes, for each ALT allele, in the same order as listed">
INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency, for each ALT allele, in the same order as listed">
INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles in called genotypes">
INFO=<ID=BaseQRankSum,Number=1,Type=Float,Description="Z-score from Wilcoxon rank sum test of Alt Vs. Ref base qualities">
INFO=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth; some reads may have been filtered">
INFO=<ID=DS,Number=0,Type=Flag,Description="Were any of the samples downsampled?">
INFO=<ID=END,Number=1,Type=Integer,Description="Stop position of the interval">
INFO=<ID=ExcessHet,Number=1,Type=Float,Description="Phred-scaled p-value for exact test of excess heterozygosity">
INFO=<ID=FS,Number=1,Type=Float,Description="Phred-scaled p-value using Fisher's exact test to detect strand bias">
INFO=<ID=InbreedingCoeff,Number=1,Type=Float,Description="Inbreeding coefficient as estimated from the genotype likelihoods per-sample when compared against the Hardy-Weinberg expectation">
INFO=<ID=MLEAC,Number=A,Type=Integer,Description="Maximum likelihood expectation (MLE) for the allele counts (not necessarily the same as the AC), for each ALT allele, in the same order as listed">
INFO=<ID=MLEAF,Number=A,Type=Float,Description="Maximum likelihood expectation (MLE) for the allele frequency (not necessarily the same as the AF), for each ALT allele, in the same order as listed">
INFO=<ID=MQ,Number=1,Type=Float,Description="RMS Mapping Quality">
INFO=<ID=MQRankSum,Number=1,Type=Float,Description="Z-score From Wilcoxon rank sum test of Alt vs. Ref read mapping qualities">
INFO=<ID=QD,Number=1,Type=Float,Description="Variant Confidence/Quality by Depth">
INFO=<ID=RAW_MQandDP,Number=2,Type=Integer,Description="Raw data (sum of squared MQ and total depth) for improved RMS Mapping Quality calculation. Incompatible with deprecated RAW_MQ formulation.">
INFO=<ID=ReadPosRankSum,Number=1,Type=Float,Description="Z-score from Wilcoxon rank sum test of Alt vs. Ref read position bias">
INFO=<ID=SOR,Number=1,Type=Float,Description="Symmetric Odds Ratio of 2x2 contingency table to detect strand bias">
``` 


## 1.5 QC metrics upon VCF processing 


VCF file metrics produced by the Snakemake [snakevcf](https://github.com/igepp-egi/snakemake_vcf):   
- Total number of variants (SNPs and InDels) at each step.   
- Total number of individuals (samples).

VCF metrics to plot:  
- SNP density across the genome (before/after filtering).  
- Transitions/Transversions ratio across the genome.   

# 2. Population metrics

## 2.1 LD decay and distribution across the genome

 


