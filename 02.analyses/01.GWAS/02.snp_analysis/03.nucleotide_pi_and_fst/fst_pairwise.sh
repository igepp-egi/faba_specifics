#!/bin/bash 

# Define the clusters
clusters=("cluster1" "cluster2" "cluster3" "cluster4" "cluster5" "cluster6")

# Input VCF file
vcf_file="vfaba_hedin_peamust_SNP_GATK.vcf.gz"

# Loop over each pair of clusters
for ((i = 0; i < ${#clusters[@]}; i++)); do
    for ((j = i + 1; j < ${#clusters[@]}; j++)); do
        # Define the input population files
        pop1="${clusters[i]}.txt"
        pop2="${clusters[j]}.txt"

        # Define the output file name
        output_file="${clusters[i]}_vs_${clusters[j]}"

        # Run VCFtools to calculate Weir and Cockerham's F_ST
        vcftools --gzvcf "$vcf_file" \
                 --weir-fst-pop "$pop1" \
                 --weir-fst-pop "$pop2" \
                 --out "$output_file"

        echo "Comparison between ${clusters[i]} and ${clusters[j]} completed."
    done
done
