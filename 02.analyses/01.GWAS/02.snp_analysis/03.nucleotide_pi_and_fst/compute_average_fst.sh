#!/bin/bash

# Define the clusters
clusters=("cluster1" "cluster2" "cluster3" "cluster4" "cluster5" "cluster6")

# Loop over each pair of clusters
for ((i = 0; i < ${#clusters[@]}; i++)); do
    for ((j = i + 1; j < ${#clusters[@]}; j++)); do
        # Define the output file name
        output_file="${clusters[i]}_vs_${clusters[j]}.weir.fst"

        # Check if the output file exists
        if [ -f "$output_file" ]; then
            # Use awk to calculate the average F_ST value
            average_fst=$(awk -F'\t' 'NR>1 && $3 != "-nan" && $3 != "nan" {sum+=$3} END {print sum/(NR-1)}' "$output_file")

            # Print the result
            echo "Average FST between ${clusters[i]} and ${clusters[j]}: $average_fst"
        else
            echo "Output file $output_file does not exist."
        fi
    done
done

