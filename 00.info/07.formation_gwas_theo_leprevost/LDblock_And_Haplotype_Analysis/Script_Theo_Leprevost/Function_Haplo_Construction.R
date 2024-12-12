########### Author: Theo Leprevost
########### Date : 09/22/2023
### Haplotype construction
#Input data : A genotyping dataframe with the lines in rownames and the markers in colnames
#Output data : A haplot_detail and a haplo_summary file with the detail and the summary of the haplotypes, respectively


Haplo_Construction<-function(data,file) {
  
  # Transform data as a data.frame
  genotyping_df<-as.data.frame(data)
  
  # Calculate haplotypes for each line
  Haplotype <- apply(genotyping_df, 1, paste, collapse = "")
  genotyping_df<-cbind(genotyping_df,data.frame(Haplotype))
  
  # Count the number of lines sharing each haplotype
  haplotype_counts <- table(Haplotype)
  
  # Calculate the percentage of lines sharing each haplotype
  haplotype_percentages <- prop.table(haplotype_counts) * 100
  
  # Create a dataframe with haplotype information
  haplotype_df <- data.frame(
    Percentage_of_Lines = haplotype_percentages,
    Number_of_Lines = as.integer(haplotype_counts),
    stringsAsFactors = FALSE
  )
  
  # Add column names
  colnames(haplotype_df) <- c("Haplotype", "Percentage_of_Lines","Number_of_Lines")
  
  # Store the rownames of genotyping_df in a separate column
  genotyping_df$Line <- rownames(genotyping_df)
  
  # Merge the haplotype information back into the original genotyping dataframe
  merged_genotyping_df <- merge(genotyping_df, haplotype_df, by.x = "Haplotype", by.y = "Haplotype", all.x = TRUE)
  
  # Restore the rownames to the merged dataframe
  rownames(merged_genotyping_df) <- merged_genotyping_df$Line
  merged_genotyping_df$Line <- NULL
  
  # Function to check for heterozygous sites
  has_heterozygous <- function(row) {
    any(row != "AA" & row != "TT" & row != "CC" & row != "GG")
  }
  
  # Add Heterozygous_Haplotype column to merged_genotyping_df indicating heterozygous haplotypes
  merged_genotyping_df$Heterozygous_Haplotype <- apply(merged_genotyping_df[, -c(1,ncol(merged_genotyping_df)-1,ncol(merged_genotyping_df))], 1, has_heterozygous)
  
  # Sort the haplotype dataframe by Heterozygous_Haplotype and then by the number of lines (descending order)
  merged_genotyping_df <- merged_genotyping_df[order(merged_genotyping_df$Heterozygous_Haplotype), ]
  merged_genotyping_df <- merged_genotyping_df[order(-merged_genotyping_df$Number_of_Lines), ]
  
  
  # Sort the merged_genotyping_df dataframe by unique combination of haplotypes
  merged_haplotype_df <- merged_genotyping_df[!duplicated(merged_genotyping_df$Haplotype), ]
  merged_haplotype_df <- merged_haplotype_df[order(-merged_haplotype_df$Number_of_Lines), ]
  rownames(merged_haplotype_df) <- 1:nrow(merged_haplotype_df)
  
  # Assign lowercase letters to haplotypes based on their rank
  LETTERS702 <- c(LETTERS, sapply(LETTERS, function(x) paste0(x, LETTERS)))
  Haplotype_ID<-LETTERS702[seq_len(nrow(merged_haplotype_df))]
  merged_haplotype_df <- cbind(merged_haplotype_df,Haplotype_ID)
  
  # Merge the haplotype ID information  into the merged_genotyping_df genotyping dataframe
  merged_genotyping_df$Haplotype_ID <- merged_haplotype_df$Haplotype_ID[match(merged_genotyping_df$Haplotype, merged_haplotype_df$Haplotype)]  
  
  # Save the merged dataframe with haplotype information
  Code<-row.names(merged_genotyping_df)
  merged_genotyping_df<-cbind(Code,merged_genotyping_df)
  rownames(merged_genotyping_df) <- 1:nrow(merged_genotyping_df)
  write.table(merged_genotyping_df,file=paste(file,"_haplo_detail.csv",sep=""),sep=";",row.names = FALSE)
  
  # Save the summary dataframe
  write.table(merged_haplotype_df,file=paste(file,"_haplo_summary.csv",sep=""),sep=";",row.names = FALSE)
}
  
  
  
  
  
  
  
  
  