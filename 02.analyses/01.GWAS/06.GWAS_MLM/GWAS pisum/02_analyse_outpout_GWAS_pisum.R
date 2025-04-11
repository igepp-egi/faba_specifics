#########################

# Load necessary libraries
library(dplyr)
library(tidyr)

df <- read.csv("../../../Documents/QTL_south_fabae.csv", sep = ",")



# Create the new dataset
df_new <- df %>%
  mutate(SNPb = paste(SNP,GWASmethod, sep = "-")) %>%
  # Remove duplicates to avoid counting the same SNP multiple times for the same combination of filter and repetition
  distinct(filter, nrep, SNPb) %>%
  # Use pivot_wider to create a column for each unique SNP
  pivot_wider(names_from = SNPb, values_from = SNPb,
              values_fn = ~1, 
              values_fill = list(SNP = 0))

# View the new dataset
print(df_new)



# Separate first two columns and the rest
first_two_columns <- df_new[, 1:2]
remaining_columns <- df_new[, 3:ncol(df_new)]


df<-remaining_columns


# Identify all columns in the dataset
cols <- colnames(df)

# Create missing columns if one of the pair (Farm.cpu / Blink) is missing
for (col in cols) {
  if (grepl("Farm.cpu", col)) {
    # Create the corresponding Blink column
    blink_col <- gsub("Farm.cpu", "Blink", col)
    
    # Check if the Blink column exists
    if (!(blink_col %in% cols)) {
      # If it does not exist, create it (but we won't assign any values, just create the column)
      df[[blink_col]] <- NA
    }
  } else if (grepl("Blink", col)) {
    # Create the corresponding Farm.cpu column
    farm_col <- gsub("Blink", "Farm.cpu", col)
    
    # Check if the Farm.cpu column exists
    if (!(farm_col %in% cols)) {
      # If it does not exist, create it (but we won't assign any values, just create the column)
      df[[farm_col]] <- NA
    }
  }
}


remaining_columns <- df





# Sort the remaining columns alphabetically
sorted_remaining_columns <- remaining_columns[, order(names(remaining_columns))]

# Combine the first two columns with the sorted ones
sorted_df <- cbind(first_two_columns, sorted_remaining_columns)





file_path <- "~/Documents/gitrepo/faba_specifics/02.analyses/01.GWAS/06.GWAS_MLM/outpoutGWAS_Fabae_south.csv"
# Save the data frame to the specified location
write.csv(sorted_df, file = file_path, row.names = FALSE)






