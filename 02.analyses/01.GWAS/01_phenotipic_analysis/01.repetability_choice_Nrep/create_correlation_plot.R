require(dplyr)
require(tidyr)
require(corrplot)


# Define the function to process data and generate correlation plots
create_correlation_plot <- function(data, 
                                    group_cols, 
                                    value_col, 
                                    save_plot = FALSE, 
                                    plot_filename = "correlation_plot.png") {
  
  # Step 1: Summarize duplicates by calculating the mean of the specified value column
  data_summarized <- data %>%
    group_by(across(all_of(group_cols))) %>%  # Group by the relevant columns (e.g., genotype, session)
    summarise(
      mean_value = mean(!!sym(value_col), na.rm = TRUE), # Calculate the mean for the specified value column
      .groups = "drop" # Remove grouping after summarizing
    )
  
  # Step 2: Pivot the summarized dataset to widen the data (one column per unique value in 'session' or other categorical column)
  data_wide <- data_summarized %>%
    pivot_wider(
      names_from = group_cols[2],              # Use the second column in 'group_cols' to create new columns (e.g., session)
      values_from = mean_value                  # Fill new columns with the calculated mean values
    ) %>%
    select(-group_cols[1])  # Remove the first grouping column (e.g., genotype)
  
  # Step 3: Compute the correlation matrix between the session columns
  
  # Compute correlation coefficients and p-values using rcorr
  rest <- rcorr(as.matrix(data_wide), type = "pearson")
  
  # Create the correlation matrix for plotting
  M2 <- rest$r
  diag(M2) <- NA  # Set the diagonal of the correlation matrix to NA
  
  # Step 4: Generate the correlation plot
  corr_plot <- corrplot(M2,
                        title = 'Correlation between Sessions',
                        na.label = "NA",
                        type = "upper",
                        method = "square",
                        p.mat = rest$P,
                        sig.level = 0.05,
                        tl.cex = 1.2,
                        tl.srt = 90,
                        tl.col = "black",
                        addCoef.col = "black",
                        mar = c(10, 0, 1, 0))
  
  # If save_plot is TRUE, save the plot as an image
  if (save_plot) {
    # Save the plot as a PNG file (or any other format you prefer)
    ggsave(plot_filename, plot = corr_plot, width = 8, height = 6)
    cat("Plot saved as", plot_filename, "\n")
  }
  
  # Return the correlation plot
  return(corr_plot)
}