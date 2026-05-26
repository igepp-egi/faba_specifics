###packages----
library(dplyr)
library(pheatmap)
library(ggplot2)
library(ggrepel)
library(readr)
library(readxl)
#----
setwd("~/Downloads/Paul Fremaux/untargeted/Negative mode")
##batch_3_negative----

batch3_negative <- read_excel("final_filtered_files/batch_3_Intensity_filtered_notame.xlsx", col_names = FALSE)

batch3_negative <- batch3_negative %>%
  mutate(across(34, ~ ifelse(!is.na(batch3_negative[[3]]) & batch3_negative[[3]] != "", batch3_negative[[3]], .)))
batch3_negative <- batch3_negative[,-c(1:33)]
batch3_negative <- t(batch3_negative)
batch3_negative <- as.data.frame(batch3_negative)
colnames(batch3_negative) <- batch3_negative[1, ]
batch3_negative <- batch3_negative[-1, ]
batch3_negative_row_name <- batch3_negative %>%
  filter(Group != "Blank") %>%
  separate(Group, into = c("Batch","treatmentA","hpi"), sep = "-") %>%
  mutate(treatment = paste(Batch, treatmentA, sep = "-")) %>%
  select(FEATURE_ID,treatment, hpi)
batch3_negative <- merge(batch3_negative_row_name,batch3_negative, by ="FEATURE_ID")
batch3_negative[is.na(batch3_negative)] <- 0 


#batch3_negative_LSR1----
# Keep rows where treatment is "E3-LSR1" OR another value (e.g., "something" or "something else")
batch3_negative_LSR1 <- batch3_negative[
  batch3_negative$treatment %in% c("E3-LSR1", "E3-control"),
]

# Remove unwanted columns
batch3_negative_LSR1 <- batch3_negative_LSR1[
  , !(names(batch3_negative_LSR1) %in% c("FEATURE_ID", "hpi", "Sample_ID", "Injection_order", "Group", "QC"))
]

# Transpose the dataframe
batch3_negative_LSR1 <- t(batch3_negative_LSR1)


# 2. Extract treatment and data

treatment <- as.character(unlist(batch3_negative_LSR1[1, ])) 
data <- batch3_negative_LSR1[-1, ]                         
data <- apply(data, 2, as.numeric)            # convert to numeric
rownames(data) <- rownames(batch3_negative_LSR1)[-1]

# 3. Compute log2 fold change

mean_R <- rowMeans(data[, treatment == "E3-LSR1"])
mean_S <- rowMeans(data[, treatment == "E3-control"])

log2FC <- log2(mean_R / mean_S)   


# 4. Compute p-values (t-test)

pvalues <- apply(data, 1, function(x) {
  t.test(x[treatment == "E3-LSR1"], x[treatment == "E3-control"])$p.value
})


# 5. Adjust p-values (FDR)

padj <- p.adjust(pvalues, method = "fdr")


# 6. Create results_batch3_negative_LSR1 table

results_batch3_negative_LSR1 <- data.frame(
  metabolite = rownames(data),
  log2FC = log2FC,
  pvalue = pvalues,
  padj = padj
)

results_batch3_negative_LSR1$negLog10P <- -log10(results_batch3_negative_LSR1$padj)


# 7. Volcano plot

# Define significance
results_batch3_negative_LSR1$significant <- "Not Significant"
results_batch3_negative_LSR1$significant[results_batch3_negative_LSR1$log2FC > 1 & results_batch3_negative_LSR1$padj < 0.05] <- "Up"
results_batch3_negative_LSR1$significant[results_batch3_negative_LSR1$log2FC < -1 & results_batch3_negative_LSR1$padj < 0.05] <- "Down"

Volcano_batch3_negative_LSR1 <- ggplot(results_batch3_negative_LSR1, aes(x = log2FC, y = negLog10P, color = significant)) +
  geom_point(size = 3, alpha = 0.8) +
  scale_color_manual(values = c("blue", "grey", "red")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  theme_minimal() +
  labs(
    x = "Log2 Fold Change (E1-AF / E1-C)",
    y = "-Log10(FDR)",
    title = "Volcano Plot"
  ) +
  geom_text_repel(
    data = subset(results_batch3_negative_LSR1, significant != "Not Significant"),
    aes(label = metabolite)
  )
results_batch3_negative_LSR1$significant<-as.factor(results_batch3_negative_LSR1$significant)
summary(results_batch3_negative_LSR1)



#batch3_negative_P123----
# Keep rows where treatment is "E3-P123" OR another value (e.g., "something" or "something else")
batch3_negative_P123 <- batch3_negative[
  batch3_negative$treatment %in% c("E3-P123", "E3-control"),
]

# Remove unwanted columns
batch3_negative_P123 <- batch3_negative_P123[
  , !(names(batch3_negative_P123) %in% c("FEATURE_ID", "hpi", "Sample_ID", "Injection_order", "Group", "QC"))
]

# Transpose the dataframe
batch3_negative_P123 <- t(batch3_negative_P123)


# 2. Extract treatment and data

treatment <- as.character(unlist(batch3_negative_P123[1, ])) 
data <- batch3_negative_P123[-1, ]                         
data <- apply(data, 2, as.numeric)            # convert to numeric
rownames(data) <- rownames(batch3_negative_P123)[-1]

# 3. Compute log2 fold change

mean_R <- rowMeans(data[, treatment == "E3-P123"])
mean_S <- rowMeans(data[, treatment == "E3-control"])

log2FC <- log2(mean_R / mean_S)   


# 4. Compute p-values (t-test)

pvalues <- apply(data, 1, function(x) {
  t.test(x[treatment == "E3-P123"], x[treatment == "E3-control"])$p.value
})


# 5. Adjust p-values (FDR)

padj <- p.adjust(pvalues, method = "fdr")


# 6. Create results_batch3_negative_P123 table

results_batch3_negative_P123 <- data.frame(
  metabolite = rownames(data),
  log2FC = log2FC,
  pvalue = pvalues,
  padj = padj
)

results_batch3_negative_P123$negLog10P <- -log10(results_batch3_negative_P123$padj)


# 7. Volcano plot

# Define significance
results_batch3_negative_P123$significant <- "Not Significant"
results_batch3_negative_P123$significant[results_batch3_negative_P123$log2FC > 1 & results_batch3_negative_P123$padj < 0.05] <- "Up"
results_batch3_negative_P123$significant[results_batch3_negative_P123$log2FC < -1 & results_batch3_negative_P123$padj < 0.05] <- "Down"

Volcano_batch3_negative_P123 <- ggplot(results_batch3_negative_P123, aes(x = log2FC, y = negLog10P, color = significant)) +
  geom_point(size = 3, alpha = 0.8) +
  scale_color_manual(values = c("blue", "grey", "red")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  theme_minimal() +
  labs(
    x = "Log2 Fold Change (E1-AF / E1-C)",
    y = "-Log10(FDR)",
    title = "Volcano Plot"
  ) +
  geom_text_repel(
    data = subset(results_batch3_negative_P123, significant != "Not Significant"),
    aes(label = metabolite)
  )
results_batch3_negative_P123$significant<-as.factor(results_batch3_negative_P123$significant)
summary(results_batch3_negative_P123)


#batch3_negative_YR2----
# Keep rows where treatment is "E3-YR2" OR another value (e.g., "something" or "something else")
batch3_negative_YR2 <- batch3_negative[
  batch3_negative$treatment %in% c("E3-YR2", "E3-control"),
]

# Remove unwanted columns
batch3_negative_YR2 <- batch3_negative_YR2[
  , !(names(batch3_negative_YR2) %in% c("FEATURE_ID", "hpi", "Sample_ID", "Injection_order", "Group", "QC"))
]

# Transpose the dataframe
batch3_negative_YR2 <- t(batch3_negative_YR2)


# 2. Extract treatment and data

treatment <- as.character(unlist(batch3_negative_YR2[1, ])) 
data <- batch3_negative_YR2[-1, ]                         
data <- apply(data, 2, as.numeric)            # convert to numeric
rownames(data) <- rownames(batch3_negative_YR2)[-1]

# 3. Compute log2 fold change

mean_R <- rowMeans(data[, treatment == "E3-YR2"])
mean_S <- rowMeans(data[, treatment == "E3-control"])

log2FC <- log2(mean_R / mean_S)   


# 4. Compute p-values (t-test)

pvalues <- apply(data, 1, function(x) {
  t.test(x[treatment == "E3-YR2"], x[treatment == "E3-control"])$p.value
})


# 5. Adjust p-values (FDR)

padj <- p.adjust(pvalues, method = "fdr")


# 6. Create results_batch3_negative_YR2 table

results_batch3_negative_YR2 <- data.frame(
  metabolite = rownames(data),
  log2FC = log2FC,
  pvalue = pvalues,
  padj = padj
)

results_batch3_negative_YR2$negLog10P <- -log10(results_batch3_negative_YR2$padj)


# 7. Volcano plot

# Define significance
results_batch3_negative_YR2$significant <- "Not Significant"
results_batch3_negative_YR2$significant[results_batch3_negative_YR2$log2FC > 1 & results_batch3_negative_YR2$padj < 0.05] <- "Up"
results_batch3_negative_YR2$significant[results_batch3_negative_YR2$log2FC < -1 & results_batch3_negative_YR2$padj < 0.05] <- "Down"

Volcano_batch3_negative_YR2 <- ggplot(results_batch3_negative_YR2, aes(x = log2FC, y = negLog10P, color = significant)) +
  geom_point(size = 3, alpha = 0.8) +
  scale_color_manual(values = c("blue", "grey", "red")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  theme_minimal() +
  labs(
    x = "Log2 Fold Change (E1-AF / E1-C)",
    y = "-Log10(FDR)",
    title = "Volcano Plot"
  ) +
  geom_text_repel(
    data = subset(results_batch3_negative_YR2, significant != "Not Significant"),
    aes(label = metabolite)
  )
results_batch3_negative_YR2$significant<-as.factor(results_batch3_negative_YR2$significant)
summary(results_batch3_negative_YR2)




#----
setwd("~/Downloads/Paul Fremaux/untargeted/Positive mode")
##batch_3_positive----

batch3_positive <- read_excel("final_filtered_files/batch_3_positive_Intensity_filtered_notame.xlsx", col_names = FALSE)

batch3_positive <- batch3_positive %>%
  mutate(across(36, ~ ifelse(!is.na(batch3_positive[[3]]) & batch3_positive[[3]] != "", batch3_positive[[3]], .)))
batch3_positive <- batch3_positive[,-c(1:35)]
batch3_positive <- t(batch3_positive)
batch3_positive <- as.data.frame(batch3_positive)
colnames(batch3_positive) <- batch3_positive[1, ]
batch3_positive <- batch3_positive[-1, ]
batch3_positive_row_name <- batch3_positive %>%
  filter(Group != "Blank") %>%
  separate(Group, into = c("Batch","treatmentA","hpi"), sep = "-") %>%
  mutate(treatment = paste(Batch, treatmentA, sep = "-")) %>%
  select(FEATURE_ID,treatment, hpi)
batch3_positive <- merge(batch3_positive_row_name,batch3_positive, by ="FEATURE_ID")
batch3_positive[is.na(batch3_positive)] <- 0 


#batch3_positive_LSR1----
# Keep rows where treatment is "E3-LSR1" OR another value (e.g., "something" or "something else")
batch3_positive_LSR1 <- batch3_positive[
  batch3_positive$treatment %in% c("E3-LSR1", "E3-control"),
]

# Remove unwanted columns
batch3_positive_LSR1 <- batch3_positive_LSR1[
  , !(names(batch3_positive_LSR1) %in% c("FEATURE_ID", "hpi", "Sample_ID", "Injection_order", "Group", "QC"))
]

# Transpose the dataframe
batch3_positive_LSR1 <- t(batch3_positive_LSR1)


# 2. Extract treatment and data

treatment <- as.character(unlist(batch3_positive_LSR1[1, ])) 
data <- batch3_positive_LSR1[-1, ]                         
data <- apply(data, 2, as.numeric)            # convert to numeric
rownames(data) <- rownames(batch3_positive_LSR1)[-1]

# 3. Compute log2 fold change

mean_R <- rowMeans(data[, treatment == "E3-LSR1"])
mean_S <- rowMeans(data[, treatment == "E3-control"])

log2FC <- log2(mean_R / mean_S)   


# 4. Compute p-values (t-test)

pvalues <- apply(data, 1, function(x) {
  t.test(x[treatment == "E3-LSR1"], x[treatment == "E3-control"])$p.value
})


# 5. Adjust p-values (FDR)

padj <- p.adjust(pvalues, method = "fdr")


# 6. Create results_batch3_positive_LSR1 table

results_batch3_positive_LSR1 <- data.frame(
  metabolite = rownames(data),
  log2FC = log2FC,
  pvalue = pvalues,
  padj = padj
)

results_batch3_positive_LSR1$negLog10P <- -log10(results_batch3_positive_LSR1$padj)


# 7. Volcano plot

# Define significance
results_batch3_positive_LSR1$significant <- "Not Significant"
results_batch3_positive_LSR1$significant[results_batch3_positive_LSR1$log2FC > 1 & results_batch3_positive_LSR1$padj < 0.05] <- "Up"
results_batch3_positive_LSR1$significant[results_batch3_positive_LSR1$log2FC < -1 & results_batch3_positive_LSR1$padj < 0.05] <- "Down"

Volcano_batch3_positive_LSR1 <- ggplot(results_batch3_positive_LSR1, aes(x = log2FC, y = negLog10P, color = significant)) +
  geom_point(size = 3, alpha = 0.8) +
  scale_color_manual(values = c("grey", "red")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  theme_minimal() +
  labs(
    x = "Log2 Fold Change (E1-AF / E1-C)",
    y = "-Log10(FDR)",
    title = "Volcano Plot"
  ) +
  geom_text_repel(
    data = subset(results_batch3_positive_LSR1, significant != "Not Significant"),
    aes(label = metabolite)
  )
results_batch3_positive_LSR1$significant<-as.factor(results_batch3_positive_LSR1$significant)
summary(results_batch3_positive_LSR1)



#batch3_positive_P123----
# Keep rows where treatment is "E3-P123" OR another value (e.g., "something" or "something else")
batch3_positive_P123 <- batch3_positive[
  batch3_positive$treatment %in% c("E3-P123", "E3-control"),
]

# Remove unwanted columns
batch3_positive_P123 <- batch3_positive_P123[
  , !(names(batch3_positive_P123) %in% c("FEATURE_ID", "hpi", "Sample_ID", "Injection_order", "Group", "QC"))
]

# Transpose the dataframe
batch3_positive_P123 <- t(batch3_positive_P123)


# 2. Extract treatment and data

treatment <- as.character(unlist(batch3_positive_P123[1, ])) 
data <- batch3_positive_P123[-1, ]                         
data <- apply(data, 2, as.numeric)            # convert to numeric
rownames(data) <- rownames(batch3_positive_P123)[-1]

# 3. Compute log2 fold change

mean_R <- rowMeans(data[, treatment == "E3-P123"])
mean_S <- rowMeans(data[, treatment == "E3-control"])

log2FC <- log2(mean_R / mean_S)   


# 4. Compute p-values (t-test)

pvalues <- apply(data, 1, function(x) {
  t.test(x[treatment == "E3-P123"], x[treatment == "E3-control"])$p.value
})


# 5. Adjust p-values (FDR)

padj <- p.adjust(pvalues, method = "fdr")


# 6. Create results_batch3_positive_P123 table

results_batch3_positive_P123 <- data.frame(
  metabolite = rownames(data),
  log2FC = log2FC,
  pvalue = pvalues,
  padj = padj
)

results_batch3_positive_P123$negLog10P <- -log10(results_batch3_positive_P123$padj)


# 7. Volcano plot

# Define significance
results_batch3_positive_P123$significant <- "Not Significant"
results_batch3_positive_P123$significant[results_batch3_positive_P123$log2FC > 1 & results_batch3_positive_P123$padj < 0.05] <- "Up"
results_batch3_positive_P123$significant[results_batch3_positive_P123$log2FC < -1 & results_batch3_positive_P123$padj < 0.05] <- "Down"

Volcano_batch3_positive_P123 <- ggplot(results_batch3_positive_P123, aes(x = log2FC, y = negLog10P, color = significant)) +
  geom_point(size = 3, alpha = 0.8) +
  scale_color_manual(values = c("grey", "red")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  theme_minimal() +
  labs(
    x = "Log2 Fold Change (E1-AF / E1-C)",
    y = "-Log10(FDR)",
    title = "Volcano Plot"
  ) +
  geom_text_repel(
    data = subset(results_batch3_positive_P123, significant != "Not Significant"),
    aes(label = metabolite)
  )
results_batch3_positive_P123$significant<-as.factor(results_batch3_positive_P123$significant)
summary(results_batch3_positive_P123)


#batch3_positive_YR2----
# Keep rows where treatment is "E3-YR2" OR another value (e.g., "something" or "something else")
batch3_positive_YR2 <- batch3_positive[
  batch3_positive$treatment %in% c("E3-YR2", "E3-control"),
]

# Remove unwanted columns
batch3_positive_YR2 <- batch3_positive_YR2[
  , !(names(batch3_positive_YR2) %in% c("FEATURE_ID", "hpi", "Sample_ID", "Injection_order", "Group", "QC"))
]

# Transpose the dataframe
batch3_positive_YR2 <- t(batch3_positive_YR2)


# 2. Extract treatment and data

treatment <- as.character(unlist(batch3_positive_YR2[1, ])) 
data <- batch3_positive_YR2[-1, ]                         
data <- apply(data, 2, as.numeric)            # convert to numeric
rownames(data) <- rownames(batch3_positive_YR2)[-1]

# 3. Compute log2 fold change

mean_R <- rowMeans(data[, treatment == "E3-YR2"])
mean_S <- rowMeans(data[, treatment == "E3-control"])

log2FC <- log2(mean_R / mean_S)   


# 4. Compute p-values (t-test)

pvalues <- apply(data, 1, function(x) {
  t.test(x[treatment == "E3-YR2"], x[treatment == "E3-control"])$p.value
})


# 5. Adjust p-values (FDR)

padj <- p.adjust(pvalues, method = "fdr")


# 6. Create results_batch3_positive_YR2 table

results_batch3_positive_YR2 <- data.frame(
  metabolite = rownames(data),
  log2FC = log2FC,
  pvalue = pvalues,
  padj = padj
)

results_batch3_positive_YR2$negLog10P <- -log10(results_batch3_positive_YR2$padj)


# 7. Volcano plot

# Define significance
results_batch3_positive_YR2$significant <- "Not Significant"
results_batch3_positive_YR2$significant[results_batch3_positive_YR2$log2FC > 1 & results_batch3_positive_YR2$padj < 0.05] <- "Up"
results_batch3_positive_YR2$significant[results_batch3_positive_YR2$log2FC < -1 & results_batch3_positive_YR2$padj < 0.05] <- "Down"

Volcano_batch3_positive_YR2 <- ggplot(results_batch3_positive_YR2, aes(x = log2FC, y = negLog10P, color = significant)) +
  geom_point(size = 3, alpha = 0.8) +
  scale_color_manual(values = c("blue", "grey", "red")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  theme_minimal() +
  labs(
    x = "Log2 Fold Change (E1-AF / E1-C)",
    y = "-Log10(FDR)",
    title = "Volcano Plot"
  ) +
  geom_text_repel(
    data = subset(results_batch3_positive_YR2, significant != "Not Significant"),
    aes(label = metabolite)
  )
results_batch3_positive_YR2$significant<-as.factor(results_batch3_positive_YR2$significant)
summary(results_batch3_positive_YR2)





ID_positive <- read_tsv('metabolites_identification/SIRIUS_Positive_mode/canopus_formula_summary.tsv')
colnames(ID_positive)[colnames(ID_positive) == "mappingFeatureId"] <- "metabolite"

ID_negative <- read_tsv('../Negative mode/metabolites identification/SIRIUS_Negative_mode/canopus_formula_summary.tsv')
colnames(ID_negative)[colnames(ID_negative) == "mappingFeatureId"] <- "metabolite"


results_batch3_negative_LSR1_sub <- results_batch3_negative_LSR1[results_batch3_negative_LSR1$significant != "Not Significant", ]
colnames(results_batch3_negative_LSR1_sub)[colnames(results_batch3_negative_LSR1_sub) == "significant"] <- "significant_LSR1"

results_batch3_negative_YR2_sub <- results_batch3_negative_YR2[results_batch3_negative_YR2$significant != "Not Significant", ]
colnames(results_batch3_negative_YR2_sub)[colnames(results_batch3_negative_YR2_sub) == "significant"] <- "significant_YR2"

results_batch3_negative_P123_sub <- results_batch3_negative_P123[results_batch3_negative_P123$significant != "Not Significant", ]
colnames(results_batch3_negative_P123_sub)[colnames(results_batch3_negative_P123_sub) == "significant"] <- "significant_P123"



results_batch3_negative_YR2_LSR1_sub <- merge(results_batch3_negative_LSR1_sub,results_batch3_negative_YR2_sub, by ="metabolite",all=TRUE)
results_batch3_negative_YR2_LSR1_P123_sub <- merge(results_batch3_negative_YR2_LSR1_sub,results_batch3_negative_P123_sub, by ="metabolite",all=TRUE)


results_batch3_negative_YR2_LSR1_P123_sub_ID <- merge(results_batch3_negative_YR2_LSR1_P123_sub, ID_negative, by ="metabolite",all=FALSE)

results_batch3_positive_LSR1_sub <- results_batch3_positive_LSR1[results_batch3_positive_LSR1$significant != "Not Significant", ]
colnames(results_batch3_positive_LSR1_sub)[colnames(results_batch3_positive_LSR1_sub) == "significant"] <- "significant_LSR1"

results_batch3_positive_YR2_sub <- results_batch3_positive_YR2[results_batch3_positive_YR2$significant != "Not Significant", ]
colnames(results_batch3_positive_YR2_sub)[colnames(results_batch3_positive_YR2_sub) == "significant"] <- "significant_YR2"

results_batch3_positive_P123_sub <- results_batch3_positive_P123[results_batch3_positive_P123$significant != "Not Significant", ]
colnames(results_batch3_positive_P123_sub)[colnames(results_batch3_positive_P123_sub) == "significant"] <- "significant_P123"



results_batch3_positive_YR2_LSR1_sub <- merge(results_batch3_positive_LSR1_sub,results_batch3_positive_YR2_sub, by ="metabolite",all=TRUE)
results_batch3_positive_YR2_LSR1_P123_sub <- merge(results_batch3_positive_YR2_LSR1_sub,results_batch3_positive_P123_sub, by ="metabolite",all=TRUE)


results_batch3_positive_YR2_LSR1_P123_sub_ID <- merge(results_batch3_positive_YR2_LSR1_P123_sub, ID_positive, by ="metabolite",all=FALSE)






# Create a list of metabolite names for each group
up_P123 <- results_batch3_positive_YR2_LSR1_P123_sub_ID %>%
  filter(significant_P123 == "Up") %>%
  pull(metabolite)


up_YR2 <- results_batch3_positive_YR2_LSR1_P123_sub_ID %>%
  filter(significant_YR2 == "Up") %>%
  pull(metabolite)

up_LSR1 <- results_batch3_positive_YR2_LSR1_P123_sub_ID %>%
  filter(significant_LSR1 == "Up") %>%
  pull(metabolite)


# Name the lists
list_sets <- list(
  "Up YR2" = up_YR2,
  "up LSR1" = up_LSR1,
  "Up P123" = up_P123
)


batch3_positive_venn<-ggvenn(list_sets, fill_color = c("#1b9e77", "#d95f02", "#7570b3", "#e7298a"),
                             stroke_size = 0.5, set_name_size = 4.5,show_percentage = FALSE)




# Create a list of metabolite names for each group
up_P123 <- results_batch3_negative_YR2_LSR1_P123_sub_ID %>%
  filter(significant_P123 == "Up") %>%
  pull(metabolite)


up_YR2 <- results_batch3_negative_YR2_LSR1_P123_sub_ID %>%
  filter(significant_YR2 == "Up") %>%
  pull(metabolite)

up_LSR1 <- results_batch3_negative_YR2_LSR1_P123_sub_ID %>%
  filter(significant_LSR1 == "Up") %>%
  pull(metabolite)


# Name the lists
list_sets <- list(
  "Up YR2" = up_YR2,
  "up LSR1" = up_LSR1,
  "Up P123" = up_P123
)


batch3_negative_venn<-ggvenn(list_sets, fill_color = c("#1b9e77", "#d95f02", "#7570b3", "#e7298a"),
                             stroke_size = 0.5, set_name_size = 4.5,show_percentage = FALSE)

