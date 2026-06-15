
###packages----
library(dplyr)
library(pheatmap)
library(ggplot2)
library(ggrepel)
library(readr)
library(readxl)
library(tidyr)
library(ggvenn)
library(dplyr)
library(purrr)
library(pheatmap)
library(stringr)
#----
setwd("~/Documents/metabolomic_too_heavy_for_git/untargeted/Negative mode/")
###batch_1_negative----

batch1_negative <- read_excel("final_filtered_files/batch_1_Intensity_filtered_notame.xlsx", col_names = FALSE)

batch1_negative <- batch1_negative %>%
  mutate(across(34, ~ ifelse(!is.na(batch1_negative[[3]]) & batch1_negative[[3]] != "", batch1_negative[[3]], .)))
batch1_negative <- batch1_negative[,-c(1:33)]
batch1_negative <- t(batch1_negative)
batch1_negative <- as.data.frame(batch1_negative)
colnames(batch1_negative) <- batch1_negative[1, ]
batch1_negative <- batch1_negative[-1, ]
batch1_negative_row_name <- batch1_negative %>%
  filter(Group != "Blank") %>%
  separate(Group, into = c("Batch","pheno","treatmentA","hpi"), sep = "-") %>%
  mutate(genotype = paste(Batch, pheno, sep = "-")) %>%
  mutate(treatment = paste(Batch, treatmentA, sep = "-")) %>%
  select(FEATURE_ID,genotype,treatment, hpi)
batch1_negative <- merge(batch1_negative_row_name,batch1_negative, by ="FEATURE_ID")
batch1_negative[is.na(batch1_negative)] <- 0 


##batch1_negative_susceptible----
#batch1_negative_susceptible----
# Keep only rows where genotype is NOT "E1-R"
batch1_negative_susceptible <- batch1_negative[batch1_negative$genotype != "E1-R", ]
batch1_negative_susceptible <- batch1_negative_susceptible[, !(names(batch1_negative_susceptible) %in% c("FEATURE_ID", "genotype", "hpi", "Sample_ID", "Injection_order", "Group", "QC"))]
batch1_negative_susceptible <- t(batch1_negative_susceptible)



# 2. Extract genotype and data

genotype <- as.character(unlist(batch1_negative_susceptible[1, ])) # first row = genotype
data <- batch1_negative_susceptible[-1, ]                          # remaining rows = metabolites
data <- apply(data, 2, as.numeric)            # convert to numeric
rownames(data) <- rownames(batch1_negative_susceptible)[-1]

# 3. Compute log2 fold change

mean_R <- rowMeans(data[, genotype == "E1-AF"])
mean_S <- rowMeans(data[, genotype == "E1-C"])

log2FC <- log2(mean_R / mean_S)   


# 4. Compute p-values (t-test)

pvalues <- apply(data, 1, function(x) {
  t.test(x[genotype == "E1-AF"], x[genotype == "E1-C"])$p.value
})


# 5. Adjust p-values (FDR)

padj <- p.adjust(pvalues, method = "fdr")


# 6. Create results_batch1_negative_susceptible table

results_batch1_negative_susceptible <- data.frame(
  metabolite = rownames(data),
  log2FC = log2FC,
  pvalue = pvalues,
  padj = padj
)

results_batch1_negative_susceptible$negLog10P <- -log10(results_batch1_negative_susceptible$padj)


# 7. Volcano plot

# Define significance
results_batch1_negative_susceptible$significant <- "Not Significant"
results_batch1_negative_susceptible$significant[results_batch1_negative_susceptible$log2FC > 1 & results_batch1_negative_susceptible$padj < 0.05] <- "Up"
results_batch1_negative_susceptible$significant[results_batch1_negative_susceptible$log2FC < -1 & results_batch1_negative_susceptible$padj < 0.05] <- "Down"

Volcano_batch1_negative_susceptible <- ggplot(results_batch1_negative_susceptible, aes(x = log2FC, y = negLog10P, color = significant)) +
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
    data = subset(results_batch1_negative_susceptible, significant != "Not Significant"),
    aes(label = metabolite)
  )+
  # Set fixed x and y axis limits
  scale_x_continuous(limits = c(-5, 5)) +  # Adjust min/max for log2FC
  scale_y_continuous(limits = c(0, 20))    # Adjust min/max for -log10(FDR)

results_batch1_negative_susceptible$significant<-as.factor(results_batch1_negative_susceptible$significant)
summary(results_batch1_negative_susceptible)




#early----
# Keep only rows where genotype is NOT "E1-R"
batch1_negative_susceptible <- batch1_negative[batch1_negative$genotype != "E1-R", ]
batch1_negative_susceptible_early <- batch1_negative_susceptible[
  batch1_negative_susceptible$hpi %in% c("3", "24"),]
batch1_negative_susceptible_early <- batch1_negative_susceptible_early[, !(names(batch1_negative_susceptible_early) %in% c("FEATURE_ID", "genotype", "hpi", "Sample_ID", "Injection_order", "Group", "QC"))]
batch1_negative_susceptible_early <- t(batch1_negative_susceptible_early)


# 2. Extract genotype and data

genotype <- as.character(unlist(batch1_negative_susceptible_early[1, ])) # first row = genotype
data <- batch1_negative_susceptible_early[-1, ]                          # remaining rows = metabolites
data <- apply(data, 2, as.numeric)            # convert to numeric
rownames(data) <- rownames(batch1_negative_susceptible_early)[-1]

# 3. Compute log2 fold change

mean_R <- rowMeans(data[, genotype == "E1-AF"])
mean_S <- rowMeans(data[, genotype == "E1-C"])

log2FC <- log2(mean_R / mean_S)   


# 4. Compute p-values (t-test)

pvalues <- apply(data, 1, function(x) {
  t.test(x[genotype == "E1-AF"], x[genotype == "E1-C"])$p.value
})


# 5. Adjust p-values (FDR)

padj <- p.adjust(pvalues, method = "fdr")


# 6. Create results_batch1_negative_susceptible_early table

results_batch1_negative_susceptible_early <- data.frame(
  metabolite = rownames(data),
  log2FC = log2FC,
  pvalue = pvalues,
  padj = padj
)

results_batch1_negative_susceptible_early$negLog10P <- -log10(results_batch1_negative_susceptible_early$padj)


# 7. Volcano plot

# Define significance
results_batch1_negative_susceptible_early$significant <- "Not Significant"
results_batch1_negative_susceptible_early$significant[results_batch1_negative_susceptible_early$log2FC > 1 & results_batch1_negative_susceptible_early$padj < 0.05] <- "Up"
results_batch1_negative_susceptible_early$significant[results_batch1_negative_susceptible_early$log2FC < -1 & results_batch1_negative_susceptible_early$padj < 0.05] <- "Down"

Volcano_batch1_negative_susceptible_early <- ggplot(results_batch1_negative_susceptible_early, aes(x = log2FC, y = negLog10P, color = significant)) +
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
    data = subset(results_batch1_negative_susceptible_early, significant != "Not Significant"),
    aes(label = metabolite)
  )+
  # Set fixed x and y axis limits
  scale_x_continuous(limits = c(-5, 5)) +  # Adjust min/max for log2FC
  scale_y_continuous(limits = c(0, 20))    # Adjust min/max for -log10(FDR)

results_batch1_negative_susceptible_early$significant<-as.factor(results_batch1_negative_susceptible_early$significant)
summary(results_batch1_negative_susceptible_early)





#late----
# Keep only rows where genotype is NOT "E1-R"
batch1_negative_susceptible <- batch1_negative[batch1_negative$genotype != "E1-R", ]
batch1_negative_susceptible_late <- batch1_negative_susceptible[
  batch1_negative_susceptible$hpi %in% c("72", "120"),
]
batch1_negative_susceptible_late <- batch1_negative_susceptible_late[, !(names(batch1_negative_susceptible_late) %in% c("FEATURE_ID", "genotype", "hpi", "Sample_ID", "Injection_order", "Group", "QC"))]
batch1_negative_susceptible_late <- t(batch1_negative_susceptible_late)



# 2. Extract genotype and data

genotype <- as.character(unlist(batch1_negative_susceptible_late[1, ])) # first row = genotype
data <- batch1_negative_susceptible_late[-1, ]                          # remaining rows = metabolites
data <- apply(data, 2, as.numeric)            # convert to numeric
rownames(data) <- rownames(batch1_negative_susceptible_late)[-1]

# 3. Compute log2 fold change

mean_R <- rowMeans(data[, genotype == "E1-AF"])
mean_S <- rowMeans(data[, genotype == "E1-C"])

log2FC <- log2(mean_R / mean_S)   


# 4. Compute p-values (t-test)

pvalues <- apply(data, 1, function(x) {
  t.test(x[genotype == "E1-AF"], x[genotype == "E1-C"])$p.value
})


# 5. Adjust p-values (FDR)

padj <- p.adjust(pvalues, method = "fdr")


# 6. Create results_batch1_negative_susceptible_late table

results_batch1_negative_susceptible_late <- data.frame(
  metabolite = rownames(data),
  log2FC = log2FC,
  pvalue = pvalues,
  padj = padj
)

results_batch1_negative_susceptible_late$negLog10P <- -log10(results_batch1_negative_susceptible_late$padj)


# 7. Volcano plot

# Define significance
results_batch1_negative_susceptible_late$significant <- "Not Significant"
results_batch1_negative_susceptible_late$significant[results_batch1_negative_susceptible_late$log2FC > 1 & results_batch1_negative_susceptible_late$padj < 0.05] <- "Up"
results_batch1_negative_susceptible_late$significant[results_batch1_negative_susceptible_late$log2FC < -1 & results_batch1_negative_susceptible_late$padj < 0.05] <- "Down"

Volcano_batch1_negative_susceptible_late <- ggplot(results_batch1_negative_susceptible_late, aes(x = log2FC, y = negLog10P, color = significant)) +
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
    data = subset(results_batch1_negative_susceptible_late, significant != "Not Significant"),
    aes(label = metabolite)
  )+
  # Set fixed x and y axis limits
  scale_x_continuous(limits = c(-5, 5)) +  # Adjust min/max for log2FC
  scale_y_continuous(limits = c(0, 20))    # Adjust min/max for -log10(FDR)

results_batch1_negative_susceptible_late$significant<-as.factor(results_batch1_negative_susceptible_late$significant)
summary(results_batch1_negative_susceptible_late)




#3hpi----
# Keep only rows where genotype is NOT "E1-R"
batch1_negative_susceptible <- batch1_negative[batch1_negative$genotype != "E1-R", ]
batch1_negative_susceptible_3hpi <- batch1_negative_susceptible[batch1_negative_susceptible$hpi == "3", ]
batch1_negative_susceptible_3hpi <- batch1_negative_susceptible_3hpi[, !(names(batch1_negative_susceptible_3hpi) %in% c("FEATURE_ID", "genotype", "hpi", "Sample_ID", "Injection_order", "Group", "QC"))]
batch1_negative_susceptible_3hpi <- t(batch1_negative_susceptible_3hpi)



# 2. Extract genotype and data

genotype <- as.character(unlist(batch1_negative_susceptible_3hpi[1, ])) # first row = genotype
data <- batch1_negative_susceptible_3hpi[-1, ]                          # remaining rows = metabolites
data <- apply(data, 2, as.numeric)            # convert to numeric
rownames(data) <- rownames(batch1_negative_susceptible_3hpi)[-1]

# 3. Compute log2 fold change

mean_R <- rowMeans(data[, genotype == "E1-AF"])
mean_S <- rowMeans(data[, genotype == "E1-C"])

log2FC <- log2(mean_R / mean_S)   


# 4. Compute p-values (t-test)

pvalues <- apply(data, 1, function(x) {
  t.test(x[genotype == "E1-AF"], x[genotype == "E1-C"])$p.value
})


# 5. Adjust p-values (FDR)

padj <- p.adjust(pvalues, method = "fdr")


# 6. Create results_batch1_negative_susceptible_3hpi table

results_batch1_negative_susceptible_3hpi <- data.frame(
  metabolite = rownames(data),
  log2FC = log2FC,
  pvalue = pvalues,
  padj = padj
)

results_batch1_negative_susceptible_3hpi$negLog10P <- -log10(results_batch1_negative_susceptible_3hpi$padj)


# 7. Volcano plot

# Define significance
results_batch1_negative_susceptible_3hpi$significant <- "Not Significant"
results_batch1_negative_susceptible_3hpi$significant[results_batch1_negative_susceptible_3hpi$log2FC > 1 & results_batch1_negative_susceptible_3hpi$padj < 0.05] <- "Up"
results_batch1_negative_susceptible_3hpi$significant[results_batch1_negative_susceptible_3hpi$log2FC < -1 & results_batch1_negative_susceptible_3hpi$padj < 0.05] <- "Down"

Volcano_batch1_negative_susceptible_3hpi <- ggplot(results_batch1_negative_susceptible_3hpi, aes(x = log2FC, y = negLog10P, color = significant)) +
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
    data = subset(results_batch1_negative_susceptible_3hpi, significant != "Not Significant"),
    aes(label = metabolite)
  )+
  # Set fixed x and y axis limits
  scale_x_continuous(limits = c(-5, 5)) +  # Adjust min/max for log2FC
  scale_y_continuous(limits = c(0, 20))    # Adjust min/max for -log10(FDR)

results_batch1_negative_susceptible_3hpi$significant<-as.factor(results_batch1_negative_susceptible_3hpi$significant)
summary(results_batch1_negative_susceptible_3hpi)




#24hpi----
# Keep only rows where genotype is NOT "E1-R"
batch1_negative_susceptible <- batch1_negative[batch1_negative$genotype != "E1-R", ]
batch1_negative_susceptible_24hpi <- batch1_negative_susceptible[batch1_negative_susceptible$hpi == "24", ]
batch1_negative_susceptible_24hpi <- batch1_negative_susceptible_24hpi[, !(names(batch1_negative_susceptible_24hpi) %in% c("FEATURE_ID", "genotype", "hpi", "Sample_ID", "Injection_order", "Group", "QC"))]
batch1_negative_susceptible_24hpi <- t(batch1_negative_susceptible_24hpi)



# 2. Extract genotype and data

genotype <- as.character(unlist(batch1_negative_susceptible_24hpi[1, ])) # first row = genotype
data <- batch1_negative_susceptible_24hpi[-1, ]                          # remaining rows = metabolites
data <- apply(data, 2, as.numeric)            # convert to numeric
rownames(data) <- rownames(batch1_negative_susceptible_24hpi)[-1]

# 3. Compute log2 fold change

mean_R <- rowMeans(data[, genotype == "E1-AF"])
mean_S <- rowMeans(data[, genotype == "E1-C"])

log2FC <- log2(mean_R / mean_S)   


# 4. Compute p-values (t-test)

pvalues <- apply(data, 1, function(x) {
  t.test(x[genotype == "E1-AF"], x[genotype == "E1-C"])$p.value
})


# 5. Adjust p-values (FDR)

padj <- p.adjust(pvalues, method = "fdr")


# 6. Create results_batch1_negative_susceptible_24hpi table

results_batch1_negative_susceptible_24hpi <- data.frame(
  metabolite = rownames(data),
  log2FC = log2FC,
  pvalue = pvalues,
  padj = padj
)

results_batch1_negative_susceptible_24hpi$negLog10P <- -log10(results_batch1_negative_susceptible_24hpi$padj)


# 7. Volcano plot

# Define significance
results_batch1_negative_susceptible_24hpi$significant <- "Not Significant"
results_batch1_negative_susceptible_24hpi$significant[results_batch1_negative_susceptible_24hpi$log2FC > 1 & results_batch1_negative_susceptible_24hpi$padj < 0.05] <- "Up"
results_batch1_negative_susceptible_24hpi$significant[results_batch1_negative_susceptible_24hpi$log2FC < -1 & results_batch1_negative_susceptible_24hpi$padj < 0.05] <- "Down"

Volcano_batch1_negative_susceptible_24hpi <- ggplot(results_batch1_negative_susceptible_24hpi, aes(x = log2FC, y = negLog10P, color = significant)) +
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
    data = subset(results_batch1_negative_susceptible_24hpi, significant != "Not Significant"),
    aes(label = metabolite)
  )+
  # Set fixed x and y axis limits
  scale_x_continuous(limits = c(-5, 5)) +  # Adjust min/max for log2FC
  scale_y_continuous(limits = c(0, 20))    # Adjust min/max for -log10(FDR)

results_batch1_negative_susceptible_24hpi$significant<-as.factor(results_batch1_negative_susceptible_24hpi$significant)
summary(results_batch1_negative_susceptible_24hpi)




#72hpi----
# Keep only rows where genotype is NOT "E1-R"
batch1_negative_susceptible <- batch1_negative[batch1_negative$genotype != "E1-R", ]
batch1_negative_susceptible_72hpi <- batch1_negative_susceptible[batch1_negative_susceptible$hpi == "72", ]
batch1_negative_susceptible_72hpi <- batch1_negative_susceptible_72hpi[, !(names(batch1_negative_susceptible_72hpi) %in% c("FEATURE_ID", "genotype", "hpi", "Sample_ID", "Injection_order", "Group", "QC"))]
batch1_negative_susceptible_72hpi <- t(batch1_negative_susceptible_72hpi)



# 2. Extract genotype and data

genotype <- as.character(unlist(batch1_negative_susceptible_72hpi[1, ])) # first row = genotype
data <- batch1_negative_susceptible_72hpi[-1, ]                          # remaining rows = metabolites
data <- apply(data, 2, as.numeric)            # convert to numeric
rownames(data) <- rownames(batch1_negative_susceptible_72hpi)[-1]

# 3. Compute log2 fold change

mean_R <- rowMeans(data[, genotype == "E1-AF"])
mean_S <- rowMeans(data[, genotype == "E1-C"])

log2FC <- log2(mean_R / mean_S)   


# 4. Compute p-values (t-test)

pvalues <- apply(data, 1, function(x) {
  t.test(x[genotype == "E1-AF"], x[genotype == "E1-C"])$p.value
})


# 5. Adjust p-values (FDR)

padj <- p.adjust(pvalues, method = "fdr")


# 6. Create results_batch1_negative_susceptible_72hpi table

results_batch1_negative_susceptible_72hpi <- data.frame(
  metabolite = rownames(data),
  log2FC = log2FC,
  pvalue = pvalues,
  padj = padj
)

results_batch1_negative_susceptible_72hpi$negLog10P <- -log10(results_batch1_negative_susceptible_72hpi$padj)


# 7. Volcano plot

# Define significance
results_batch1_negative_susceptible_72hpi$significant <- "Not Significant"
results_batch1_negative_susceptible_72hpi$significant[results_batch1_negative_susceptible_72hpi$log2FC > 1 & results_batch1_negative_susceptible_72hpi$padj < 0.05] <- "Up"
results_batch1_negative_susceptible_72hpi$significant[results_batch1_negative_susceptible_72hpi$log2FC < -1 & results_batch1_negative_susceptible_72hpi$padj < 0.05] <- "Down"

Volcano_batch1_negative_susceptible_72hpi <- ggplot(results_batch1_negative_susceptible_72hpi, aes(x = log2FC, y = negLog10P, color = significant)) +
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
    data = subset(results_batch1_negative_susceptible_72hpi, significant != "Not Significant"),
    aes(label = metabolite)
  )+
  # Set fixed x and y axis limits
  scale_x_continuous(limits = c(-5, 5)) +  # Adjust min/max for log2FC
  scale_y_continuous(limits = c(0, 20))    # Adjust min/max for -log10(FDR)

results_batch1_negative_susceptible_72hpi$significant<-as.factor(results_batch1_negative_susceptible_72hpi$significant)
summary(results_batch1_negative_susceptible_72hpi)

#120hpi----
# Keep only rows where genotype is NOT "E1-R"
batch1_negative_susceptible <- batch1_negative[batch1_negative$genotype != "E1-R", ]
batch1_negative_susceptible_120hpi <- batch1_negative_susceptible[batch1_negative_susceptible$hpi == "120", ]
batch1_negative_susceptible_120hpi <- batch1_negative_susceptible_120hpi[, !(names(batch1_negative_susceptible_120hpi) %in% c("FEATURE_ID", "genotype", "hpi", "Sample_ID", "Injection_order", "Group", "QC"))]
batch1_negative_susceptible_120hpi <- t(batch1_negative_susceptible_120hpi)



# 2. Extract genotype and data

genotype <- as.character(unlist(batch1_negative_susceptible_120hpi[1, ])) # first row = genotype
data <- batch1_negative_susceptible_120hpi[-1, ]                          # remaining rows = metabolites
data <- apply(data, 2, as.numeric)            # convert to numeric
rownames(data) <- rownames(batch1_negative_susceptible_120hpi)[-1]

# 3. Compute log2 fold change

mean_R <- rowMeans(data[, genotype == "E1-AF"])
mean_S <- rowMeans(data[, genotype == "E1-C"])

log2FC <- log2(mean_R / mean_S)   


# 4. Compute p-values (t-test)

pvalues <- apply(data, 1, function(x) {
  t.test(x[genotype == "E1-AF"], x[genotype == "E1-C"])$p.value
})


# 5. Adjust p-values (FDR)

padj <- p.adjust(pvalues, method = "fdr")


# 6. Create results_batch1_negative_susceptible_120hpi table

results_batch1_negative_susceptible_120hpi <- data.frame(
  metabolite = rownames(data),
  log2FC = log2FC,
  pvalue = pvalues,
  padj = padj
)

results_batch1_negative_susceptible_120hpi$negLog10P <- -log10(results_batch1_negative_susceptible_120hpi$padj)


# 7. Volcano plot

# Define significance
results_batch1_negative_susceptible_120hpi$significant <- "Not Significant"
results_batch1_negative_susceptible_120hpi$significant[results_batch1_negative_susceptible_120hpi$log2FC > 1 & results_batch1_negative_susceptible_120hpi$padj < 0.05] <- "Up"
results_batch1_negative_susceptible_120hpi$significant[results_batch1_negative_susceptible_120hpi$log2FC < -1 & results_batch1_negative_susceptible_120hpi$padj < 0.05] <- "Down"

Volcano_batch1_negative_susceptible_120hpi <- ggplot(results_batch1_negative_susceptible_120hpi, aes(x = log2FC, y = negLog10P, color = significant)) +
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
    data = subset(results_batch1_negative_susceptible_120hpi, significant != "Not Significant"),
    aes(label = metabolite)
  )+
  # Set fixed x and y axis limits
  scale_x_continuous(limits = c(-5, 5)) +  # Adjust min/max for log2FC
  scale_y_continuous(limits = c(0, 20))    # Adjust min/max for -log10(FDR)

results_batch1_negative_susceptible_120hpi$significant<-as.factor(results_batch1_negative_susceptible_120hpi$significant)
summary(results_batch1_negative_susceptible_120hpi)







##batch1_negative_resistant----
#batch1_negative_resistant----
# Keep only rows where genotype is NOT "E1-S"
batch1_negative_resistant <- batch1_negative[batch1_negative$genotype != "E1-S", ]
batch1_negative_resistant <- batch1_negative_resistant[, !(names(batch1_negative_resistant) %in% c("FEATURE_ID", "genotype", "hpi", "Sample_ID", "Injection_order", "Group", "QC"))]
batch1_negative_resistant <- t(batch1_negative_resistant)



# 2. Extract genotype and data

genotype <- as.character(unlist(batch1_negative_resistant[1, ])) # first row = genotype
data <- batch1_negative_resistant[-1, ]                          # remaining rows = metabolites
data <- apply(data, 2, as.numeric)            # convert to numeric
rownames(data) <- rownames(batch1_negative_resistant)[-1]

# 3. Compute log2 fold change

mean_R <- rowMeans(data[, genotype == "E1-AF"])
mean_S <- rowMeans(data[, genotype == "E1-C"])

log2FC <- log2(mean_R / mean_S)   


# 4. Compute p-values (t-test)

pvalues <- apply(data, 1, function(x) {
  t.test(x[genotype == "E1-AF"], x[genotype == "E1-C"])$p.value
})


# 5. Adjust p-values (FDR)

padj <- p.adjust(pvalues, method = "fdr")


# 6. Create results_batch1_negative_resistant table

results_batch1_negative_resistant <- data.frame(
  metabolite = rownames(data),
  log2FC = log2FC,
  pvalue = pvalues,
  padj = padj
)

results_batch1_negative_resistant$negLog10P <- -log10(results_batch1_negative_resistant$padj)


# 7. Volcano plot

# Define significance
results_batch1_negative_resistant$significant <- "Not Significant"
results_batch1_negative_resistant$significant[results_batch1_negative_resistant$log2FC > 1 & results_batch1_negative_resistant$padj < 0.05] <- "Up"
results_batch1_negative_resistant$significant[results_batch1_negative_resistant$log2FC < -1 & results_batch1_negative_resistant$padj < 0.05] <- "Down"

Volcano_batch1_negative_resistant <- ggplot(results_batch1_negative_resistant, aes(x = log2FC, y = negLog10P, color = significant)) +
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
    data = subset(results_batch1_negative_resistant, significant != "Not Significant"),
    aes(label = metabolite)
  )+
  # Set fixed x and y axis limits
  scale_x_continuous(limits = c(-5, 5)) +  # Adjust min/max for log2FC
  scale_y_continuous(limits = c(0, 20))    # Adjust min/max for -log10(FDR)

results_batch1_negative_resistant$significant<-as.factor(results_batch1_negative_resistant$significant)
summary(results_batch1_negative_resistant)




#early----
# Keep only rows where genotype is NOT "E1-S"
batch1_negative_resistant <- batch1_negative[batch1_negative$genotype != "E1-S", ]
batch1_negative_resistant_early <- batch1_negative_resistant[
  batch1_negative_resistant$hpi %in% c("3", "24"),]
batch1_negative_resistant_early <- batch1_negative_resistant_early[, !(names(batch1_negative_resistant_early) %in% c("FEATURE_ID", "genotype", "hpi", "Sample_ID", "Injection_order", "Group", "QC"))]
batch1_negative_resistant_early <- t(batch1_negative_resistant_early)



# 2. Extract genotype and data

genotype <- as.character(unlist(batch1_negative_resistant_early[1, ])) # first row = genotype
data <- batch1_negative_resistant_early[-1, ]                          # remaining rows = metabolites
data <- apply(data, 2, as.numeric)            # convert to numeric
rownames(data) <- rownames(batch1_negative_resistant_early)[-1]

# 3. Compute log2 fold change

mean_R <- rowMeans(data[, genotype == "E1-AF"])
mean_S <- rowMeans(data[, genotype == "E1-C"])

log2FC <- log2(mean_R / mean_S)   


# 4. Compute p-values (t-test)

pvalues <- apply(data, 1, function(x) {
  t.test(x[genotype == "E1-AF"], x[genotype == "E1-C"])$p.value
})


# 5. Adjust p-values (FDR)

padj <- p.adjust(pvalues, method = "fdr")


# 6. Create results_batch1_negative_resistant_early table

results_batch1_negative_resistant_early <- data.frame(
  metabolite = rownames(data),
  log2FC = log2FC,
  pvalue = pvalues,
  padj = padj
)

results_batch1_negative_resistant_early$negLog10P <- -log10(results_batch1_negative_resistant_early$padj)


# 7. Volcano plot

# Define significance
results_batch1_negative_resistant_early$significant <- "Not Significant"
results_batch1_negative_resistant_early$significant[results_batch1_negative_resistant_early$log2FC > 1 & results_batch1_negative_resistant_early$padj < 0.05] <- "Up"
results_batch1_negative_resistant_early$significant[results_batch1_negative_resistant_early$log2FC < -1 & results_batch1_negative_resistant_early$padj < 0.05] <- "Down"

Volcano_batch1_negative_resistant_early <- ggplot(results_batch1_negative_resistant_early, aes(x = log2FC, y = negLog10P, color = significant)) +
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
    data = subset(results_batch1_negative_resistant_early, significant != "Not Significant"),
    aes(label = metabolite)
  )+
  # Set fixed x and y axis limits
  scale_x_continuous(limits = c(-5, 5)) +  # Adjust min/max for log2FC
  scale_y_continuous(limits = c(0, 20))    # Adjust min/max for -log10(FDR)

results_batch1_negative_resistant_early$significant<-as.factor(results_batch1_negative_resistant_early$significant)
summary(results_batch1_negative_resistant_early)





#late----
# Keep only rows where genotype is NOT "E1-S"
batch1_negative_resistant <- batch1_negative[batch1_negative$genotype != "E1-S", ]
batch1_negative_resistant_late <- batch1_negative_resistant[
  batch1_negative_resistant$hpi %in% c("72", "120"),
]
batch1_negative_resistant_late <- batch1_negative_resistant_late[, !(names(batch1_negative_resistant_late) %in% c("FEATURE_ID", "genotype", "hpi", "Sample_ID", "Injection_order", "Group", "QC"))]
batch1_negative_resistant_late <- t(batch1_negative_resistant_late)



# 2. Extract genotype and data

genotype <- as.character(unlist(batch1_negative_resistant_late[1, ])) # first row = genotype
data <- batch1_negative_resistant_late[-1, ]                          # remaining rows = metabolites
data <- apply(data, 2, as.numeric)            # convert to numeric
rownames(data) <- rownames(batch1_negative_resistant_late)[-1]

# 3. Compute log2 fold change

mean_R <- rowMeans(data[, genotype == "E1-AF"])
mean_S <- rowMeans(data[, genotype == "E1-C"])

log2FC <- log2(mean_R / mean_S)   


# 4. Compute p-values (t-test)

pvalues <- apply(data, 1, function(x) {
  t.test(x[genotype == "E1-AF"], x[genotype == "E1-C"])$p.value
})


# 5. Adjust p-values (FDR)

padj <- p.adjust(pvalues, method = "fdr")


# 6. Create results_batch1_negative_resistant_late table

results_batch1_negative_resistant_late <- data.frame(
  metabolite = rownames(data),
  log2FC = log2FC,
  pvalue = pvalues,
  padj = padj
)

results_batch1_negative_resistant_late$negLog10P <- -log10(results_batch1_negative_resistant_late$padj)


# 7. Volcano plot

# Define significance
results_batch1_negative_resistant_late$significant <- "Not Significant"
results_batch1_negative_resistant_late$significant[results_batch1_negative_resistant_late$log2FC > 1 & results_batch1_negative_resistant_late$padj < 0.05] <- "Up"
results_batch1_negative_resistant_late$significant[results_batch1_negative_resistant_late$log2FC < -1 & results_batch1_negative_resistant_late$padj < 0.05] <- "Down"

Volcano_batch1_negative_resistant_late <- ggplot(results_batch1_negative_resistant_late, aes(x = log2FC, y = negLog10P, color = significant)) +
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
    data = subset(results_batch1_negative_resistant_late, significant != "Not Significant"),
    aes(label = metabolite)
  )+
  # Set fixed x and y axis limits
  scale_x_continuous(limits = c(-5, 5)) +  # Adjust min/max for log2FC
  scale_y_continuous(limits = c(0, 20))    # Adjust min/max for -log10(FDR)

results_batch1_negative_resistant_late$significant<-as.factor(results_batch1_negative_resistant_late$significant)
summary(results_batch1_negative_resistant_late)




#3hpi----
# Keep only rows where genotype is NOT "E1-S"
batch1_negative_resistant <- batch1_negative[batch1_negative$genotype != "E1-S", ]
batch1_negative_resistant_3hpi <- batch1_negative_resistant[batch1_negative_resistant$hpi == "3", ]
batch1_negative_resistant_3hpi <- batch1_negative_resistant_3hpi[, !(names(batch1_negative_resistant_3hpi) %in% c("FEATURE_ID", "genotype", "hpi", "Sample_ID", "Injection_order", "Group", "QC"))]
batch1_negative_resistant_3hpi <- t(batch1_negative_resistant_3hpi)



# 2. Extract genotype and data

genotype <- as.character(unlist(batch1_negative_resistant_3hpi[1, ])) # first row = genotype
data <- batch1_negative_resistant_3hpi[-1, ]                          # remaining rows = metabolites
data <- apply(data, 2, as.numeric)            # convert to numeric
rownames(data) <- rownames(batch1_negative_resistant_3hpi)[-1]

# 3. Compute log2 fold change

mean_R <- rowMeans(data[, genotype == "E1-AF"])
mean_S <- rowMeans(data[, genotype == "E1-C"])

log2FC <- log2(mean_R / mean_S)   


# 4. Compute p-values (t-test)

pvalues <- apply(data, 1, function(x) {
  t.test(x[genotype == "E1-AF"], x[genotype == "E1-C"])$p.value
})


# 5. Adjust p-values (FDR)

padj <- p.adjust(pvalues, method = "fdr")


# 6. Create results_batch1_negative_resistant_3hpi table

results_batch1_negative_resistant_3hpi <- data.frame(
  metabolite = rownames(data),
  log2FC = log2FC,
  pvalue = pvalues,
  padj = padj
)

results_batch1_negative_resistant_3hpi$negLog10P <- -log10(results_batch1_negative_resistant_3hpi$padj)


# 7. Volcano plot

# Define significance
results_batch1_negative_resistant_3hpi$significant <- "Not Significant"
results_batch1_negative_resistant_3hpi$significant[results_batch1_negative_resistant_3hpi$log2FC > 1 & results_batch1_negative_resistant_3hpi$padj < 0.05] <- "Up"
results_batch1_negative_resistant_3hpi$significant[results_batch1_negative_resistant_3hpi$log2FC < -1 & results_batch1_negative_resistant_3hpi$padj < 0.05] <- "Down"

Volcano_batch1_negative_resistant_3hpi <- ggplot(results_batch1_negative_resistant_3hpi, aes(x = log2FC, y = negLog10P, color = significant)) +
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
    data = subset(results_batch1_negative_resistant_3hpi, significant != "Not Significant"),
    aes(label = metabolite)
  )+
  # Set fixed x and y axis limits
  scale_x_continuous(limits = c(-5, 5)) +  # Adjust min/max for log2FC
  scale_y_continuous(limits = c(0, 20))    # Adjust min/max for -log10(FDR)

results_batch1_negative_resistant_3hpi$significant<-as.factor(results_batch1_negative_resistant_3hpi$significant)
summary(results_batch1_negative_resistant_3hpi)




#24hpi----
# Keep only rows where genotype is NOT "E1-S"
batch1_negative_resistant <- batch1_negative[batch1_negative$genotype != "E1-S", ]
batch1_negative_resistant_24hpi <- batch1_negative_resistant[batch1_negative_resistant$hpi == "24", ]
batch1_negative_resistant_24hpi <- batch1_negative_resistant_24hpi[, !(names(batch1_negative_resistant_24hpi) %in% c("FEATURE_ID", "genotype", "hpi", "Sample_ID", "Injection_order", "Group", "QC"))]
batch1_negative_resistant_24hpi <- t(batch1_negative_resistant_24hpi)



# 2. Extract genotype and data

genotype <- as.character(unlist(batch1_negative_resistant_24hpi[1, ])) # first row = genotype
data <- batch1_negative_resistant_24hpi[-1, ]                          # remaining rows = metabolites
data <- apply(data, 2, as.numeric)            # convert to numeric
rownames(data) <- rownames(batch1_negative_resistant_24hpi)[-1]

# 3. Compute log2 fold change

mean_R <- rowMeans(data[, genotype == "E1-AF"])
mean_S <- rowMeans(data[, genotype == "E1-C"])

log2FC <- log2(mean_R / mean_S)   


# 4. Compute p-values (t-test)

pvalues <- apply(data, 1, function(x) {
  t.test(x[genotype == "E1-AF"], x[genotype == "E1-C"])$p.value
})


# 5. Adjust p-values (FDR)

padj <- p.adjust(pvalues, method = "fdr")


# 6. Create results_batch1_negative_resistant_24hpi table

results_batch1_negative_resistant_24hpi <- data.frame(
  metabolite = rownames(data),
  log2FC = log2FC,
  pvalue = pvalues,
  padj = padj
)

results_batch1_negative_resistant_24hpi$negLog10P <- -log10(results_batch1_negative_resistant_24hpi$padj)


# 7. Volcano plot

# Define significance
results_batch1_negative_resistant_24hpi$significant <- "Not Significant"
results_batch1_negative_resistant_24hpi$significant[results_batch1_negative_resistant_24hpi$log2FC > 1 & results_batch1_negative_resistant_24hpi$padj < 0.05] <- "Up"
results_batch1_negative_resistant_24hpi$significant[results_batch1_negative_resistant_24hpi$log2FC < -1 & results_batch1_negative_resistant_24hpi$padj < 0.05] <- "Down"

Volcano_batch1_negative_resistant_24hpi <- ggplot(results_batch1_negative_resistant_24hpi, aes(x = log2FC, y = negLog10P, color = significant)) +
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
    data = subset(results_batch1_negative_resistant_24hpi, significant != "Not Significant"),
    aes(label = metabolite)
  )+
  # Set fixed x and y axis limits
  scale_x_continuous(limits = c(-5, 5)) +  # Adjust min/max for log2FC
  scale_y_continuous(limits = c(0, 20))    # Adjust min/max for -log10(FDR)

results_batch1_negative_resistant_24hpi$significant<-as.factor(results_batch1_negative_resistant_24hpi$significant)
summary(results_batch1_negative_resistant_24hpi)




#72hpi----
# Keep only rows where genotype is NOT "E1-S"
batch1_negative_resistant <- batch1_negative[batch1_negative$genotype != "E1-S", ]
batch1_negative_resistant_72hpi <- batch1_negative_resistant[batch1_negative_resistant$hpi == "72", ]
batch1_negative_resistant_72hpi <- batch1_negative_resistant_72hpi[, !(names(batch1_negative_resistant_72hpi) %in% c("FEATURE_ID", "genotype", "hpi", "Sample_ID", "Injection_order", "Group", "QC"))]
batch1_negative_resistant_72hpi <- t(batch1_negative_resistant_72hpi)



# 2. Extract genotype and data

genotype <- as.character(unlist(batch1_negative_resistant_72hpi[1, ])) # first row = genotype
data <- batch1_negative_resistant_72hpi[-1, ]                          # remaining rows = metabolites
data <- apply(data, 2, as.numeric)            # convert to numeric
rownames(data) <- rownames(batch1_negative_resistant_72hpi)[-1]

# 3. Compute log2 fold change

mean_R <- rowMeans(data[, genotype == "E1-AF"])
mean_S <- rowMeans(data[, genotype == "E1-C"])

log2FC <- log2(mean_R / mean_S)   


# 4. Compute p-values (t-test)

pvalues <- apply(data, 1, function(x) {
  t.test(x[genotype == "E1-AF"], x[genotype == "E1-C"])$p.value
})


# 5. Adjust p-values (FDR)

padj <- p.adjust(pvalues, method = "fdr")


# 6. Create results_batch1_negative_resistant_72hpi table

results_batch1_negative_resistant_72hpi <- data.frame(
  metabolite = rownames(data),
  log2FC = log2FC,
  pvalue = pvalues,
  padj = padj
)

results_batch1_negative_resistant_72hpi$negLog10P <- -log10(results_batch1_negative_resistant_72hpi$padj)


# 7. Volcano plot

# Define significance
results_batch1_negative_resistant_72hpi$significant <- "Not Significant"
results_batch1_negative_resistant_72hpi$significant[results_batch1_negative_resistant_72hpi$log2FC > 1 & results_batch1_negative_resistant_72hpi$padj < 0.05] <- "Up"
results_batch1_negative_resistant_72hpi$significant[results_batch1_negative_resistant_72hpi$log2FC < -1 & results_batch1_negative_resistant_72hpi$padj < 0.05] <- "Down"

Volcano_batch1_negative_resistant_72hpi <- ggplot(results_batch1_negative_resistant_72hpi, aes(x = log2FC, y = negLog10P, color = significant)) +
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
    data = subset(results_batch1_negative_resistant_72hpi, significant != "Not Significant"),
    aes(label = metabolite)
  )+
  # Set fixed x and y axis limits
  scale_x_continuous(limits = c(-5, 5)) +  # Adjust min/max for log2FC
  scale_y_continuous(limits = c(0, 20))    # Adjust min/max for -log10(FDR)

results_batch1_negative_resistant_72hpi$significant<-as.factor(results_batch1_negative_resistant_72hpi$significant)
summary(results_batch1_negative_resistant_72hpi)

#120hpi----
# Keep only rows where genotype is NOT "E1-S"
batch1_negative_resistant <- batch1_negative[batch1_negative$genotype != "E1-S", ]
batch1_negative_resistant_120hpi <- batch1_negative_resistant[batch1_negative_resistant$hpi == "120", ]
batch1_negative_resistant_120hpi <- batch1_negative_resistant_120hpi[, !(names(batch1_negative_resistant_120hpi) %in% c("FEATURE_ID", "genotype", "hpi", "Sample_ID", "Injection_order", "Group", "QC"))]
batch1_negative_resistant_120hpi <- t(batch1_negative_resistant_120hpi)



# 2. Extract genotype and data

genotype <- as.character(unlist(batch1_negative_resistant_120hpi[1, ])) # first row = genotype
data <- batch1_negative_resistant_120hpi[-1, ]                          # remaining rows = metabolites
data <- apply(data, 2, as.numeric)            # convert to numeric
rownames(data) <- rownames(batch1_negative_resistant_120hpi)[-1]

# 3. Compute log2 fold change

mean_R <- rowMeans(data[, genotype == "E1-AF"])
mean_S <- rowMeans(data[, genotype == "E1-C"])

log2FC <- log2(mean_R / mean_S)   


# 4. Compute p-values (t-test)

pvalues <- apply(data, 1, function(x) {
  t.test(x[genotype == "E1-AF"], x[genotype == "E1-C"])$p.value
})


# 5. Adjust p-values (FDR)

padj <- p.adjust(pvalues, method = "fdr")


# 6. Create results_batch1_negative_resistant_120hpi table

results_batch1_negative_resistant_120hpi <- data.frame(
  metabolite = rownames(data),
  log2FC = log2FC,
  pvalue = pvalues,
  padj = padj
)

results_batch1_negative_resistant_120hpi$negLog10P <- -log10(results_batch1_negative_resistant_120hpi$padj)


# 7. Volcano plot

# Define significance
results_batch1_negative_resistant_120hpi$significant <- "Not Significant"
results_batch1_negative_resistant_120hpi$significant[results_batch1_negative_resistant_120hpi$log2FC > 1 & results_batch1_negative_resistant_120hpi$padj < 0.05] <- "Up"
results_batch1_negative_resistant_120hpi$significant[results_batch1_negative_resistant_120hpi$log2FC < -1 & results_batch1_negative_resistant_120hpi$padj < 0.05] <- "Down"

Volcano_batch1_negative_resistant_120hpi <- ggplot(results_batch1_negative_resistant_120hpi, aes(x = log2FC, y = negLog10P, color = significant)) +
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
    data = subset(results_batch1_negative_resistant_120hpi, significant != "Not Significant"),
    aes(label = metabolite)
  )+
  # Set fixed x and y axis limits
  scale_x_continuous(limits = c(-5, 5)) +  # Adjust min/max for log2FC
  scale_y_continuous(limits = c(0, 20))    # Adjust min/max for -log10(FDR)

results_batch1_negative_resistant_120hpi$significant<-as.factor(results_batch1_negative_resistant_120hpi$significant)
summary(results_batch1_negative_resistant_120hpi)





#####

results_batch1_negative_susceptible<-results_batch1_negative_susceptible[,c(1,6)]
colnames(results_batch1_negative_susceptible)[colnames(results_batch1_negative_susceptible) == "significant"] <- "batch1_negative_susceptible"
results_batch1_negative_resistant<-results_batch1_negative_resistant[,c(1,6)]
colnames(results_batch1_negative_resistant)[colnames(results_batch1_negative_resistant) == "significant"] <- "batch1_negative_resistant"

results_batch1_negative_susceptible_early<-results_batch1_negative_susceptible_early[,c(1,6)]
colnames(results_batch1_negative_susceptible_early)[colnames(results_batch1_negative_susceptible_early) == "significant"] <- "batch1_negative_susceptible_early"

results_batch1_negative_resistant_early<-results_batch1_negative_resistant_early[,c(1,6)]
colnames(results_batch1_negative_resistant_early)[colnames(results_batch1_negative_resistant_early) == "significant"] <- "batch1_negative_resistant_early"


results_batch1_negative_susceptible_late<-results_batch1_negative_susceptible_late[,c(1,6)]
colnames(results_batch1_negative_susceptible_late)[colnames(results_batch1_negative_susceptible_late) == "significant"] <- "batch1_negative_susceptible_late"
results_batch1_negative_resistant_late<-results_batch1_negative_resistant_late[,c(1,6)]
colnames(results_batch1_negative_resistant_late)[colnames(results_batch1_negative_resistant_late) == "significant"] <- "batch1_negative_resistant_late"
results_batch1_negative_susceptible_3hpi<-results_batch1_negative_susceptible_3hpi[,c(1,6)]
colnames(results_batch1_negative_susceptible_3hpi)[colnames(results_batch1_negative_susceptible_3hpi) == "significant"] <- "batch1_negative_susceptible_3hpi"
results_batch1_negative_resistant_3hpi<-results_batch1_negative_resistant_3hpi[,c(1,6)]
colnames(results_batch1_negative_resistant_3hpi)[colnames(results_batch1_negative_resistant_3hpi) == "significant"] <- "batch1_negative_resistant_3hpi"
results_batch1_negative_susceptible_24hpi<-results_batch1_negative_susceptible_24hpi[,c(1,6)]
colnames(results_batch1_negative_susceptible_24hpi)[colnames(results_batch1_negative_susceptible_24hpi) == "significant"] <- "batch1_negative_susceptible_24hpi"
results_batch1_negative_resistant_24hpi<-results_batch1_negative_resistant_24hpi[,c(1,6)]
colnames(results_batch1_negative_resistant_24hpi)[colnames(results_batch1_negative_resistant_24hpi) == "significant"] <- "batch1_negative_resistant_24hpi"
results_batch1_negative_susceptible_72hpi<-results_batch1_negative_susceptible_72hpi[,c(1,6)]
colnames(results_batch1_negative_susceptible_72hpi)[colnames(results_batch1_negative_susceptible_72hpi) == "significant"] <- "batch1_negative_susceptible_72hpi"
results_batch1_negative_resistant_72hpi<-results_batch1_negative_resistant_72hpi[,c(1,6)]
colnames(results_batch1_negative_resistant_72hpi)[colnames(results_batch1_negative_resistant_72hpi) == "significant"] <- "batch1_negative_resistant_72hpi"
results_batch1_negative_susceptible_120hpi<-results_batch1_negative_susceptible_120hpi[,c(1,6)]
colnames(results_batch1_negative_susceptible_120hpi)[colnames(results_batch1_negative_susceptible_120hpi) == "significant"] <- "batch1_negative_susceptible_120hpi"
results_batch1_negative_resistant_120hpi<-results_batch1_negative_resistant_120hpi[,c(1,6)]
colnames(results_batch1_negative_resistant_120hpi)[colnames(results_batch1_negative_resistant_120hpi) == "significant"] <- "batch1_negative_resistant_120hpi"




all_results <- list(
  results_batch1_negative_susceptible,
  results_batch1_negative_resistant,
  results_batch1_negative_susceptible_early,
  results_batch1_negative_resistant_early,
  results_batch1_negative_susceptible_late,
  results_batch1_negative_resistant_late,
  results_batch1_negative_susceptible_3hpi,
  results_batch1_negative_resistant_3hpi,
  results_batch1_negative_susceptible_24hpi,
  results_batch1_negative_resistant_24hpi,
  results_batch1_negative_susceptible_72hpi,
  results_batch1_negative_resistant_72hpi,
  results_batch1_negative_susceptible_120hpi,
  results_batch1_negative_resistant_120hpi
)

merged_results_results_batch1_negative <- reduce(all_results, full_join, by = "metabolite")
merged_results_results_batch1_negative <- merged_results_results_batch1_negative[
  rowSums(merged_results_results_batch1_negative[, -1] != "Not Significant", na.rm = TRUE) > 0,
]

merged_results_results_batch1_negative$metabolite <- str_c(
  merged_results_results_batch1_negative$metabolite,
  "neg"
)


heatmap_batch1_negative <- merged_results_results_batch1_negative


heatmap_batch1_negative <- as.data.frame(lapply(heatmap_batch1_negative, as.character))
rownames(heatmap_batch1_negative) <- heatmap_batch1_negative$metabolite
heatmap_batch1_negative <- heatmap_batch1_negative[, -1]

heatmap_batch1_negative[] <- lapply(heatmap_batch1_negative, function(x) {
  x[x == "Up"] <- 1
  x[x == "Down"] <- -1
  x[x == "Not Significant"] <- 0
  return(as.numeric(x))
})

# Convert to matrix
heatmap_matrix <- as.matrix(heatmap_batch1_negative)

# Convert to numeric
mode(heatmap_matrix) <- "numeric"

# compute row sums (total signal per metabolite)
row_score <- rowSums(heatmap_matrix, na.rm = TRUE)

# order from highest to lowest
heatmap_matrix_batch1_negative <- heatmap_matrix[order(row_score, decreasing = TRUE), ]


# Plot heatmap
pheatmap(
  heatmap_matrix_batch1_negative,
  color = c("blue", "white", "red"),
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  border_color = "grey80"
)

###batch_2_negative----

batch2_negative <- read_excel("final_filtered_files/batch_2_Intensity_filtered_notame.xlsx", col_names = FALSE)

batch2_negative <- batch2_negative %>%
  mutate(across(34, ~ ifelse(!is.na(batch2_negative[[3]]) & batch2_negative[[3]] != "", batch2_negative[[3]], .)))
batch2_negative <- batch2_negative[,-c(1:33)]
batch2_negative <- t(batch2_negative)
batch2_negative <- as.data.frame(batch2_negative)
colnames(batch2_negative) <- batch2_negative[1, ]
batch2_negative <- batch2_negative[-1, ]
batch2_negative_row_name <- batch2_negative %>%
  filter(Group != "Blank") %>%
  separate(Group, into = c("Batch","pheno","treatmentA","hpi"), sep = "-") %>%
  mutate(genotype = paste(Batch, pheno, sep = "-")) %>%
  mutate(treatment = paste(Batch, treatmentA, sep = "-")) %>%
  select(FEATURE_ID,genotype,treatment, hpi)
batch2_negative <- merge(batch2_negative_row_name,batch2_negative, by ="FEATURE_ID")
batch2_negative[is.na(batch2_negative)] <- 0 


##batch2_negative_susceptible----
#batch2_negative_susceptible----
# Keep only rows where genotype is NOT "E2-R"
batch2_negative_susceptible <- batch2_negative[batch2_negative$genotype != "E2-R", ]
batch2_negative_susceptible <- batch2_negative_susceptible[, !(names(batch2_negative_susceptible) %in% c("FEATURE_ID", "genotype", "hpi", "Sample_ID", "Injection_order", "Group", "QC"))]
batch2_negative_susceptible <- t(batch2_negative_susceptible)



# 2. Extract genotype and data

genotype <- as.character(unlist(batch2_negative_susceptible[1, ])) # first row = genotype
data <- batch2_negative_susceptible[-1, ]                          # remaining rows = metabolites
data <- apply(data, 2, as.numeric)            # convert to numeric
rownames(data) <- rownames(batch2_negative_susceptible)[-1]

# 3. Compute log2 fold change

mean_R <- rowMeans(data[, genotype == "E2-AP"])
mean_S <- rowMeans(data[, genotype == "E2-C"])

log2FC <- log2(mean_R / mean_S)   


# 4. Compute p-values (t-test)

pvalues <- apply(data, 1, function(x) {
  t.test(x[genotype == "E2-AP"], x[genotype == "E2-C"])$p.value
})


# 5. Adjust p-values (FDR)

padj <- p.adjust(pvalues, method = "fdr")


# 6. Create results_batch2_negative_susceptible table

results_batch2_negative_susceptible <- data.frame(
  metabolite = rownames(data),
  log2FC = log2FC,
  pvalue = pvalues,
  padj = padj
)

results_batch2_negative_susceptible$negLog10P <- -log10(results_batch2_negative_susceptible$padj)


# 7. Volcano plot

# Define significance
results_batch2_negative_susceptible$significant <- "Not Significant"
results_batch2_negative_susceptible$significant[results_batch2_negative_susceptible$log2FC > 1 & results_batch2_negative_susceptible$padj < 0.05] <- "Up"
results_batch2_negative_susceptible$significant[results_batch2_negative_susceptible$log2FC < -1 & results_batch2_negative_susceptible$padj < 0.05] <- "Down"

Volcano_batch2_negative_susceptible <- ggplot(results_batch2_negative_susceptible, aes(x = log2FC, y = negLog10P, color = significant)) +
  geom_point(size = 3, alpha = 0.8) +
  scale_color_manual(values = c( "grey", "red")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  theme_minimal() +
  labs(
    x = "Log2 Fold Change (E2-AP / E2-C)",
    y = "-Log10(FDR)",
    title = "Volcano Plot"
  ) +
  geom_text_repel(
    data = subset(results_batch2_negative_susceptible, significant != "Not Significant"),
    aes(label = metabolite)
  )+
  # Set fixed x and y axis limits
  scale_x_continuous(limits = c(-5, 5)) +  # Adjust min/max for log2FC
  scale_y_continuous(limits = c(0, 5))    # Adjust min/max for -log10(FDR)

results_batch2_negative_susceptible$significant<-as.factor(results_batch2_negative_susceptible$significant)
summary(results_batch2_negative_susceptible)




#early----
# Keep only rows where genotype is NOT "E2-R"
batch2_negative_susceptible <- batch2_negative[batch2_negative$genotype != "E2-R", ]
batch2_negative_susceptible_early <- batch2_negative_susceptible[
  batch2_negative_susceptible$hpi %in% c("3", "24"),]
batch2_negative_susceptible_early <- batch2_negative_susceptible_early[, !(names(batch2_negative_susceptible_early) %in% c("FEATURE_ID", "genotype", "hpi", "Sample_ID", "Injection_order", "Group", "QC"))]
batch2_negative_susceptible_early <- t(batch2_negative_susceptible_early)


# 2. Extract genotype and data

genotype <- as.character(unlist(batch2_negative_susceptible_early[1, ])) # first row = genotype
data <- batch2_negative_susceptible_early[-1, ]                          # remaining rows = metabolites
data <- apply(data, 2, as.numeric)            # convert to numeric
rownames(data) <- rownames(batch2_negative_susceptible_early)[-1]

# 3. Compute log2 fold change

mean_R <- rowMeans(data[, genotype == "E2-AP"])
mean_S <- rowMeans(data[, genotype == "E2-C"])

log2FC <- log2(mean_R / mean_S)   


# 4. Compute p-values (t-test)

pvalues <- apply(data, 1, function(x) {
  t.test(x[genotype == "E2-AP"], x[genotype == "E2-C"])$p.value
})


# 5. Adjust p-values (FDR)

padj <- p.adjust(pvalues, method = "fdr")


# 6. Create results_batch2_negative_susceptible_early table

results_batch2_negative_susceptible_early <- data.frame(
  metabolite = rownames(data),
  log2FC = log2FC,
  pvalue = pvalues,
  padj = padj
)

results_batch2_negative_susceptible_early$negLog10P <- -log10(results_batch2_negative_susceptible_early$padj)


# 7. Volcano plot

# Define significance
results_batch2_negative_susceptible_early$significant <- "Not Significant"
results_batch2_negative_susceptible_early$significant[results_batch2_negative_susceptible_early$log2FC > 1 & results_batch2_negative_susceptible_early$padj < 0.05] <- "Up"
results_batch2_negative_susceptible_early$significant[results_batch2_negative_susceptible_early$log2FC < -1 & results_batch2_negative_susceptible_early$padj < 0.05] <- "Down"

Volcano_batch2_negative_susceptible_early <- ggplot(results_batch2_negative_susceptible_early, aes(x = log2FC, y = negLog10P, color = significant)) +
  geom_point(size = 3, alpha = 0.8) +
  scale_color_manual(values = c("blue", "grey", "red")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  theme_minimal() +
  labs(
    x = "Log2 Fold Change (E2-AP / E2-C)",
    y = "-Log10(FDR)",
    title = "Volcano Plot"
  ) +
  geom_text_repel(
    data = subset(results_batch2_negative_susceptible_early, significant != "Not Significant"),
    aes(label = metabolite)
  )+
  # Set fixed x and y axis limits
  scale_x_continuous(limits = c(-5, 5)) +  # Adjust min/max for log2FC
  scale_y_continuous(limits = c(0, 20))    # Adjust min/max for -log10(FDR)

results_batch2_negative_susceptible_early$significant<-as.factor(results_batch2_negative_susceptible_early$significant)
summary(results_batch2_negative_susceptible_early)





#late----
# Keep only rows where genotype is NOT "E2-R"
batch2_negative_susceptible <- batch2_negative[batch2_negative$genotype != "E2-R", ]
batch2_negative_susceptible_late <- batch2_negative_susceptible[
  batch2_negative_susceptible$hpi %in% c("72", "120"),
]
batch2_negative_susceptible_late <- batch2_negative_susceptible_late[, !(names(batch2_negative_susceptible_late) %in% c("FEATURE_ID", "genotype", "hpi", "Sample_ID", "Injection_order", "Group", "QC"))]
batch2_negative_susceptible_late <- t(batch2_negative_susceptible_late)



# 2. Extract genotype and data

genotype <- as.character(unlist(batch2_negative_susceptible_late[1, ])) # first row = genotype
data <- batch2_negative_susceptible_late[-1, ]                          # remaining rows = metabolites
data <- apply(data, 2, as.numeric)            # convert to numeric
rownames(data) <- rownames(batch2_negative_susceptible_late)[-1]

# 3. Compute log2 fold change

mean_R <- rowMeans(data[, genotype == "E2-AP"])
mean_S <- rowMeans(data[, genotype == "E2-C"])

log2FC <- log2(mean_R / mean_S)   


# 4. Compute p-values (t-test)

pvalues <- apply(data, 1, function(x) {
  t.test(x[genotype == "E2-AP"], x[genotype == "E2-C"])$p.value
})


# 5. Adjust p-values (FDR)

padj <- p.adjust(pvalues, method = "fdr")


# 6. Create results_batch2_negative_susceptible_late table

results_batch2_negative_susceptible_late <- data.frame(
  metabolite = rownames(data),
  log2FC = log2FC,
  pvalue = pvalues,
  padj = padj
)

results_batch2_negative_susceptible_late$negLog10P <- -log10(results_batch2_negative_susceptible_late$padj)


# 7. Volcano plot

# Define significance
results_batch2_negative_susceptible_late$significant <- "Not Significant"
results_batch2_negative_susceptible_late$significant[results_batch2_negative_susceptible_late$log2FC > 1 & results_batch2_negative_susceptible_late$padj < 0.05] <- "Up"
results_batch2_negative_susceptible_late$significant[results_batch2_negative_susceptible_late$log2FC < -1 & results_batch2_negative_susceptible_late$padj < 0.05] <- "Down"

Volcano_batch2_negative_susceptible_late <- ggplot(results_batch2_negative_susceptible_late, aes(x = log2FC, y = negLog10P, color = significant)) +
  geom_point(size = 3, alpha = 0.8) +
  scale_color_manual(values = c("blue", "grey", "red")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  theme_minimal() +
  labs(
    x = "Log2 Fold Change (E2-AP / E2-C)",
    y = "-Log10(FDR)",
    title = "Volcano Plot"
  ) +
  geom_text_repel(
    data = subset(results_batch2_negative_susceptible_late, significant != "Not Significant"),
    aes(label = metabolite)
  )+
  # Set fixed x and y axis limits
  scale_x_continuous(limits = c(-5, 5)) +  # Adjust min/max for log2FC
  scale_y_continuous(limits = c(0, 20))    # Adjust min/max for -log10(FDR)

results_batch2_negative_susceptible_late$significant<-as.factor(results_batch2_negative_susceptible_late$significant)
summary(results_batch2_negative_susceptible_late)




#3hpi----
# Keep only rows where genotype is NOT "E2-R"
batch2_negative_susceptible <- batch2_negative[batch2_negative$genotype != "E2-R", ]
batch2_negative_susceptible_3hpi <- batch2_negative_susceptible[batch2_negative_susceptible$hpi == "3", ]
batch2_negative_susceptible_3hpi <- batch2_negative_susceptible_3hpi[, !(names(batch2_negative_susceptible_3hpi) %in% c("FEATURE_ID", "genotype", "hpi", "Sample_ID", "Injection_order", "Group", "QC"))]
batch2_negative_susceptible_3hpi <- t(batch2_negative_susceptible_3hpi)



# 2. Extract genotype and data

genotype <- as.character(unlist(batch2_negative_susceptible_3hpi[1, ])) # first row = genotype
data <- batch2_negative_susceptible_3hpi[-1, ]                          # remaining rows = metabolites
data <- apply(data, 2, as.numeric)            # convert to numeric
rownames(data) <- rownames(batch2_negative_susceptible_3hpi)[-1]

# 3. Compute log2 fold change

mean_R <- rowMeans(data[, genotype == "E2-AP"])
mean_S <- rowMeans(data[, genotype == "E2-C"])

log2FC <- log2(mean_R / mean_S)   


# 4. Compute p-values (t-test)

pvalues <- apply(data, 1, function(x) {
  t.test(x[genotype == "E2-AP"], x[genotype == "E2-C"])$p.value
})


# 5. Adjust p-values (FDR)

padj <- p.adjust(pvalues, method = "fdr")


# 6. Create results_batch2_negative_susceptible_3hpi table

results_batch2_negative_susceptible_3hpi <- data.frame(
  metabolite = rownames(data),
  log2FC = log2FC,
  pvalue = pvalues,
  padj = padj
)

results_batch2_negative_susceptible_3hpi$negLog10P <- -log10(results_batch2_negative_susceptible_3hpi$padj)


# 7. Volcano plot

# Define significance
results_batch2_negative_susceptible_3hpi$significant <- "Not Significant"
results_batch2_negative_susceptible_3hpi$significant[results_batch2_negative_susceptible_3hpi$log2FC > 1 & results_batch2_negative_susceptible_3hpi$padj < 0.05] <- "Up"
results_batch2_negative_susceptible_3hpi$significant[results_batch2_negative_susceptible_3hpi$log2FC < -1 & results_batch2_negative_susceptible_3hpi$padj < 0.05] <- "Down"

Volcano_batch2_negative_susceptible_3hpi <- ggplot(results_batch2_negative_susceptible_3hpi, aes(x = log2FC, y = negLog10P, color = significant)) +
  geom_point(size = 3, alpha = 0.8) +
  scale_color_manual(values = c("blue", "grey", "red")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  theme_minimal() +
  labs(
    x = "Log2 Fold Change (E2-AP / E2-C)",
    y = "-Log10(FDR)",
    title = "Volcano Plot"
  ) +
  geom_text_repel(
    data = subset(results_batch2_negative_susceptible_3hpi, significant != "Not Significant"),
    aes(label = metabolite)
  )+
  # Set fixed x and y axis limits
  scale_x_continuous(limits = c(-5, 5)) +  # Adjust min/max for log2FC
  scale_y_continuous(limits = c(0, 20))    # Adjust min/max for -log10(FDR)

results_batch2_negative_susceptible_3hpi$significant<-as.factor(results_batch2_negative_susceptible_3hpi$significant)
summary(results_batch2_negative_susceptible_3hpi)




#24hpi----
# Keep only rows where genotype is NOT "E2-R"
batch2_negative_susceptible <- batch2_negative[batch2_negative$genotype != "E2-R", ]
batch2_negative_susceptible_24hpi <- batch2_negative_susceptible[batch2_negative_susceptible$hpi == "24", ]
batch2_negative_susceptible_24hpi <- batch2_negative_susceptible_24hpi[, !(names(batch2_negative_susceptible_24hpi) %in% c("FEATURE_ID", "genotype", "hpi", "Sample_ID", "Injection_order", "Group", "QC"))]
batch2_negative_susceptible_24hpi <- t(batch2_negative_susceptible_24hpi)



# 2. Extract genotype and data

genotype <- as.character(unlist(batch2_negative_susceptible_24hpi[1, ])) # first row = genotype
data <- batch2_negative_susceptible_24hpi[-1, ]                          # remaining rows = metabolites
data <- apply(data, 2, as.numeric)            # convert to numeric
rownames(data) <- rownames(batch2_negative_susceptible_24hpi)[-1]

# 3. Compute log2 fold change

mean_R <- rowMeans(data[, genotype == "E2-AP"])
mean_S <- rowMeans(data[, genotype == "E2-C"])

log2FC <- log2(mean_R / mean_S)   


# 4. Compute p-values (t-test)

pvalues <- apply(data, 1, function(x) {
  t.test(x[genotype == "E2-AP"], x[genotype == "E2-C"])$p.value
})


# 5. Adjust p-values (FDR)

padj <- p.adjust(pvalues, method = "fdr")


# 6. Create results_batch2_negative_susceptible_24hpi table

results_batch2_negative_susceptible_24hpi <- data.frame(
  metabolite = rownames(data),
  log2FC = log2FC,
  pvalue = pvalues,
  padj = padj
)

results_batch2_negative_susceptible_24hpi$negLog10P <- -log10(results_batch2_negative_susceptible_24hpi$padj)


# 7. Volcano plot

# Define significance
results_batch2_negative_susceptible_24hpi$significant <- "Not Significant"
results_batch2_negative_susceptible_24hpi$significant[results_batch2_negative_susceptible_24hpi$log2FC > 1 & results_batch2_negative_susceptible_24hpi$padj < 0.05] <- "Up"
results_batch2_negative_susceptible_24hpi$significant[results_batch2_negative_susceptible_24hpi$log2FC < -1 & results_batch2_negative_susceptible_24hpi$padj < 0.05] <- "Down"

Volcano_batch2_negative_susceptible_24hpi <- ggplot(results_batch2_negative_susceptible_24hpi, aes(x = log2FC, y = negLog10P, color = significant)) +
  geom_point(size = 3, alpha = 0.8) +
  scale_color_manual(values = c("blue", "grey", "red")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  theme_minimal() +
  labs(
    x = "Log2 Fold Change (E2-AP / E2-C)",
    y = "-Log10(FDR)",
    title = "Volcano Plot"
  ) +
  geom_text_repel(
    data = subset(results_batch2_negative_susceptible_24hpi, significant != "Not Significant"),
    aes(label = metabolite)
  )+
  # Set fixed x and y axis limits
  scale_x_continuous(limits = c(-5, 5)) +  # Adjust min/max for log2FC
  scale_y_continuous(limits = c(0, 20))    # Adjust min/max for -log10(FDR)

results_batch2_negative_susceptible_24hpi$significant<-as.factor(results_batch2_negative_susceptible_24hpi$significant)
summary(results_batch2_negative_susceptible_24hpi)




#72hpi----
# Keep only rows where genotype is NOT "E2-R"
batch2_negative_susceptible <- batch2_negative[batch2_negative$genotype != "E2-R", ]
batch2_negative_susceptible_72hpi <- batch2_negative_susceptible[batch2_negative_susceptible$hpi == "72", ]
batch2_negative_susceptible_72hpi <- batch2_negative_susceptible_72hpi[, !(names(batch2_negative_susceptible_72hpi) %in% c("FEATURE_ID", "genotype", "hpi", "Sample_ID", "Injection_order", "Group", "QC"))]
batch2_negative_susceptible_72hpi <- t(batch2_negative_susceptible_72hpi)



# 2. Extract genotype and data

genotype <- as.character(unlist(batch2_negative_susceptible_72hpi[1, ])) # first row = genotype
data <- batch2_negative_susceptible_72hpi[-1, ]                          # remaining rows = metabolites
data <- apply(data, 2, as.numeric)            # convert to numeric
rownames(data) <- rownames(batch2_negative_susceptible_72hpi)[-1]

# 3. Compute log2 fold change

mean_R <- rowMeans(data[, genotype == "E2-AP"])
mean_S <- rowMeans(data[, genotype == "E2-C"])

log2FC <- log2(mean_R / mean_S)   


# 4. Compute p-values (t-test)

pvalues <- apply(data, 1, function(x) {
  t.test(x[genotype == "E2-AP"], x[genotype == "E2-C"])$p.value
})


# 5. Adjust p-values (FDR)

padj <- p.adjust(pvalues, method = "fdr")


# 6. Create results_batch2_negative_susceptible_72hpi table

results_batch2_negative_susceptible_72hpi <- data.frame(
  metabolite = rownames(data),
  log2FC = log2FC,
  pvalue = pvalues,
  padj = padj
)

results_batch2_negative_susceptible_72hpi$negLog10P <- -log10(results_batch2_negative_susceptible_72hpi$padj)


# 7. Volcano plot

# Define significance
results_batch2_negative_susceptible_72hpi$significant <- "Not Significant"
results_batch2_negative_susceptible_72hpi$significant[results_batch2_negative_susceptible_72hpi$log2FC > 1 & results_batch2_negative_susceptible_72hpi$padj < 0.05] <- "Up"
results_batch2_negative_susceptible_72hpi$significant[results_batch2_negative_susceptible_72hpi$log2FC < -1 & results_batch2_negative_susceptible_72hpi$padj < 0.05] <- "Down"

Volcano_batch2_negative_susceptible_72hpi <- ggplot(results_batch2_negative_susceptible_72hpi, aes(x = log2FC, y = negLog10P, color = significant)) +
  geom_point(size = 3, alpha = 0.8) +
  scale_color_manual(values = c("blue", "grey", "red")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  theme_minimal() +
  labs(
    x = "Log2 Fold Change (E2-AP / E2-C)",
    y = "-Log10(FDR)",
    title = "Volcano Plot"
  ) +
  geom_text_repel(
    data = subset(results_batch2_negative_susceptible_72hpi, significant != "Not Significant"),
    aes(label = metabolite)
  )+
  # Set fixed x and y axis limits
  scale_x_continuous(limits = c(-5, 5)) +  # Adjust min/max for log2FC
  scale_y_continuous(limits = c(0, 20))    # Adjust min/max for -log10(FDR)

results_batch2_negative_susceptible_72hpi$significant<-as.factor(results_batch2_negative_susceptible_72hpi$significant)
summary(results_batch2_negative_susceptible_72hpi)

#120hpi----
# Keep only rows where genotype is NOT "E2-R"
batch2_negative_susceptible <- batch2_negative[batch2_negative$genotype != "E2-R", ]
batch2_negative_susceptible_120hpi <- batch2_negative_susceptible[batch2_negative_susceptible$hpi == "120", ]
batch2_negative_susceptible_120hpi <- batch2_negative_susceptible_120hpi[, !(names(batch2_negative_susceptible_120hpi) %in% c("FEATURE_ID", "genotype", "hpi", "Sample_ID", "Injection_order", "Group", "QC"))]
batch2_negative_susceptible_120hpi <- t(batch2_negative_susceptible_120hpi)



# 2. Extract genotype and data

genotype <- as.character(unlist(batch2_negative_susceptible_120hpi[1, ])) # first row = genotype
data <- batch2_negative_susceptible_120hpi[-1, ]                          # remaining rows = metabolites
data <- apply(data, 2, as.numeric)            # convert to numeric
rownames(data) <- rownames(batch2_negative_susceptible_120hpi)[-1]

# 3. Compute log2 fold change

mean_R <- rowMeans(data[, genotype == "E2-AP"])
mean_S <- rowMeans(data[, genotype == "E2-C"])

log2FC <- log2(mean_R / mean_S)   


# 4. Compute p-values (t-test)

pvalues <- apply(data, 1, function(x) {
  t.test(x[genotype == "E2-AP"], x[genotype == "E2-C"])$p.value
})


# 5. Adjust p-values (FDR)

padj <- p.adjust(pvalues, method = "fdr")


# 6. Create results_batch2_negative_susceptible_120hpi table

results_batch2_negative_susceptible_120hpi <- data.frame(
  metabolite = rownames(data),
  log2FC = log2FC,
  pvalue = pvalues,
  padj = padj
)

results_batch2_negative_susceptible_120hpi$negLog10P <- -log10(results_batch2_negative_susceptible_120hpi$padj)


# 7. Volcano plot

# Define significance
results_batch2_negative_susceptible_120hpi$significant <- "Not Significant"
results_batch2_negative_susceptible_120hpi$significant[results_batch2_negative_susceptible_120hpi$log2FC > 1 & results_batch2_negative_susceptible_120hpi$padj < 0.05] <- "Up"
results_batch2_negative_susceptible_120hpi$significant[results_batch2_negative_susceptible_120hpi$log2FC < -1 & results_batch2_negative_susceptible_120hpi$padj < 0.05] <- "Down"

Volcano_batch2_negative_susceptible_120hpi <- ggplot(results_batch2_negative_susceptible_120hpi, aes(x = log2FC, y = negLog10P, color = significant)) +
  geom_point(size = 3, alpha = 0.8) +
  scale_color_manual(values = c("blue", "grey", "red")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  theme_minimal() +
  labs(
    x = "Log2 Fold Change (E2-AP / E2-C)",
    y = "-Log10(FDR)",
    title = "Volcano Plot"
  ) +
  geom_text_repel(
    data = subset(results_batch2_negative_susceptible_120hpi, significant != "Not Significant"),
    aes(label = metabolite)
  )+
  # Set fixed x and y axis limits
  scale_x_continuous(limits = c(-5, 5)) +  # Adjust min/max for log2FC
  scale_y_continuous(limits = c(0, 20))    # Adjust min/max for -log10(FDR)

results_batch2_negative_susceptible_120hpi$significant<-as.factor(results_batch2_negative_susceptible_120hpi$significant)
summary(results_batch2_negative_susceptible_120hpi)







##batch2_negative_resistant----
#batch2_negative_resistant----
# Keep only rows where genotype is NOT "E2-S"
batch2_negative_resistant <- batch2_negative[batch2_negative$genotype != "E2-S", ]
batch2_negative_resistant <- batch2_negative_resistant[, !(names(batch2_negative_resistant) %in% c("FEATURE_ID", "genotype", "hpi", "Sample_ID", "Injection_order", "Group", "QC"))]
batch2_negative_resistant <- t(batch2_negative_resistant)



# 2. Extract genotype and data

genotype <- as.character(unlist(batch2_negative_resistant[1, ])) # first row = genotype
data <- batch2_negative_resistant[-1, ]                          # remaining rows = metabolites
data <- apply(data, 2, as.numeric)            # convert to numeric
rownames(data) <- rownames(batch2_negative_resistant)[-1]

# 3. Compute log2 fold change

mean_R <- rowMeans(data[, genotype == "E2-AP"])
mean_S <- rowMeans(data[, genotype == "E2-C"])

log2FC <- log2(mean_R / mean_S)   


# 4. Compute p-values (t-test)

pvalues <- apply(data, 1, function(x) {
  t.test(x[genotype == "E2-AP"], x[genotype == "E2-C"])$p.value
})


# 5. Adjust p-values (FDR)

padj <- p.adjust(pvalues, method = "fdr")


# 6. Create results_batch2_negative_resistant table

results_batch2_negative_resistant <- data.frame(
  metabolite = rownames(data),
  log2FC = log2FC,
  pvalue = pvalues,
  padj = padj
)

results_batch2_negative_resistant$negLog10P <- -log10(results_batch2_negative_resistant$padj)


# 7. Volcano plot

# Define significance
results_batch2_negative_resistant$significant <- "Not Significant"
results_batch2_negative_resistant$significant[results_batch2_negative_resistant$log2FC > 1 & results_batch2_negative_resistant$padj < 0.05] <- "Up"
results_batch2_negative_resistant$significant[results_batch2_negative_resistant$log2FC < -1 & results_batch2_negative_resistant$padj < 0.05] <- "Down"

Volcano_batch2_negative_resistant <- ggplot(results_batch2_negative_resistant, aes(x = log2FC, y = negLog10P, color = significant)) +
  geom_point(size = 3, alpha = 0.8) +
  scale_color_manual(values = c("grey", "red")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  theme_minimal() +
  labs(
    x = "Log2 Fold Change (E2-AP / E2-C)",
    y = "-Log10(FDR)",
    title = "Volcano Plot"
  ) +
  geom_text_repel(
    data = subset(results_batch2_negative_resistant, significant != "Not Significant"),
    aes(label = metabolite)
  )+
  # Set fixed x and y axis limits
  scale_x_continuous(limits = c(-5, 5)) +  # Adjust min/max for log2FC
  scale_y_continuous(limits = c(0, 5))    # Adjust min/max for -log10(FDR)

results_batch2_negative_resistant$significant<-as.factor(results_batch2_negative_resistant$significant)
summary(results_batch2_negative_resistant)




#early----
# Keep only rows where genotype is NOT "E2-S"
batch2_negative_resistant <- batch2_negative[batch2_negative$genotype != "E2-S", ]
batch2_negative_resistant_early <- batch2_negative_resistant[
  batch2_negative_resistant$hpi %in% c("3", "24"),]
batch2_negative_resistant_early <- batch2_negative_resistant_early[, !(names(batch2_negative_resistant_early) %in% c("FEATURE_ID", "genotype", "hpi", "Sample_ID", "Injection_order", "Group", "QC"))]
batch2_negative_resistant_early <- t(batch2_negative_resistant_early)



# 2. Extract genotype and data

genotype <- as.character(unlist(batch2_negative_resistant_early[1, ])) # first row = genotype
data <- batch2_negative_resistant_early[-1, ]                          # remaining rows = metabolites
data <- apply(data, 2, as.numeric)            # convert to numeric
rownames(data) <- rownames(batch2_negative_resistant_early)[-1]

# 3. Compute log2 fold change

mean_R <- rowMeans(data[, genotype == "E2-AP"])
mean_S <- rowMeans(data[, genotype == "E2-C"])

log2FC <- log2(mean_R / mean_S)   


# 4. Compute p-values (t-test)

pvalues <- apply(data, 1, function(x) {
  t.test(x[genotype == "E2-AP"], x[genotype == "E2-C"])$p.value
})


# 5. Adjust p-values (FDR)

padj <- p.adjust(pvalues, method = "fdr")


# 6. Create results_batch2_negative_resistant_early table

results_batch2_negative_resistant_early <- data.frame(
  metabolite = rownames(data),
  log2FC = log2FC,
  pvalue = pvalues,
  padj = padj
)

results_batch2_negative_resistant_early$negLog10P <- -log10(results_batch2_negative_resistant_early$padj)


# 7. Volcano plot

# Define significance
results_batch2_negative_resistant_early$significant <- "Not Significant"
results_batch2_negative_resistant_early$significant[results_batch2_negative_resistant_early$log2FC > 1 & results_batch2_negative_resistant_early$padj < 0.05] <- "Up"
results_batch2_negative_resistant_early$significant[results_batch2_negative_resistant_early$log2FC < -1 & results_batch2_negative_resistant_early$padj < 0.05] <- "Down"

Volcano_batch2_negative_resistant_early <- ggplot(results_batch2_negative_resistant_early, aes(x = log2FC, y = negLog10P, color = significant)) +
  geom_point(size = 3, alpha = 0.8) +
  scale_color_manual(values = c("blue", "grey", "red")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  theme_minimal() +
  labs(
    x = "Log2 Fold Change (E2-AP / E2-C)",
    y = "-Log10(FDR)",
    title = "Volcano Plot"
  ) +
  geom_text_repel(
    data = subset(results_batch2_negative_resistant_early, significant != "Not Significant"),
    aes(label = metabolite)
  )+
  # Set fixed x and y axis limits
  scale_x_continuous(limits = c(-5, 5)) +  # Adjust min/max for log2FC
  scale_y_continuous(limits = c(0, 20))    # Adjust min/max for -log10(FDR)

results_batch2_negative_resistant_early$significant<-as.factor(results_batch2_negative_resistant_early$significant)
summary(results_batch2_negative_resistant_early)





#late----
# Keep only rows where genotype is NOT "E2-S"
batch2_negative_resistant <- batch2_negative[batch2_negative$genotype != "E2-S", ]
batch2_negative_resistant_late <- batch2_negative_resistant[
  batch2_negative_resistant$hpi %in% c("72", "120"),
]
batch2_negative_resistant_late <- batch2_negative_resistant_late[, !(names(batch2_negative_resistant_late) %in% c("FEATURE_ID", "genotype", "hpi", "Sample_ID", "Injection_order", "Group", "QC"))]
batch2_negative_resistant_late <- t(batch2_negative_resistant_late)



# 2. Extract genotype and data

genotype <- as.character(unlist(batch2_negative_resistant_late[1, ])) # first row = genotype
data <- batch2_negative_resistant_late[-1, ]                          # remaining rows = metabolites
data <- apply(data, 2, as.numeric)            # convert to numeric
rownames(data) <- rownames(batch2_negative_resistant_late)[-1]

# 3. Compute log2 fold change

mean_R <- rowMeans(data[, genotype == "E2-AP"])
mean_S <- rowMeans(data[, genotype == "E2-C"])

log2FC <- log2(mean_R / mean_S)   


# 4. Compute p-values (t-test)

pvalues <- apply(data, 1, function(x) {
  t.test(x[genotype == "E2-AP"], x[genotype == "E2-C"])$p.value
})


# 5. Adjust p-values (FDR)

padj <- p.adjust(pvalues, method = "fdr")


# 6. Create results_batch2_negative_resistant_late table

results_batch2_negative_resistant_late <- data.frame(
  metabolite = rownames(data),
  log2FC = log2FC,
  pvalue = pvalues,
  padj = padj
)

results_batch2_negative_resistant_late$negLog10P <- -log10(results_batch2_negative_resistant_late$padj)


# 7. Volcano plot

# Define significance
results_batch2_negative_resistant_late$significant <- "Not Significant"
results_batch2_negative_resistant_late$significant[results_batch2_negative_resistant_late$log2FC > 1 & results_batch2_negative_resistant_late$padj < 0.05] <- "Up"
results_batch2_negative_resistant_late$significant[results_batch2_negative_resistant_late$log2FC < -1 & results_batch2_negative_resistant_late$padj < 0.05] <- "Down"

Volcano_batch2_negative_resistant_late <- ggplot(results_batch2_negative_resistant_late, aes(x = log2FC, y = negLog10P, color = significant)) +
  geom_point(size = 3, alpha = 0.8) +
  scale_color_manual(values = c("blue", "grey", "red")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  theme_minimal() +
  labs(
    x = "Log2 Fold Change (E2-AP / E2-C)",
    y = "-Log10(FDR)",
    title = "Volcano Plot"
  ) +
  geom_text_repel(
    data = subset(results_batch2_negative_resistant_late, significant != "Not Significant"),
    aes(label = metabolite)
  )+
  # Set fixed x and y axis limits
  scale_x_continuous(limits = c(-5, 5)) +  # Adjust min/max for log2FC
  scale_y_continuous(limits = c(0, 20))    # Adjust min/max for -log10(FDR)

results_batch2_negative_resistant_late$significant<-as.factor(results_batch2_negative_resistant_late$significant)
summary(results_batch2_negative_resistant_late)




#3hpi----
# Keep only rows where genotype is NOT "E2-S"
batch2_negative_resistant <- batch2_negative[batch2_negative$genotype != "E2-S", ]
batch2_negative_resistant_3hpi <- batch2_negative_resistant[batch2_negative_resistant$hpi == "3", ]
batch2_negative_resistant_3hpi <- batch2_negative_resistant_3hpi[, !(names(batch2_negative_resistant_3hpi) %in% c("FEATURE_ID", "genotype", "hpi", "Sample_ID", "Injection_order", "Group", "QC"))]
batch2_negative_resistant_3hpi <- t(batch2_negative_resistant_3hpi)



# 2. Extract genotype and data

genotype <- as.character(unlist(batch2_negative_resistant_3hpi[1, ])) # first row = genotype
data <- batch2_negative_resistant_3hpi[-1, ]                          # remaining rows = metabolites
data <- apply(data, 2, as.numeric)            # convert to numeric
rownames(data) <- rownames(batch2_negative_resistant_3hpi)[-1]

# 3. Compute log2 fold change

mean_R <- rowMeans(data[, genotype == "E2-AP"])
mean_S <- rowMeans(data[, genotype == "E2-C"])

log2FC <- log2(mean_R / mean_S)   


# 4. Compute p-values (t-test)

pvalues <- apply(data, 1, function(x) {
  t.test(x[genotype == "E2-AP"], x[genotype == "E2-C"])$p.value
})


# 5. Adjust p-values (FDR)

padj <- p.adjust(pvalues, method = "fdr")


# 6. Create results_batch2_negative_resistant_3hpi table

results_batch2_negative_resistant_3hpi <- data.frame(
  metabolite = rownames(data),
  log2FC = log2FC,
  pvalue = pvalues,
  padj = padj
)

results_batch2_negative_resistant_3hpi$negLog10P <- -log10(results_batch2_negative_resistant_3hpi$padj)


# 7. Volcano plot

# Define significance
results_batch2_negative_resistant_3hpi$significant <- "Not Significant"
results_batch2_negative_resistant_3hpi$significant[results_batch2_negative_resistant_3hpi$log2FC > 1 & results_batch2_negative_resistant_3hpi$padj < 0.05] <- "Up"
results_batch2_negative_resistant_3hpi$significant[results_batch2_negative_resistant_3hpi$log2FC < -1 & results_batch2_negative_resistant_3hpi$padj < 0.05] <- "Down"

Volcano_batch2_negative_resistant_3hpi <- ggplot(results_batch2_negative_resistant_3hpi, aes(x = log2FC, y = negLog10P, color = significant)) +
  geom_point(size = 3, alpha = 0.8) +
  scale_color_manual(values = c("blue", "grey", "red")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  theme_minimal() +
  labs(
    x = "Log2 Fold Change (E2-AP / E2-C)",
    y = "-Log10(FDR)",
    title = "Volcano Plot"
  ) +
  geom_text_repel(
    data = subset(results_batch2_negative_resistant_3hpi, significant != "Not Significant"),
    aes(label = metabolite)
  )+
  # Set fixed x and y axis limits
  scale_x_continuous(limits = c(-5, 5)) +  # Adjust min/max for log2FC
  scale_y_continuous(limits = c(0, 20))    # Adjust min/max for -log10(FDR)

results_batch2_negative_resistant_3hpi$significant<-as.factor(results_batch2_negative_resistant_3hpi$significant)
summary(results_batch2_negative_resistant_3hpi)




#24hpi----
# Keep only rows where genotype is NOT "E2-S"
batch2_negative_resistant <- batch2_negative[batch2_negative$genotype != "E2-S", ]
batch2_negative_resistant_24hpi <- batch2_negative_resistant[batch2_negative_resistant$hpi == "24", ]
batch2_negative_resistant_24hpi <- batch2_negative_resistant_24hpi[, !(names(batch2_negative_resistant_24hpi) %in% c("FEATURE_ID", "genotype", "hpi", "Sample_ID", "Injection_order", "Group", "QC"))]
batch2_negative_resistant_24hpi <- t(batch2_negative_resistant_24hpi)



# 2. Extract genotype and data

genotype <- as.character(unlist(batch2_negative_resistant_24hpi[1, ])) # first row = genotype
data <- batch2_negative_resistant_24hpi[-1, ]                          # remaining rows = metabolites
data <- apply(data, 2, as.numeric)            # convert to numeric
rownames(data) <- rownames(batch2_negative_resistant_24hpi)[-1]

# 3. Compute log2 fold change

mean_R <- rowMeans(data[, genotype == "E2-AP"])
mean_S <- rowMeans(data[, genotype == "E2-C"])

log2FC <- log2(mean_R / mean_S)   


# 4. Compute p-values (t-test)

pvalues <- apply(data, 1, function(x) {
  t.test(x[genotype == "E2-AP"], x[genotype == "E2-C"])$p.value
})


# 5. Adjust p-values (FDR)

padj <- p.adjust(pvalues, method = "fdr")


# 6. Create results_batch2_negative_resistant_24hpi table

results_batch2_negative_resistant_24hpi <- data.frame(
  metabolite = rownames(data),
  log2FC = log2FC,
  pvalue = pvalues,
  padj = padj
)

results_batch2_negative_resistant_24hpi$negLog10P <- -log10(results_batch2_negative_resistant_24hpi$padj)


# 7. Volcano plot

# Define significance
results_batch2_negative_resistant_24hpi$significant <- "Not Significant"
results_batch2_negative_resistant_24hpi$significant[results_batch2_negative_resistant_24hpi$log2FC > 1 & results_batch2_negative_resistant_24hpi$padj < 0.05] <- "Up"
results_batch2_negative_resistant_24hpi$significant[results_batch2_negative_resistant_24hpi$log2FC < -1 & results_batch2_negative_resistant_24hpi$padj < 0.05] <- "Down"

Volcano_batch2_negative_resistant_24hpi <- ggplot(results_batch2_negative_resistant_24hpi, aes(x = log2FC, y = negLog10P, color = significant)) +
  geom_point(size = 3, alpha = 0.8) +
  scale_color_manual(values = c("blue", "grey", "red")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  theme_minimal() +
  labs(
    x = "Log2 Fold Change (E2-AP / E2-C)",
    y = "-Log10(FDR)",
    title = "Volcano Plot"
  ) +
  geom_text_repel(
    data = subset(results_batch2_negative_resistant_24hpi, significant != "Not Significant"),
    aes(label = metabolite)
  )+
  # Set fixed x and y axis limits
  scale_x_continuous(limits = c(-5, 5)) +  # Adjust min/max for log2FC
  scale_y_continuous(limits = c(0, 20))    # Adjust min/max for -log10(FDR)

results_batch2_negative_resistant_24hpi$significant<-as.factor(results_batch2_negative_resistant_24hpi$significant)
summary(results_batch2_negative_resistant_24hpi)




#72hpi----
# Keep only rows where genotype is NOT "E2-S"
batch2_negative_resistant <- batch2_negative[batch2_negative$genotype != "E2-S", ]
batch2_negative_resistant_72hpi <- batch2_negative_resistant[batch2_negative_resistant$hpi == "72", ]
batch2_negative_resistant_72hpi <- batch2_negative_resistant_72hpi[, !(names(batch2_negative_resistant_72hpi) %in% c("FEATURE_ID", "genotype", "hpi", "Sample_ID", "Injection_order", "Group", "QC"))]
batch2_negative_resistant_72hpi <- t(batch2_negative_resistant_72hpi)



# 2. Extract genotype and data

genotype <- as.character(unlist(batch2_negative_resistant_72hpi[1, ])) # first row = genotype
data <- batch2_negative_resistant_72hpi[-1, ]                          # remaining rows = metabolites
data <- apply(data, 2, as.numeric)            # convert to numeric
rownames(data) <- rownames(batch2_negative_resistant_72hpi)[-1]

# 3. Compute log2 fold change

mean_R <- rowMeans(data[, genotype == "E2-AP"])
mean_S <- rowMeans(data[, genotype == "E2-C"])

log2FC <- log2(mean_R / mean_S)   


# 4. Compute p-values (t-test)

pvalues <- apply(data, 1, function(x) {
  t.test(x[genotype == "E2-AP"], x[genotype == "E2-C"])$p.value
})


# 5. Adjust p-values (FDR)

padj <- p.adjust(pvalues, method = "fdr")


# 6. Create results_batch2_negative_resistant_72hpi table

results_batch2_negative_resistant_72hpi <- data.frame(
  metabolite = rownames(data),
  log2FC = log2FC,
  pvalue = pvalues,
  padj = padj
)

results_batch2_negative_resistant_72hpi$negLog10P <- -log10(results_batch2_negative_resistant_72hpi$padj)


# 7. Volcano plot

# Define significance
results_batch2_negative_resistant_72hpi$significant <- "Not Significant"
results_batch2_negative_resistant_72hpi$significant[results_batch2_negative_resistant_72hpi$log2FC > 1 & results_batch2_negative_resistant_72hpi$padj < 0.05] <- "Up"
results_batch2_negative_resistant_72hpi$significant[results_batch2_negative_resistant_72hpi$log2FC < -1 & results_batch2_negative_resistant_72hpi$padj < 0.05] <- "Down"

Volcano_batch2_negative_resistant_72hpi <- ggplot(results_batch2_negative_resistant_72hpi, aes(x = log2FC, y = negLog10P, color = significant)) +
  geom_point(size = 3, alpha = 0.8) +
  scale_color_manual(values = c("blue", "grey", "red")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  theme_minimal() +
  labs(
    x = "Log2 Fold Change (E2-AP / E2-C)",
    y = "-Log10(FDR)",
    title = "Volcano Plot"
  ) +
  geom_text_repel(
    data = subset(results_batch2_negative_resistant_72hpi, significant != "Not Significant"),
    aes(label = metabolite)
  )+
  # Set fixed x and y axis limits
  scale_x_continuous(limits = c(-5, 5)) +  # Adjust min/max for log2FC
  scale_y_continuous(limits = c(0, 20))    # Adjust min/max for -log10(FDR)

results_batch2_negative_resistant_72hpi$significant<-as.factor(results_batch2_negative_resistant_72hpi$significant)
summary(results_batch2_negative_resistant_72hpi)

#120hpi----
# Keep only rows where genotype is NOT "E2-S"
batch2_negative_resistant <- batch2_negative[batch2_negative$genotype != "E2-S", ]
batch2_negative_resistant_120hpi <- batch2_negative_resistant[batch2_negative_resistant$hpi == "120", ]
batch2_negative_resistant_120hpi <- batch2_negative_resistant_120hpi[, !(names(batch2_negative_resistant_120hpi) %in% c("FEATURE_ID", "genotype", "hpi", "Sample_ID", "Injection_order", "Group", "QC"))]
batch2_negative_resistant_120hpi <- t(batch2_negative_resistant_120hpi)



# 2. Extract genotype and data

genotype <- as.character(unlist(batch2_negative_resistant_120hpi[1, ])) # first row = genotype
data <- batch2_negative_resistant_120hpi[-1, ]                          # remaining rows = metabolites
data <- apply(data, 2, as.numeric)            # convert to numeric
rownames(data) <- rownames(batch2_negative_resistant_120hpi)[-1]

# 3. Compute log2 fold change

mean_R <- rowMeans(data[, genotype == "E2-AP"])
mean_S <- rowMeans(data[, genotype == "E2-C"])

log2FC <- log2(mean_R / mean_S)   


# 4. Compute p-values (t-test)

pvalues <- apply(data, 1, function(x) {
  t.test(x[genotype == "E2-AP"], x[genotype == "E2-C"])$p.value
})


# 5. Adjust p-values (FDR)

padj <- p.adjust(pvalues, method = "fdr")


# 6. Create results_batch2_negative_resistant_120hpi table

results_batch2_negative_resistant_120hpi <- data.frame(
  metabolite = rownames(data),
  log2FC = log2FC,
  pvalue = pvalues,
  padj = padj
)

results_batch2_negative_resistant_120hpi$negLog10P <- -log10(results_batch2_negative_resistant_120hpi$padj)


# 7. Volcano plot

# Define significance
results_batch2_negative_resistant_120hpi$significant <- "Not Significant"
results_batch2_negative_resistant_120hpi$significant[results_batch2_negative_resistant_120hpi$log2FC > 1 & results_batch2_negative_resistant_120hpi$padj < 0.05] <- "Up"
results_batch2_negative_resistant_120hpi$significant[results_batch2_negative_resistant_120hpi$log2FC < -1 & results_batch2_negative_resistant_120hpi$padj < 0.05] <- "Down"

Volcano_batch2_negative_resistant_120hpi <- ggplot(results_batch2_negative_resistant_120hpi, aes(x = log2FC, y = negLog10P, color = significant)) +
  geom_point(size = 3, alpha = 0.8) +
  scale_color_manual(values = c("blue", "grey", "red")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  theme_minimal() +
  labs(
    x = "Log2 Fold Change (E2-AP / E2-C)",
    y = "-Log10(FDR)",
    title = "Volcano Plot"
  ) +
  geom_text_repel(
    data = subset(results_batch2_negative_resistant_120hpi, significant != "Not Significant"),
    aes(label = metabolite)
  )+
  # Set fixed x and y axis limits
  scale_x_continuous(limits = c(-5, 5)) +  # Adjust min/max for log2FC
  scale_y_continuous(limits = c(0, 20))    # Adjust min/max for -log10(FDR)

results_batch2_negative_resistant_120hpi$significant<-as.factor(results_batch2_negative_resistant_120hpi$significant)
summary(results_batch2_negative_resistant_120hpi)





#####

results_batch2_negative_susceptible<-results_batch2_negative_susceptible[,c(1,6)]
colnames(results_batch2_negative_susceptible)[colnames(results_batch2_negative_susceptible) == "significant"] <- "batch2_negative_susceptible"
results_batch2_negative_resistant<-results_batch2_negative_resistant[,c(1,6)]
colnames(results_batch2_negative_resistant)[colnames(results_batch2_negative_resistant) == "significant"] <- "batch2_negative_resistant"

results_batch2_negative_susceptible_early<-results_batch2_negative_susceptible_early[,c(1,6)]
colnames(results_batch2_negative_susceptible_early)[colnames(results_batch2_negative_susceptible_early) == "significant"] <- "batch2_negative_susceptible_early"

results_batch2_negative_resistant_early<-results_batch2_negative_resistant_early[,c(1,6)]
colnames(results_batch2_negative_resistant_early)[colnames(results_batch2_negative_resistant_early) == "significant"] <- "batch2_negative_resistant_early"


results_batch2_negative_susceptible_late<-results_batch2_negative_susceptible_late[,c(1,6)]
colnames(results_batch2_negative_susceptible_late)[colnames(results_batch2_negative_susceptible_late) == "significant"] <- "batch2_negative_susceptible_late"
results_batch2_negative_resistant_late<-results_batch2_negative_resistant_late[,c(1,6)]
colnames(results_batch2_negative_resistant_late)[colnames(results_batch2_negative_resistant_late) == "significant"] <- "batch2_negative_resistant_late"
results_batch2_negative_susceptible_3hpi<-results_batch2_negative_susceptible_3hpi[,c(1,6)]
colnames(results_batch2_negative_susceptible_3hpi)[colnames(results_batch2_negative_susceptible_3hpi) == "significant"] <- "batch2_negative_susceptible_3hpi"
results_batch2_negative_resistant_3hpi<-results_batch2_negative_resistant_3hpi[,c(1,6)]
colnames(results_batch2_negative_resistant_3hpi)[colnames(results_batch2_negative_resistant_3hpi) == "significant"] <- "batch2_negative_resistant_3hpi"
results_batch2_negative_susceptible_24hpi<-results_batch2_negative_susceptible_24hpi[,c(1,6)]
colnames(results_batch2_negative_susceptible_24hpi)[colnames(results_batch2_negative_susceptible_24hpi) == "significant"] <- "batch2_negative_susceptible_24hpi"
results_batch2_negative_resistant_24hpi<-results_batch2_negative_resistant_24hpi[,c(1,6)]
colnames(results_batch2_negative_resistant_24hpi)[colnames(results_batch2_negative_resistant_24hpi) == "significant"] <- "batch2_negative_resistant_24hpi"
results_batch2_negative_susceptible_72hpi<-results_batch2_negative_susceptible_72hpi[,c(1,6)]
colnames(results_batch2_negative_susceptible_72hpi)[colnames(results_batch2_negative_susceptible_72hpi) == "significant"] <- "batch2_negative_susceptible_72hpi"
results_batch2_negative_resistant_72hpi<-results_batch2_negative_resistant_72hpi[,c(1,6)]
colnames(results_batch2_negative_resistant_72hpi)[colnames(results_batch2_negative_resistant_72hpi) == "significant"] <- "batch2_negative_resistant_72hpi"
results_batch2_negative_susceptible_120hpi<-results_batch2_negative_susceptible_120hpi[,c(1,6)]
colnames(results_batch2_negative_susceptible_120hpi)[colnames(results_batch2_negative_susceptible_120hpi) == "significant"] <- "batch2_negative_susceptible_120hpi"
results_batch2_negative_resistant_120hpi<-results_batch2_negative_resistant_120hpi[,c(1,6)]
colnames(results_batch2_negative_resistant_120hpi)[colnames(results_batch2_negative_resistant_120hpi) == "significant"] <- "batch2_negative_resistant_120hpi"




all_results <- list(
  results_batch2_negative_susceptible,
  results_batch2_negative_resistant,
  results_batch2_negative_susceptible_early,
  results_batch2_negative_resistant_early,
  results_batch2_negative_susceptible_late,
  results_batch2_negative_resistant_late,
  results_batch2_negative_susceptible_3hpi,
  results_batch2_negative_resistant_3hpi,
  results_batch2_negative_susceptible_24hpi,
  results_batch2_negative_resistant_24hpi,
  results_batch2_negative_susceptible_72hpi,
  results_batch2_negative_resistant_72hpi,
  results_batch2_negative_susceptible_120hpi,
  results_batch2_negative_resistant_120hpi
)

merged_results_results_batch2_negative <- reduce(all_results, full_join, by = "metabolite")
merged_results_results_batch2_negative <- merged_results_results_batch2_negative[
  rowSums(merged_results_results_batch2_negative[, -1] != "Not Significant", na.rm = TRUE) > 0,
]

merged_results_results_batch2_negative$metabolite <- str_c(
  merged_results_results_batch2_negative$metabolite,
  "neg"
)

heatmap_batch2_negative <- merged_results_results_batch2_negative



heatmap_batch2_negative <- as.data.frame(lapply(heatmap_batch2_negative, as.character))
rownames(heatmap_batch2_negative) <- heatmap_batch2_negative$metabolite
heatmap_batch2_negative <- heatmap_batch2_negative[, -1]

heatmap_batch2_negative[] <- lapply(heatmap_batch2_negative, function(x) {
  x[x == "Up"] <- 1
  x[x == "Down"] <- -1
  x[x == "Not Significant"] <- 0
  return(as.numeric(x))
})

# Convert to matrix
heatmap_matrix <- as.matrix(heatmap_batch2_negative)

# Convert to numeric
mode(heatmap_matrix) <- "numeric"

# compute row sums (total signal per metabolite)
row_score <- rowSums(heatmap_matrix, na.rm = TRUE)

# order from highest to lowest
heatmap_matrix_batch2_negative <- heatmap_matrix[order(row_score, decreasing = TRUE), ]


# Plot heatmap
pheatmap(
  heatmap_matrix_batch2_negative,
  color = c("blue", "white", "red"),
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  border_color = "grey80"
)

setwd("~/Documents/metabolomic_too_heavy_for_git/untargeted/Positive mode")
###batch_1_positive----

batch1_positive <- read_excel("final_filtered_files/batch_1_positive_Intensity_filtered_notame.xlsx", col_names = FALSE)

batch1_positive <- batch1_positive %>%
  mutate(across(34, ~ ifelse(!is.na(batch1_positive[[3]]) & batch1_positive[[3]] != "", batch1_positive[[3]], .)))
batch1_positive <- batch1_positive[,-c(1:33)]
batch1_positive <- t(batch1_positive)
batch1_positive <- as.data.frame(batch1_positive)
colnames(batch1_positive) <- batch1_positive[1, ]
batch1_positive <- batch1_positive[-1, ]
batch1_positive_row_name <- batch1_positive %>%
  filter(Group != "Blank") %>%
  separate(Group, into = c("Batch","pheno","treatmentA","hpi"), sep = "-") %>%
  mutate(genotype = paste(Batch, pheno, sep = "-")) %>%
  mutate(treatment = paste(Batch, treatmentA, sep = "-")) %>%
  select(FEATURE_ID,genotype,treatment, hpi)
batch1_positive <- merge(batch1_positive_row_name,batch1_positive, by ="FEATURE_ID")
batch1_positive[is.na(batch1_positive)] <- 0 




##batch1_positive_susceptible----
#batch1_positive_susceptible----
# Keep only rows where genotype is NOT "E1-R"
batch1_positive_susceptible <- batch1_positive[batch1_positive$genotype != "E1-R", ]
batch1_positive_susceptible <- batch1_positive_susceptible[, !(names(batch1_positive_susceptible) %in% c("FEATURE_ID", "genotype", "hpi", "Sample_ID", "Injection_order", "Group", "QC"))]
batch1_positive_susceptible <- t(batch1_positive_susceptible)



# 2. Extract genotype and data

genotype <- as.character(unlist(batch1_positive_susceptible[1, ])) # first row = genotype
data <- batch1_positive_susceptible[-1, ]                          # remaining rows = metabolites
data <- apply(data, 2, as.numeric)            # convert to numeric
rownames(data) <- rownames(batch1_positive_susceptible)[-1]

# 3. Compute log2 fold change

mean_R <- rowMeans(data[, genotype == "E1-AF"])
mean_S <- rowMeans(data[, genotype == "E1-C"])

log2FC <- log2(mean_R / mean_S)   


# 4. Compute p-values (t-test)

pvalues <- apply(data, 1, function(x) {
  t.test(x[genotype == "E1-AF"], x[genotype == "E1-C"])$p.value
})


# 5. Adjust p-values (FDR)

padj <- p.adjust(pvalues, method = "fdr")


# 6. Create results_batch1_positive_susceptible table

results_batch1_positive_susceptible <- data.frame(
  metabolite = rownames(data),
  log2FC = log2FC,
  pvalue = pvalues,
  padj = padj
)

results_batch1_positive_susceptible$negLog10P <- -log10(results_batch1_positive_susceptible$padj)


# 7. Volcano plot

# Define significance
results_batch1_positive_susceptible$significant <- "Not Significant"
results_batch1_positive_susceptible$significant[results_batch1_positive_susceptible$log2FC > 1 & results_batch1_positive_susceptible$padj < 0.05] <- "Up"
results_batch1_positive_susceptible$significant[results_batch1_positive_susceptible$log2FC < -1 & results_batch1_positive_susceptible$padj < 0.05] <- "Down"

Volcano_batch1_positive_susceptible <- ggplot(results_batch1_positive_susceptible, aes(x = log2FC, y = negLog10P, color = significant)) +
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
    data = subset(results_batch1_positive_susceptible, significant != "Not Significant"),
    aes(label = metabolite)
  )+
  # Set fixed x and y axis limits
  scale_x_continuous(limits = c(-5, 5)) +  # Adjust min/max for log2FC
  scale_y_continuous(limits = c(0, 20))    # Adjust min/max for -log10(FDR)

results_batch1_positive_susceptible$significant<-as.factor(results_batch1_positive_susceptible$significant)
summary(results_batch1_positive_susceptible)




#early----
# Keep only rows where genotype is NOT "E1-R"
batch1_positive_susceptible <- batch1_positive[batch1_positive$genotype != "E1-R", ]
batch1_positive_susceptible_early <- batch1_positive_susceptible[
  batch1_positive_susceptible$hpi %in% c("3", "24"),]
batch1_positive_susceptible_early <- batch1_positive_susceptible_early[, !(names(batch1_positive_susceptible_early) %in% c("FEATURE_ID", "genotype", "hpi", "Sample_ID", "Injection_order", "Group", "QC"))]
batch1_positive_susceptible_early <- t(batch1_positive_susceptible_early)


# 2. Extract genotype and data

genotype <- as.character(unlist(batch1_positive_susceptible_early[1, ])) # first row = genotype
data <- batch1_positive_susceptible_early[-1, ]                          # remaining rows = metabolites
data <- apply(data, 2, as.numeric)            # convert to numeric
rownames(data) <- rownames(batch1_positive_susceptible_early)[-1]

# 3. Compute log2 fold change

mean_R <- rowMeans(data[, genotype == "E1-AF"])
mean_S <- rowMeans(data[, genotype == "E1-C"])

log2FC <- log2(mean_R / mean_S)   


# 4. Compute p-values (t-test)

pvalues <- apply(data, 1, function(x) {
  t.test(x[genotype == "E1-AF"], x[genotype == "E1-C"])$p.value
})


# 5. Adjust p-values (FDR)

padj <- p.adjust(pvalues, method = "fdr")


# 6. Create results_batch1_positive_susceptible_early table

results_batch1_positive_susceptible_early <- data.frame(
  metabolite = rownames(data),
  log2FC = log2FC,
  pvalue = pvalues,
  padj = padj
)

results_batch1_positive_susceptible_early$negLog10P <- -log10(results_batch1_positive_susceptible_early$padj)


# 7. Volcano plot

# Define significance
results_batch1_positive_susceptible_early$significant <- "Not Significant"
results_batch1_positive_susceptible_early$significant[results_batch1_positive_susceptible_early$log2FC > 1 & results_batch1_positive_susceptible_early$padj < 0.05] <- "Up"
results_batch1_positive_susceptible_early$significant[results_batch1_positive_susceptible_early$log2FC < -1 & results_batch1_positive_susceptible_early$padj < 0.05] <- "Down"

Volcano_batch1_positive_susceptible_early <- ggplot(results_batch1_positive_susceptible_early, aes(x = log2FC, y = negLog10P, color = significant)) +
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
    data = subset(results_batch1_positive_susceptible_early, significant != "Not Significant"),
    aes(label = metabolite)
  )+
  # Set fixed x and y axis limits
  scale_x_continuous(limits = c(-5, 5)) +  # Adjust min/max for log2FC
  scale_y_continuous(limits = c(0, 20))    # Adjust min/max for -log10(FDR)

results_batch1_positive_susceptible_early$significant<-as.factor(results_batch1_positive_susceptible_early$significant)
summary(results_batch1_positive_susceptible_early)





#late----
# Keep only rows where genotype is NOT "E1-R"
batch1_positive_susceptible <- batch1_positive[batch1_positive$genotype != "E1-R", ]
batch1_positive_susceptible_late <- batch1_positive_susceptible[
  batch1_positive_susceptible$hpi %in% c("72", "120"),
]
batch1_positive_susceptible_late <- batch1_positive_susceptible_late[, !(names(batch1_positive_susceptible_late) %in% c("FEATURE_ID", "genotype", "hpi", "Sample_ID", "Injection_order", "Group", "QC"))]
batch1_positive_susceptible_late <- t(batch1_positive_susceptible_late)



# 2. Extract genotype and data

genotype <- as.character(unlist(batch1_positive_susceptible_late[1, ])) # first row = genotype
data <- batch1_positive_susceptible_late[-1, ]                          # remaining rows = metabolites
data <- apply(data, 2, as.numeric)            # convert to numeric
rownames(data) <- rownames(batch1_positive_susceptible_late)[-1]

# 3. Compute log2 fold change

mean_R <- rowMeans(data[, genotype == "E1-AF"])
mean_S <- rowMeans(data[, genotype == "E1-C"])

log2FC <- log2(mean_R / mean_S)   


# 4. Compute p-values (t-test)

pvalues <- apply(data, 1, function(x) {
  t.test(x[genotype == "E1-AF"], x[genotype == "E1-C"])$p.value
})


# 5. Adjust p-values (FDR)

padj <- p.adjust(pvalues, method = "fdr")


# 6. Create results_batch1_positive_susceptible_late table

results_batch1_positive_susceptible_late <- data.frame(
  metabolite = rownames(data),
  log2FC = log2FC,
  pvalue = pvalues,
  padj = padj
)

results_batch1_positive_susceptible_late$negLog10P <- -log10(results_batch1_positive_susceptible_late$padj)


# 7. Volcano plot

# Define significance
results_batch1_positive_susceptible_late$significant <- "Not Significant"
results_batch1_positive_susceptible_late$significant[results_batch1_positive_susceptible_late$log2FC > 1 & results_batch1_positive_susceptible_late$padj < 0.05] <- "Up"
results_batch1_positive_susceptible_late$significant[results_batch1_positive_susceptible_late$log2FC < -1 & results_batch1_positive_susceptible_late$padj < 0.05] <- "Down"

Volcano_batch1_positive_susceptible_late <- ggplot(results_batch1_positive_susceptible_late, aes(x = log2FC, y = negLog10P, color = significant)) +
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
    data = subset(results_batch1_positive_susceptible_late, significant != "Not Significant"),
    aes(label = metabolite)
  )+
  # Set fixed x and y axis limits
  scale_x_continuous(limits = c(-5, 5)) +  # Adjust min/max for log2FC
  scale_y_continuous(limits = c(0, 20))    # Adjust min/max for -log10(FDR)

results_batch1_positive_susceptible_late$significant<-as.factor(results_batch1_positive_susceptible_late$significant)
summary(results_batch1_positive_susceptible_late)




#3hpi----
# Keep only rows where genotype is NOT "E1-R"
batch1_positive_susceptible <- batch1_positive[batch1_positive$genotype != "E1-R", ]
batch1_positive_susceptible_3hpi <- batch1_positive_susceptible[batch1_positive_susceptible$hpi == "3", ]
batch1_positive_susceptible_3hpi <- batch1_positive_susceptible_3hpi[, !(names(batch1_positive_susceptible_3hpi) %in% c("FEATURE_ID", "genotype", "hpi", "Sample_ID", "Injection_order", "Group", "QC"))]
batch1_positive_susceptible_3hpi <- t(batch1_positive_susceptible_3hpi)



# 2. Extract genotype and data

genotype <- as.character(unlist(batch1_positive_susceptible_3hpi[1, ])) # first row = genotype
data <- batch1_positive_susceptible_3hpi[-1, ]                          # remaining rows = metabolites
data <- apply(data, 2, as.numeric)            # convert to numeric
rownames(data) <- rownames(batch1_positive_susceptible_3hpi)[-1]

# 3. Compute log2 fold change

mean_R <- rowMeans(data[, genotype == "E1-AF"])
mean_S <- rowMeans(data[, genotype == "E1-C"])

log2FC <- log2(mean_R / mean_S)   


# 4. Compute p-values (t-test)

pvalues <- apply(data, 1, function(x) {
  t.test(x[genotype == "E1-AF"], x[genotype == "E1-C"])$p.value
})


# 5. Adjust p-values (FDR)

padj <- p.adjust(pvalues, method = "fdr")


# 6. Create results_batch1_positive_susceptible_3hpi table

results_batch1_positive_susceptible_3hpi <- data.frame(
  metabolite = rownames(data),
  log2FC = log2FC,
  pvalue = pvalues,
  padj = padj
)

results_batch1_positive_susceptible_3hpi$negLog10P <- -log10(results_batch1_positive_susceptible_3hpi$padj)


# 7. Volcano plot

# Define significance
results_batch1_positive_susceptible_3hpi$significant <- "Not Significant"
results_batch1_positive_susceptible_3hpi$significant[results_batch1_positive_susceptible_3hpi$log2FC > 1 & results_batch1_positive_susceptible_3hpi$padj < 0.05] <- "Up"
results_batch1_positive_susceptible_3hpi$significant[results_batch1_positive_susceptible_3hpi$log2FC < -1 & results_batch1_positive_susceptible_3hpi$padj < 0.05] <- "Down"

Volcano_batch1_positive_susceptible_3hpi <- ggplot(results_batch1_positive_susceptible_3hpi, aes(x = log2FC, y = negLog10P, color = significant)) +
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
    data = subset(results_batch1_positive_susceptible_3hpi, significant != "Not Significant"),
    aes(label = metabolite)
  )+
  # Set fixed x and y axis limits
  scale_x_continuous(limits = c(-5, 5)) +  # Adjust min/max for log2FC
  scale_y_continuous(limits = c(0, 20))    # Adjust min/max for -log10(FDR)

results_batch1_positive_susceptible_3hpi$significant<-as.factor(results_batch1_positive_susceptible_3hpi$significant)
summary(results_batch1_positive_susceptible_3hpi)




#24hpi----
# Keep only rows where genotype is NOT "E1-R"
batch1_positive_susceptible <- batch1_positive[batch1_positive$genotype != "E1-R", ]
batch1_positive_susceptible_24hpi <- batch1_positive_susceptible[batch1_positive_susceptible$hpi == "24", ]
batch1_positive_susceptible_24hpi <- batch1_positive_susceptible_24hpi[, !(names(batch1_positive_susceptible_24hpi) %in% c("FEATURE_ID", "genotype", "hpi", "Sample_ID", "Injection_order", "Group", "QC"))]
batch1_positive_susceptible_24hpi <- t(batch1_positive_susceptible_24hpi)



# 2. Extract genotype and data

genotype <- as.character(unlist(batch1_positive_susceptible_24hpi[1, ])) # first row = genotype
data <- batch1_positive_susceptible_24hpi[-1, ]                          # remaining rows = metabolites
data <- apply(data, 2, as.numeric)            # convert to numeric
rownames(data) <- rownames(batch1_positive_susceptible_24hpi)[-1]

# 3. Compute log2 fold change

mean_R <- rowMeans(data[, genotype == "E1-AF"])
mean_S <- rowMeans(data[, genotype == "E1-C"])

log2FC <- log2(mean_R / mean_S)   


# 4. Compute p-values (t-test)

pvalues <- apply(data, 1, function(x) {
  t.test(x[genotype == "E1-AF"], x[genotype == "E1-C"])$p.value
})


# 5. Adjust p-values (FDR)

padj <- p.adjust(pvalues, method = "fdr")


# 6. Create results_batch1_positive_susceptible_24hpi table

results_batch1_positive_susceptible_24hpi <- data.frame(
  metabolite = rownames(data),
  log2FC = log2FC,
  pvalue = pvalues,
  padj = padj
)

results_batch1_positive_susceptible_24hpi$negLog10P <- -log10(results_batch1_positive_susceptible_24hpi$padj)


# 7. Volcano plot

# Define significance
results_batch1_positive_susceptible_24hpi$significant <- "Not Significant"
results_batch1_positive_susceptible_24hpi$significant[results_batch1_positive_susceptible_24hpi$log2FC > 1 & results_batch1_positive_susceptible_24hpi$padj < 0.05] <- "Up"
results_batch1_positive_susceptible_24hpi$significant[results_batch1_positive_susceptible_24hpi$log2FC < -1 & results_batch1_positive_susceptible_24hpi$padj < 0.05] <- "Down"

Volcano_batch1_positive_susceptible_24hpi <- ggplot(results_batch1_positive_susceptible_24hpi, aes(x = log2FC, y = negLog10P, color = significant)) +
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
    data = subset(results_batch1_positive_susceptible_24hpi, significant != "Not Significant"),
    aes(label = metabolite)
  )+
  # Set fixed x and y axis limits
  scale_x_continuous(limits = c(-5, 5)) +  # Adjust min/max for log2FC
  scale_y_continuous(limits = c(0, 20))    # Adjust min/max for -log10(FDR)

results_batch1_positive_susceptible_24hpi$significant<-as.factor(results_batch1_positive_susceptible_24hpi$significant)
summary(results_batch1_positive_susceptible_24hpi)




#72hpi----
# Keep only rows where genotype is NOT "E1-R"
batch1_positive_susceptible <- batch1_positive[batch1_positive$genotype != "E1-R", ]
batch1_positive_susceptible_72hpi <- batch1_positive_susceptible[batch1_positive_susceptible$hpi == "72", ]
batch1_positive_susceptible_72hpi <- batch1_positive_susceptible_72hpi[, !(names(batch1_positive_susceptible_72hpi) %in% c("FEATURE_ID", "genotype", "hpi", "Sample_ID", "Injection_order", "Group", "QC"))]
batch1_positive_susceptible_72hpi <- t(batch1_positive_susceptible_72hpi)



# 2. Extract genotype and data

genotype <- as.character(unlist(batch1_positive_susceptible_72hpi[1, ])) # first row = genotype
data <- batch1_positive_susceptible_72hpi[-1, ]                          # remaining rows = metabolites
data <- apply(data, 2, as.numeric)            # convert to numeric
rownames(data) <- rownames(batch1_positive_susceptible_72hpi)[-1]

# 3. Compute log2 fold change

mean_R <- rowMeans(data[, genotype == "E1-AF"])
mean_S <- rowMeans(data[, genotype == "E1-C"])

log2FC <- log2(mean_R / mean_S)   


# 4. Compute p-values (t-test)

pvalues <- apply(data, 1, function(x) {
  t.test(x[genotype == "E1-AF"], x[genotype == "E1-C"])$p.value
})


# 5. Adjust p-values (FDR)

padj <- p.adjust(pvalues, method = "fdr")


# 6. Create results_batch1_positive_susceptible_72hpi table

results_batch1_positive_susceptible_72hpi <- data.frame(
  metabolite = rownames(data),
  log2FC = log2FC,
  pvalue = pvalues,
  padj = padj
)

results_batch1_positive_susceptible_72hpi$negLog10P <- -log10(results_batch1_positive_susceptible_72hpi$padj)


# 7. Volcano plot

# Define significance
results_batch1_positive_susceptible_72hpi$significant <- "Not Significant"
results_batch1_positive_susceptible_72hpi$significant[results_batch1_positive_susceptible_72hpi$log2FC > 1 & results_batch1_positive_susceptible_72hpi$padj < 0.05] <- "Up"
results_batch1_positive_susceptible_72hpi$significant[results_batch1_positive_susceptible_72hpi$log2FC < -1 & results_batch1_positive_susceptible_72hpi$padj < 0.05] <- "Down"

Volcano_batch1_positive_susceptible_72hpi <- ggplot(results_batch1_positive_susceptible_72hpi, aes(x = log2FC, y = negLog10P, color = significant)) +
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
    data = subset(results_batch1_positive_susceptible_72hpi, significant != "Not Significant"),
    aes(label = metabolite)
  )+
  # Set fixed x and y axis limits
  scale_x_continuous(limits = c(-5, 5)) +  # Adjust min/max for log2FC
  scale_y_continuous(limits = c(0, 20))    # Adjust min/max for -log10(FDR)

results_batch1_positive_susceptible_72hpi$significant<-as.factor(results_batch1_positive_susceptible_72hpi$significant)
summary(results_batch1_positive_susceptible_72hpi)

#120hpi----
# Keep only rows where genotype is NOT "E1-R"
batch1_positive_susceptible <- batch1_positive[batch1_positive$genotype != "E1-R", ]
batch1_positive_susceptible_120hpi <- batch1_positive_susceptible[batch1_positive_susceptible$hpi == "120", ]
batch1_positive_susceptible_120hpi <- batch1_positive_susceptible_120hpi[, !(names(batch1_positive_susceptible_120hpi) %in% c("FEATURE_ID", "genotype", "hpi", "Sample_ID", "Injection_order", "Group", "QC"))]
batch1_positive_susceptible_120hpi <- t(batch1_positive_susceptible_120hpi)



# 2. Extract genotype and data

genotype <- as.character(unlist(batch1_positive_susceptible_120hpi[1, ])) # first row = genotype
data <- batch1_positive_susceptible_120hpi[-1, ]                          # remaining rows = metabolites
data <- apply(data, 2, as.numeric)            # convert to numeric
rownames(data) <- rownames(batch1_positive_susceptible_120hpi)[-1]

# 3. Compute log2 fold change

mean_R <- rowMeans(data[, genotype == "E1-AF"])
mean_S <- rowMeans(data[, genotype == "E1-C"])

log2FC <- log2(mean_R / mean_S)   


# 4. Compute p-values (t-test)

pvalues <- apply(data, 1, function(x) {
  t.test(x[genotype == "E1-AF"], x[genotype == "E1-C"])$p.value
})


# 5. Adjust p-values (FDR)

padj <- p.adjust(pvalues, method = "fdr")


# 6. Create results_batch1_positive_susceptible_120hpi table

results_batch1_positive_susceptible_120hpi <- data.frame(
  metabolite = rownames(data),
  log2FC = log2FC,
  pvalue = pvalues,
  padj = padj
)

results_batch1_positive_susceptible_120hpi$negLog10P <- -log10(results_batch1_positive_susceptible_120hpi$padj)


# 7. Volcano plot

# Define significance
results_batch1_positive_susceptible_120hpi$significant <- "Not Significant"
results_batch1_positive_susceptible_120hpi$significant[results_batch1_positive_susceptible_120hpi$log2FC > 1 & results_batch1_positive_susceptible_120hpi$padj < 0.05] <- "Up"
results_batch1_positive_susceptible_120hpi$significant[results_batch1_positive_susceptible_120hpi$log2FC < -1 & results_batch1_positive_susceptible_120hpi$padj < 0.05] <- "Down"

Volcano_batch1_positive_susceptible_120hpi <- ggplot(results_batch1_positive_susceptible_120hpi, aes(x = log2FC, y = negLog10P, color = significant)) +
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
    data = subset(results_batch1_positive_susceptible_120hpi, significant != "Not Significant"),
    aes(label = metabolite)
  )+
  # Set fixed x and y axis limits
  scale_x_continuous(limits = c(-5, 5)) +  # Adjust min/max for log2FC
  scale_y_continuous(limits = c(0, 20))    # Adjust min/max for -log10(FDR)

results_batch1_positive_susceptible_120hpi$significant<-as.factor(results_batch1_positive_susceptible_120hpi$significant)
summary(results_batch1_positive_susceptible_120hpi)







##batch1_positive_resistant----
#batch1_positive_resistant----
# Keep only rows where genotype is NOT "E1-S"
batch1_positive_resistant <- batch1_positive[batch1_positive$genotype != "E1-S", ]
batch1_positive_resistant <- batch1_positive_resistant[, !(names(batch1_positive_resistant) %in% c("FEATURE_ID", "genotype", "hpi", "Sample_ID", "Injection_order", "Group", "QC"))]
batch1_positive_resistant <- t(batch1_positive_resistant)



# 2. Extract genotype and data

genotype <- as.character(unlist(batch1_positive_resistant[1, ])) # first row = genotype
data <- batch1_positive_resistant[-1, ]                          # remaining rows = metabolites
data <- apply(data, 2, as.numeric)            # convert to numeric
rownames(data) <- rownames(batch1_positive_resistant)[-1]

# 3. Compute log2 fold change

mean_R <- rowMeans(data[, genotype == "E1-AF"])
mean_S <- rowMeans(data[, genotype == "E1-C"])

log2FC <- log2(mean_R / mean_S)   


# 4. Compute p-values (t-test)

pvalues <- apply(data, 1, function(x) {
  t.test(x[genotype == "E1-AF"], x[genotype == "E1-C"])$p.value
})


# 5. Adjust p-values (FDR)

padj <- p.adjust(pvalues, method = "fdr")


# 6. Create results_batch1_positive_resistant table

results_batch1_positive_resistant <- data.frame(
  metabolite = rownames(data),
  log2FC = log2FC,
  pvalue = pvalues,
  padj = padj
)

results_batch1_positive_resistant$negLog10P <- -log10(results_batch1_positive_resistant$padj)


# 7. Volcano plot

# Define significance
results_batch1_positive_resistant$significant <- "Not Significant"
results_batch1_positive_resistant$significant[results_batch1_positive_resistant$log2FC > 1 & results_batch1_positive_resistant$padj < 0.05] <- "Up"
results_batch1_positive_resistant$significant[results_batch1_positive_resistant$log2FC < -1 & results_batch1_positive_resistant$padj < 0.05] <- "Down"

Volcano_batch1_positive_resistant <- ggplot(results_batch1_positive_resistant, aes(x = log2FC, y = negLog10P, color = significant)) +
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
    data = subset(results_batch1_positive_resistant, significant != "Not Significant"),
    aes(label = metabolite)
  )+
  # Set fixed x and y axis limits
  scale_x_continuous(limits = c(-5, 5)) +  # Adjust min/max for log2FC
  scale_y_continuous(limits = c(0, 20))    # Adjust min/max for -log10(FDR)

results_batch1_positive_resistant$significant<-as.factor(results_batch1_positive_resistant$significant)
summary(results_batch1_positive_resistant)




#early----
# Keep only rows where genotype is NOT "E1-S"
batch1_positive_resistant <- batch1_positive[batch1_positive$genotype != "E1-S", ]
batch1_positive_resistant_early <- batch1_positive_resistant[
  batch1_positive_resistant$hpi %in% c("3", "24"),]
batch1_positive_resistant_early <- batch1_positive_resistant_early[, !(names(batch1_positive_resistant_early) %in% c("FEATURE_ID", "genotype", "hpi", "Sample_ID", "Injection_order", "Group", "QC"))]
batch1_positive_resistant_early <- t(batch1_positive_resistant_early)



# 2. Extract genotype and data

genotype <- as.character(unlist(batch1_positive_resistant_early[1, ])) # first row = genotype
data <- batch1_positive_resistant_early[-1, ]                          # remaining rows = metabolites
data <- apply(data, 2, as.numeric)            # convert to numeric
rownames(data) <- rownames(batch1_positive_resistant_early)[-1]

# 3. Compute log2 fold change

mean_R <- rowMeans(data[, genotype == "E1-AF"])
mean_S <- rowMeans(data[, genotype == "E1-C"])

log2FC <- log2(mean_R / mean_S)   


# 4. Compute p-values (t-test)

pvalues <- apply(data, 1, function(x) {
  t.test(x[genotype == "E1-AF"], x[genotype == "E1-C"])$p.value
})


# 5. Adjust p-values (FDR)

padj <- p.adjust(pvalues, method = "fdr")


# 6. Create results_batch1_positive_resistant_early table

results_batch1_positive_resistant_early <- data.frame(
  metabolite = rownames(data),
  log2FC = log2FC,
  pvalue = pvalues,
  padj = padj
)

results_batch1_positive_resistant_early$negLog10P <- -log10(results_batch1_positive_resistant_early$padj)


# 7. Volcano plot

# Define significance
results_batch1_positive_resistant_early$significant <- "Not Significant"
results_batch1_positive_resistant_early$significant[results_batch1_positive_resistant_early$log2FC > 1 & results_batch1_positive_resistant_early$padj < 0.05] <- "Up"
results_batch1_positive_resistant_early$significant[results_batch1_positive_resistant_early$log2FC < -1 & results_batch1_positive_resistant_early$padj < 0.05] <- "Down"

Volcano_batch1_positive_resistant_early <- ggplot(results_batch1_positive_resistant_early, aes(x = log2FC, y = negLog10P, color = significant)) +
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
    data = subset(results_batch1_positive_resistant_early, significant != "Not Significant"),
    aes(label = metabolite)
  )+
  # Set fixed x and y axis limits
  scale_x_continuous(limits = c(-5, 5)) +  # Adjust min/max for log2FC
  scale_y_continuous(limits = c(0, 20))    # Adjust min/max for -log10(FDR)

results_batch1_positive_resistant_early$significant<-as.factor(results_batch1_positive_resistant_early$significant)
summary(results_batch1_positive_resistant_early)





#late----
# Keep only rows where genotype is NOT "E1-S"
batch1_positive_resistant <- batch1_positive[batch1_positive$genotype != "E1-S", ]
batch1_positive_resistant_late <- batch1_positive_resistant[
  batch1_positive_resistant$hpi %in% c("72", "120"),
]
batch1_positive_resistant_late <- batch1_positive_resistant_late[, !(names(batch1_positive_resistant_late) %in% c("FEATURE_ID", "genotype", "hpi", "Sample_ID", "Injection_order", "Group", "QC"))]
batch1_positive_resistant_late <- t(batch1_positive_resistant_late)



# 2. Extract genotype and data

genotype <- as.character(unlist(batch1_positive_resistant_late[1, ])) # first row = genotype
data <- batch1_positive_resistant_late[-1, ]                          # remaining rows = metabolites
data <- apply(data, 2, as.numeric)            # convert to numeric
rownames(data) <- rownames(batch1_positive_resistant_late)[-1]

# 3. Compute log2 fold change

mean_R <- rowMeans(data[, genotype == "E1-AF"])
mean_S <- rowMeans(data[, genotype == "E1-C"])

log2FC <- log2(mean_R / mean_S)   


# 4. Compute p-values (t-test)

pvalues <- apply(data, 1, function(x) {
  t.test(x[genotype == "E1-AF"], x[genotype == "E1-C"])$p.value
})


# 5. Adjust p-values (FDR)

padj <- p.adjust(pvalues, method = "fdr")


# 6. Create results_batch1_positive_resistant_late table

results_batch1_positive_resistant_late <- data.frame(
  metabolite = rownames(data),
  log2FC = log2FC,
  pvalue = pvalues,
  padj = padj
)

results_batch1_positive_resistant_late$negLog10P <- -log10(results_batch1_positive_resistant_late$padj)


# 7. Volcano plot

# Define significance
results_batch1_positive_resistant_late$significant <- "Not Significant"
results_batch1_positive_resistant_late$significant[results_batch1_positive_resistant_late$log2FC > 1 & results_batch1_positive_resistant_late$padj < 0.05] <- "Up"
results_batch1_positive_resistant_late$significant[results_batch1_positive_resistant_late$log2FC < -1 & results_batch1_positive_resistant_late$padj < 0.05] <- "Down"

Volcano_batch1_positive_resistant_late <- ggplot(results_batch1_positive_resistant_late, aes(x = log2FC, y = negLog10P, color = significant)) +
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
    data = subset(results_batch1_positive_resistant_late, significant != "Not Significant"),
    aes(label = metabolite)
  )+
  # Set fixed x and y axis limits
  scale_x_continuous(limits = c(-5, 5)) +  # Adjust min/max for log2FC
  scale_y_continuous(limits = c(0, 20))    # Adjust min/max for -log10(FDR)

results_batch1_positive_resistant_late$significant<-as.factor(results_batch1_positive_resistant_late$significant)
summary(results_batch1_positive_resistant_late)




#3hpi----
# Keep only rows where genotype is NOT "E1-S"
batch1_positive_resistant <- batch1_positive[batch1_positive$genotype != "E1-S", ]
batch1_positive_resistant_3hpi <- batch1_positive_resistant[batch1_positive_resistant$hpi == "3", ]
batch1_positive_resistant_3hpi <- batch1_positive_resistant_3hpi[, !(names(batch1_positive_resistant_3hpi) %in% c("FEATURE_ID", "genotype", "hpi", "Sample_ID", "Injection_order", "Group", "QC"))]
batch1_positive_resistant_3hpi <- t(batch1_positive_resistant_3hpi)



# 2. Extract genotype and data

genotype <- as.character(unlist(batch1_positive_resistant_3hpi[1, ])) # first row = genotype
data <- batch1_positive_resistant_3hpi[-1, ]                          # remaining rows = metabolites
data <- apply(data, 2, as.numeric)            # convert to numeric
rownames(data) <- rownames(batch1_positive_resistant_3hpi)[-1]

# 3. Compute log2 fold change

mean_R <- rowMeans(data[, genotype == "E1-AF"])
mean_S <- rowMeans(data[, genotype == "E1-C"])

log2FC <- log2(mean_R / mean_S)   


# 4. Compute p-values (t-test)

pvalues <- apply(data, 1, function(x) {
  t.test(x[genotype == "E1-AF"], x[genotype == "E1-C"])$p.value
})


# 5. Adjust p-values (FDR)

padj <- p.adjust(pvalues, method = "fdr")


# 6. Create results_batch1_positive_resistant_3hpi table

results_batch1_positive_resistant_3hpi <- data.frame(
  metabolite = rownames(data),
  log2FC = log2FC,
  pvalue = pvalues,
  padj = padj
)

results_batch1_positive_resistant_3hpi$negLog10P <- -log10(results_batch1_positive_resistant_3hpi$padj)


# 7. Volcano plot

# Define significance
results_batch1_positive_resistant_3hpi$significant <- "Not Significant"
results_batch1_positive_resistant_3hpi$significant[results_batch1_positive_resistant_3hpi$log2FC > 1 & results_batch1_positive_resistant_3hpi$padj < 0.05] <- "Up"
results_batch1_positive_resistant_3hpi$significant[results_batch1_positive_resistant_3hpi$log2FC < -1 & results_batch1_positive_resistant_3hpi$padj < 0.05] <- "Down"

Volcano_batch1_positive_resistant_3hpi <- ggplot(results_batch1_positive_resistant_3hpi, aes(x = log2FC, y = negLog10P, color = significant)) +
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
    data = subset(results_batch1_positive_resistant_3hpi, significant != "Not Significant"),
    aes(label = metabolite)
  )+
  # Set fixed x and y axis limits
  scale_x_continuous(limits = c(-5, 5)) +  # Adjust min/max for log2FC
  scale_y_continuous(limits = c(0, 20))    # Adjust min/max for -log10(FDR)

results_batch1_positive_resistant_3hpi$significant<-as.factor(results_batch1_positive_resistant_3hpi$significant)
summary(results_batch1_positive_resistant_3hpi)




#24hpi----
# Keep only rows where genotype is NOT "E1-S"
batch1_positive_resistant <- batch1_positive[batch1_positive$genotype != "E1-S", ]
batch1_positive_resistant_24hpi <- batch1_positive_resistant[batch1_positive_resistant$hpi == "24", ]
batch1_positive_resistant_24hpi <- batch1_positive_resistant_24hpi[, !(names(batch1_positive_resistant_24hpi) %in% c("FEATURE_ID", "genotype", "hpi", "Sample_ID", "Injection_order", "Group", "QC"))]
batch1_positive_resistant_24hpi <- t(batch1_positive_resistant_24hpi)



# 2. Extract genotype and data

genotype <- as.character(unlist(batch1_positive_resistant_24hpi[1, ])) # first row = genotype
data <- batch1_positive_resistant_24hpi[-1, ]                          # remaining rows = metabolites
data <- apply(data, 2, as.numeric)            # convert to numeric
rownames(data) <- rownames(batch1_positive_resistant_24hpi)[-1]

# 3. Compute log2 fold change

mean_R <- rowMeans(data[, genotype == "E1-AF"])
mean_S <- rowMeans(data[, genotype == "E1-C"])

log2FC <- log2(mean_R / mean_S)   


# 4. Compute p-values (t-test)

pvalues <- apply(data, 1, function(x) {
  t.test(x[genotype == "E1-AF"], x[genotype == "E1-C"])$p.value
})


# 5. Adjust p-values (FDR)

padj <- p.adjust(pvalues, method = "fdr")


# 6. Create results_batch1_positive_resistant_24hpi table

results_batch1_positive_resistant_24hpi <- data.frame(
  metabolite = rownames(data),
  log2FC = log2FC,
  pvalue = pvalues,
  padj = padj
)

results_batch1_positive_resistant_24hpi$negLog10P <- -log10(results_batch1_positive_resistant_24hpi$padj)


# 7. Volcano plot

# Define significance
results_batch1_positive_resistant_24hpi$significant <- "Not Significant"
results_batch1_positive_resistant_24hpi$significant[results_batch1_positive_resistant_24hpi$log2FC > 1 & results_batch1_positive_resistant_24hpi$padj < 0.05] <- "Up"
results_batch1_positive_resistant_24hpi$significant[results_batch1_positive_resistant_24hpi$log2FC < -1 & results_batch1_positive_resistant_24hpi$padj < 0.05] <- "Down"

Volcano_batch1_positive_resistant_24hpi <- ggplot(results_batch1_positive_resistant_24hpi, aes(x = log2FC, y = negLog10P, color = significant)) +
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
    data = subset(results_batch1_positive_resistant_24hpi, significant != "Not Significant"),
    aes(label = metabolite)
  )+
  # Set fixed x and y axis limits
  scale_x_continuous(limits = c(-5, 5)) +  # Adjust min/max for log2FC
  scale_y_continuous(limits = c(0, 20))    # Adjust min/max for -log10(FDR)

results_batch1_positive_resistant_24hpi$significant<-as.factor(results_batch1_positive_resistant_24hpi$significant)
summary(results_batch1_positive_resistant_24hpi)




#72hpi----
# Keep only rows where genotype is NOT "E1-S"
batch1_positive_resistant <- batch1_positive[batch1_positive$genotype != "E1-S", ]
batch1_positive_resistant_72hpi <- batch1_positive_resistant[batch1_positive_resistant$hpi == "72", ]
batch1_positive_resistant_72hpi <- batch1_positive_resistant_72hpi[, !(names(batch1_positive_resistant_72hpi) %in% c("FEATURE_ID", "genotype", "hpi", "Sample_ID", "Injection_order", "Group", "QC"))]
batch1_positive_resistant_72hpi <- t(batch1_positive_resistant_72hpi)



# 2. Extract genotype and data

genotype <- as.character(unlist(batch1_positive_resistant_72hpi[1, ])) # first row = genotype
data <- batch1_positive_resistant_72hpi[-1, ]                          # remaining rows = metabolites
data <- apply(data, 2, as.numeric)            # convert to numeric
rownames(data) <- rownames(batch1_positive_resistant_72hpi)[-1]

# 3. Compute log2 fold change

mean_R <- rowMeans(data[, genotype == "E1-AF"])
mean_S <- rowMeans(data[, genotype == "E1-C"])

log2FC <- log2(mean_R / mean_S)   


# 4. Compute p-values (t-test)

pvalues <- apply(data, 1, function(x) {
  t.test(x[genotype == "E1-AF"], x[genotype == "E1-C"])$p.value
})


# 5. Adjust p-values (FDR)

padj <- p.adjust(pvalues, method = "fdr")


# 6. Create results_batch1_positive_resistant_72hpi table

results_batch1_positive_resistant_72hpi <- data.frame(
  metabolite = rownames(data),
  log2FC = log2FC,
  pvalue = pvalues,
  padj = padj
)

results_batch1_positive_resistant_72hpi$negLog10P <- -log10(results_batch1_positive_resistant_72hpi$padj)


# 7. Volcano plot

# Define significance
results_batch1_positive_resistant_72hpi$significant <- "Not Significant"
results_batch1_positive_resistant_72hpi$significant[results_batch1_positive_resistant_72hpi$log2FC > 1 & results_batch1_positive_resistant_72hpi$padj < 0.05] <- "Up"
results_batch1_positive_resistant_72hpi$significant[results_batch1_positive_resistant_72hpi$log2FC < -1 & results_batch1_positive_resistant_72hpi$padj < 0.05] <- "Down"

Volcano_batch1_positive_resistant_72hpi <- ggplot(results_batch1_positive_resistant_72hpi, aes(x = log2FC, y = negLog10P, color = significant)) +
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
    data = subset(results_batch1_positive_resistant_72hpi, significant != "Not Significant"),
    aes(label = metabolite)
  )+
  # Set fixed x and y axis limits
  scale_x_continuous(limits = c(-5, 5)) +  # Adjust min/max for log2FC
  scale_y_continuous(limits = c(0, 20))    # Adjust min/max for -log10(FDR)

results_batch1_positive_resistant_72hpi$significant<-as.factor(results_batch1_positive_resistant_72hpi$significant)
summary(results_batch1_positive_resistant_72hpi)

#120hpi----
# Keep only rows where genotype is NOT "E1-S"
batch1_positive_resistant <- batch1_positive[batch1_positive$genotype != "E1-S", ]
batch1_positive_resistant_120hpi <- batch1_positive_resistant[batch1_positive_resistant$hpi == "120", ]
batch1_positive_resistant_120hpi <- batch1_positive_resistant_120hpi[, !(names(batch1_positive_resistant_120hpi) %in% c("FEATURE_ID", "genotype", "hpi", "Sample_ID", "Injection_order", "Group", "QC"))]
batch1_positive_resistant_120hpi <- t(batch1_positive_resistant_120hpi)



# 2. Extract genotype and data

genotype <- as.character(unlist(batch1_positive_resistant_120hpi[1, ])) # first row = genotype
data <- batch1_positive_resistant_120hpi[-1, ]                          # remaining rows = metabolites
data <- apply(data, 2, as.numeric)            # convert to numeric
rownames(data) <- rownames(batch1_positive_resistant_120hpi)[-1]

# 3. Compute log2 fold change

mean_R <- rowMeans(data[, genotype == "E1-AF"])
mean_S <- rowMeans(data[, genotype == "E1-C"])

log2FC <- log2(mean_R / mean_S)   


# 4. Compute p-values (t-test)

pvalues <- apply(data, 1, function(x) {
  t.test(x[genotype == "E1-AF"], x[genotype == "E1-C"])$p.value
})


# 5. Adjust p-values (FDR)

padj <- p.adjust(pvalues, method = "fdr")


# 6. Create results_batch1_positive_resistant_120hpi table

results_batch1_positive_resistant_120hpi <- data.frame(
  metabolite = rownames(data),
  log2FC = log2FC,
  pvalue = pvalues,
  padj = padj
)

results_batch1_positive_resistant_120hpi$negLog10P <- -log10(results_batch1_positive_resistant_120hpi$padj)


# 7. Volcano plot

# Define significance
results_batch1_positive_resistant_120hpi$significant <- "Not Significant"
results_batch1_positive_resistant_120hpi$significant[results_batch1_positive_resistant_120hpi$log2FC > 1 & results_batch1_positive_resistant_120hpi$padj < 0.05] <- "Up"
results_batch1_positive_resistant_120hpi$significant[results_batch1_positive_resistant_120hpi$log2FC < -1 & results_batch1_positive_resistant_120hpi$padj < 0.05] <- "Down"

Volcano_batch1_positive_resistant_120hpi <- ggplot(results_batch1_positive_resistant_120hpi, aes(x = log2FC, y = negLog10P, color = significant)) +
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
    data = subset(results_batch1_positive_resistant_120hpi, significant != "Not Significant"),
    aes(label = metabolite)
  )+
  # Set fixed x and y axis limits
  scale_x_continuous(limits = c(-5, 5)) +  # Adjust min/max for log2FC
  scale_y_continuous(limits = c(0, 20))    # Adjust min/max for -log10(FDR)

results_batch1_positive_resistant_120hpi$significant<-as.factor(results_batch1_positive_resistant_120hpi$significant)
summary(results_batch1_positive_resistant_120hpi)





#####

results_batch1_positive_susceptible<-results_batch1_positive_susceptible[,c(1,6)]
colnames(results_batch1_positive_susceptible)[colnames(results_batch1_positive_susceptible) == "significant"] <- "batch1_positive_susceptible"
results_batch1_positive_resistant<-results_batch1_positive_resistant[,c(1,6)]
colnames(results_batch1_positive_resistant)[colnames(results_batch1_positive_resistant) == "significant"] <- "batch1_positive_resistant"

results_batch1_positive_susceptible_early<-results_batch1_positive_susceptible_early[,c(1,6)]
colnames(results_batch1_positive_susceptible_early)[colnames(results_batch1_positive_susceptible_early) == "significant"] <- "batch1_positive_susceptible_early"

results_batch1_positive_resistant_early<-results_batch1_positive_resistant_early[,c(1,6)]
colnames(results_batch1_positive_resistant_early)[colnames(results_batch1_positive_resistant_early) == "significant"] <- "batch1_positive_resistant_early"


results_batch1_positive_susceptible_late<-results_batch1_positive_susceptible_late[,c(1,6)]
colnames(results_batch1_positive_susceptible_late)[colnames(results_batch1_positive_susceptible_late) == "significant"] <- "batch1_positive_susceptible_late"
results_batch1_positive_resistant_late<-results_batch1_positive_resistant_late[,c(1,6)]
colnames(results_batch1_positive_resistant_late)[colnames(results_batch1_positive_resistant_late) == "significant"] <- "batch1_positive_resistant_late"
results_batch1_positive_susceptible_3hpi<-results_batch1_positive_susceptible_3hpi[,c(1,6)]
colnames(results_batch1_positive_susceptible_3hpi)[colnames(results_batch1_positive_susceptible_3hpi) == "significant"] <- "batch1_positive_susceptible_3hpi"
results_batch1_positive_resistant_3hpi<-results_batch1_positive_resistant_3hpi[,c(1,6)]
colnames(results_batch1_positive_resistant_3hpi)[colnames(results_batch1_positive_resistant_3hpi) == "significant"] <- "batch1_positive_resistant_3hpi"
results_batch1_positive_susceptible_24hpi<-results_batch1_positive_susceptible_24hpi[,c(1,6)]
colnames(results_batch1_positive_susceptible_24hpi)[colnames(results_batch1_positive_susceptible_24hpi) == "significant"] <- "batch1_positive_susceptible_24hpi"
results_batch1_positive_resistant_24hpi<-results_batch1_positive_resistant_24hpi[,c(1,6)]
colnames(results_batch1_positive_resistant_24hpi)[colnames(results_batch1_positive_resistant_24hpi) == "significant"] <- "batch1_positive_resistant_24hpi"
results_batch1_positive_susceptible_72hpi<-results_batch1_positive_susceptible_72hpi[,c(1,6)]
colnames(results_batch1_positive_susceptible_72hpi)[colnames(results_batch1_positive_susceptible_72hpi) == "significant"] <- "batch1_positive_susceptible_72hpi"
results_batch1_positive_resistant_72hpi<-results_batch1_positive_resistant_72hpi[,c(1,6)]
colnames(results_batch1_positive_resistant_72hpi)[colnames(results_batch1_positive_resistant_72hpi) == "significant"] <- "batch1_positive_resistant_72hpi"
results_batch1_positive_susceptible_120hpi<-results_batch1_positive_susceptible_120hpi[,c(1,6)]
colnames(results_batch1_positive_susceptible_120hpi)[colnames(results_batch1_positive_susceptible_120hpi) == "significant"] <- "batch1_positive_susceptible_120hpi"
results_batch1_positive_resistant_120hpi<-results_batch1_positive_resistant_120hpi[,c(1,6)]
colnames(results_batch1_positive_resistant_120hpi)[colnames(results_batch1_positive_resistant_120hpi) == "significant"] <- "batch1_positive_resistant_120hpi"




all_results <- list(
  results_batch1_positive_susceptible,
  results_batch1_positive_resistant,
  results_batch1_positive_susceptible_early,
  results_batch1_positive_resistant_early,
  results_batch1_positive_susceptible_late,
  results_batch1_positive_resistant_late,
  results_batch1_positive_susceptible_3hpi,
  results_batch1_positive_resistant_3hpi,
  results_batch1_positive_susceptible_24hpi,
  results_batch1_positive_resistant_24hpi,
  results_batch1_positive_susceptible_72hpi,
  results_batch1_positive_resistant_72hpi,
  results_batch1_positive_susceptible_120hpi,
  results_batch1_positive_resistant_120hpi
)

merged_results_results_batch1_positive <- reduce(all_results, full_join, by = "metabolite")
merged_results_results_batch1_positive <- merged_results_results_batch1_positive[
  rowSums(merged_results_results_batch1_positive[, -1] != "Not Significant", na.rm = TRUE) > 0,
]

merged_results_results_batch1_positive$metabolite <- paste0(
  merged_results_results_batch1_positive$metabolite,
  "pos")

heatmap_batch1_positive <- merged_results_results_batch1_positive



heatmap_batch1_positive <- as.data.frame(lapply(heatmap_batch1_positive, as.character))
rownames(heatmap_batch1_positive) <- heatmap_batch1_positive$metabolite
heatmap_batch1_positive <- heatmap_batch1_positive[, -1]

heatmap_batch1_positive[] <- lapply(heatmap_batch1_positive, function(x) {
  x[x == "Up"] <- 1
  x[x == "Down"] <- -1
  x[x == "Not Significant"] <- 0
  return(as.numeric(x))
})

# Convert to matrix
heatmap_matrix <- as.matrix(heatmap_batch1_positive)

# Convert to numeric
mode(heatmap_matrix) <- "numeric"

# compute row sums (total signal per metabolite)
row_score <- rowSums(heatmap_matrix, na.rm = TRUE)

# order from highest to lowest
heatmap_matrix_batch1_positive <- heatmap_matrix[order(row_score, decreasing = TRUE), ]

# Plot heatmap
pheatmap(
  heatmap_matrix_batch1_positive,
  color = c("blue", "white", "red"),
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  border_color = "grey80"
)

###batch_2_positive----

batch2_positive <- read_excel("final_filtered_files/batch_2_positive_Intensity_filtered_notame.xlsx", col_names = FALSE)

batch2_positive <- batch2_positive %>%
  mutate(across(35, ~ ifelse(!is.na(batch2_positive[[3]]) & batch2_positive[[3]] != "", batch2_positive[[3]], .)))
batch2_positive <- batch2_positive[,-c(1:34)]
batch2_positive <- t(batch2_positive)
batch2_positive <- as.data.frame(batch2_positive)
colnames(batch2_positive) <- batch2_positive[1, ]
batch2_positive <- batch2_positive[-1, ]

batch2_positive_row_name <- batch2_positive %>%
  filter(Group != "Blank") %>%
  separate(Group, into = c("Batch","pheno","treatmentA","hpi"), sep = "-") %>%
  mutate(genotype = paste(Batch, pheno, sep = "-")) %>%
  mutate(treatment = paste(Batch, treatmentA, sep = "-")) %>%
  select(FEATURE_ID,genotype,treatment, hpi)
batch2_positive <- merge(batch2_positive_row_name,batch2_positive, by ="FEATURE_ID")
batch2_positive[is.na(batch2_positive)] <- 0 
batch2_positive <- batch2_positive %>%
  filter(!is.na(hpi) & !is.na(hpi))

##batch2_positive_susceptible----
#batch2_positive_susceptible----
# Keep only rows where genotype is NOT "E2-R"
batch2_positive_susceptible <- batch2_positive[batch2_positive$genotype != "E2-R", ]
batch2_positive_susceptible <- batch2_positive_susceptible[, !(names(batch2_positive_susceptible) %in% c("FEATURE_ID", "genotype", "hpi", "Sample_ID", "Injection_order", "Group", "QC"))]
batch2_positive_susceptible <- t(batch2_positive_susceptible)



# 2. Extract genotype and data

genotype <- as.character(unlist(batch2_positive_susceptible[1, ])) # first row = genotype
data <- batch2_positive_susceptible[-1, ]                          # remaining rows = metabolites
data <- apply(data, 2, as.numeric)            # convert to numeric
rownames(data) <- rownames(batch2_positive_susceptible)[-1]

# 3. Compute log2 fold change

mean_R <- rowMeans(data[, genotype == "E2-AP"])
mean_S <- rowMeans(data[, genotype == "E2-C"])

log2FC <- log2(mean_R / mean_S)   


# 4. Compute p-values (t-test)

pvalues <- apply(data, 1, function(x) {
  t.test(x[genotype == "E2-AP"], x[genotype == "E2-C"])$p.value
})


# 5. Adjust p-values (FDR)

padj <- p.adjust(pvalues, method = "fdr")


# 6. Create results_batch2_positive_susceptible table

results_batch2_positive_susceptible <- data.frame(
  metabolite = rownames(data),
  log2FC = log2FC,
  pvalue = pvalues,
  padj = padj
)

results_batch2_positive_susceptible$negLog10P <- -log10(results_batch2_positive_susceptible$padj)


# 7. Volcano plot

# Define significance
results_batch2_positive_susceptible$significant <- "Not Significant"
results_batch2_positive_susceptible$significant[results_batch2_positive_susceptible$log2FC > 1 & results_batch2_positive_susceptible$padj < 0.05] <- "Up"
results_batch2_positive_susceptible$significant[results_batch2_positive_susceptible$log2FC < -1 & results_batch2_positive_susceptible$padj < 0.05] <- "Down"

Volcano_batch2_positive_susceptible <- ggplot(results_batch2_positive_susceptible, aes(x = log2FC, y = negLog10P, color = significant)) +
  geom_point(size = 3, alpha = 0.8) +
  scale_color_manual(values = c("blue", "grey", "red")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  theme_minimal() +
  labs(
    x = "Log2 Fold Change (E2-AP / E2-C)",
    y = "-Log10(FDR)",
    title = "Volcano Plot"
  ) +
  geom_text_repel(
    data = subset(results_batch2_positive_susceptible, significant != "Not Significant"),
    aes(label = metabolite)
  )+
  # Set fixed x and y axis limits
  scale_x_continuous(limits = c(-5, 5)) +  # Adjust min/max for log2FC
  scale_y_continuous(limits = c(0, 10))    # Adjust min/max for -log10(FDR)

results_batch2_positive_susceptible$significant<-as.factor(results_batch2_positive_susceptible$significant)
summary(results_batch2_positive_susceptible)




#early----
# Keep only rows where genotype is NOT "E2-R"
batch2_positive_susceptible <- batch2_positive[batch2_positive$genotype != "E2-R", ]
batch2_positive_susceptible_early <- batch2_positive_susceptible[
  batch2_positive_susceptible$hpi %in% c("3", "24"),]
batch2_positive_susceptible_early <- batch2_positive_susceptible_early[, !(names(batch2_positive_susceptible_early) %in% c("FEATURE_ID", "genotype", "hpi", "Sample_ID", "Injection_order", "Group", "QC"))]
batch2_positive_susceptible_early <- t(batch2_positive_susceptible_early)


# 2. Extract genotype and data

genotype <- as.character(unlist(batch2_positive_susceptible_early[1, ])) # first row = genotype
data <- batch2_positive_susceptible_early[-1, ]                          # remaining rows = metabolites
data <- apply(data, 2, as.numeric)            # convert to numeric
rownames(data) <- rownames(batch2_positive_susceptible_early)[-1]

# 3. Compute log2 fold change

mean_R <- rowMeans(data[, genotype == "E2-AP"])
mean_S <- rowMeans(data[, genotype == "E2-C"])

log2FC <- log2(mean_R / mean_S)   


# 4. Compute p-values (t-test)

pvalues <- apply(data, 1, function(x) {
  t.test(x[genotype == "E2-AP"], x[genotype == "E2-C"])$p.value
})


# 5. Adjust p-values (FDR)

padj <- p.adjust(pvalues, method = "fdr")


# 6. Create results_batch2_positive_susceptible_early table

results_batch2_positive_susceptible_early <- data.frame(
  metabolite = rownames(data),
  log2FC = log2FC,
  pvalue = pvalues,
  padj = padj
)

results_batch2_positive_susceptible_early$negLog10P <- -log10(results_batch2_positive_susceptible_early$padj)


# 7. Volcano plot

# Define significance
results_batch2_positive_susceptible_early$significant <- "Not Significant"
results_batch2_positive_susceptible_early$significant[results_batch2_positive_susceptible_early$log2FC > 1 & results_batch2_positive_susceptible_early$padj < 0.05] <- "Up"
results_batch2_positive_susceptible_early$significant[results_batch2_positive_susceptible_early$log2FC < -1 & results_batch2_positive_susceptible_early$padj < 0.05] <- "Down"

Volcano_batch2_positive_susceptible_early <- ggplot(results_batch2_positive_susceptible_early, aes(x = log2FC, y = negLog10P, color = significant)) +
  geom_point(size = 3, alpha = 0.8) +
  scale_color_manual(values = c("blue", "grey", "red")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  theme_minimal() +
  labs(
    x = "Log2 Fold Change (E2-AP / E2-C)",
    y = "-Log10(FDR)",
    title = "Volcano Plot"
  ) +
  geom_text_repel(
    data = subset(results_batch2_positive_susceptible_early, significant != "Not Significant"),
    aes(label = metabolite)
  )+
  # Set fixed x and y axis limits
  scale_x_continuous(limits = c(-5, 5)) +  # Adjust min/max for log2FC
  scale_y_continuous(limits = c(0, 10))    # Adjust min/max for -log10(FDR)

results_batch2_positive_susceptible_early$significant<-as.factor(results_batch2_positive_susceptible_early$significant)
summary(results_batch2_positive_susceptible_early)





#late----
# Keep only rows where genotype is NOT "E2-R"
batch2_positive_susceptible <- batch2_positive[batch2_positive$genotype != "E2-R", ]
batch2_positive_susceptible_late <- batch2_positive_susceptible[
  batch2_positive_susceptible$hpi %in% c("72", "120"),
]
batch2_positive_susceptible_late <- batch2_positive_susceptible_late[, !(names(batch2_positive_susceptible_late) %in% c("FEATURE_ID", "genotype", "hpi", "Sample_ID", "Injection_order", "Group", "QC"))]
batch2_positive_susceptible_late <- t(batch2_positive_susceptible_late)



# 2. Extract genotype and data

genotype <- as.character(unlist(batch2_positive_susceptible_late[1, ])) # first row = genotype
data <- batch2_positive_susceptible_late[-1, ]                          # remaining rows = metabolites
data <- apply(data, 2, as.numeric)            # convert to numeric
rownames(data) <- rownames(batch2_positive_susceptible_late)[-1]

# 3. Compute log2 fold change

mean_R <- rowMeans(data[, genotype == "E2-AP"])
mean_S <- rowMeans(data[, genotype == "E2-C"])

log2FC <- log2(mean_R / mean_S)   


# 4. Compute p-values (t-test)

pvalues <- apply(data, 1, function(x) {
  t.test(x[genotype == "E2-AP"], x[genotype == "E2-C"])$p.value
})


# 5. Adjust p-values (FDR)

padj <- p.adjust(pvalues, method = "fdr")


# 6. Create results_batch2_positive_susceptible_late table

results_batch2_positive_susceptible_late <- data.frame(
  metabolite = rownames(data),
  log2FC = log2FC,
  pvalue = pvalues,
  padj = padj
)

results_batch2_positive_susceptible_late$negLog10P <- -log10(results_batch2_positive_susceptible_late$padj)


# 7. Volcano plot

# Define significance
results_batch2_positive_susceptible_late$significant <- "Not Significant"
results_batch2_positive_susceptible_late$significant[results_batch2_positive_susceptible_late$log2FC > 1 & results_batch2_positive_susceptible_late$padj < 0.05] <- "Up"
results_batch2_positive_susceptible_late$significant[results_batch2_positive_susceptible_late$log2FC < -1 & results_batch2_positive_susceptible_late$padj < 0.05] <- "Down"

Volcano_batch2_positive_susceptible_late <- ggplot(results_batch2_positive_susceptible_late, aes(x = log2FC, y = negLog10P, color = significant)) +
  geom_point(size = 3, alpha = 0.8) +
  scale_color_manual(values = c("blue", "grey", "red")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  theme_minimal() +
  labs(
    x = "Log2 Fold Change (E2-AP / E2-C)",
    y = "-Log10(FDR)",
    title = "Volcano Plot"
  ) +
  geom_text_repel(
    data = subset(results_batch2_positive_susceptible_late, significant != "Not Significant"),
    aes(label = metabolite)
  )+
  # Set fixed x and y axis limits
  scale_x_continuous(limits = c(-5, 5)) +  # Adjust min/max for log2FC
  scale_y_continuous(limits = c(0, 10))    # Adjust min/max for -log10(FDR)

results_batch2_positive_susceptible_late$significant<-as.factor(results_batch2_positive_susceptible_late$significant)
summary(results_batch2_positive_susceptible_late)




#3hpi----
# Keep only rows where genotype is NOT "E2-R"
batch2_positive_susceptible <- batch2_positive[batch2_positive$genotype != "E2-R", ]
batch2_positive_susceptible_3hpi <- batch2_positive_susceptible[batch2_positive_susceptible$hpi == "3", ]
batch2_positive_susceptible_3hpi <- batch2_positive_susceptible_3hpi[, !(names(batch2_positive_susceptible_3hpi) %in% c("FEATURE_ID", "genotype", "hpi", "Sample_ID", "Injection_order", "Group", "QC"))]
batch2_positive_susceptible_3hpi <- t(batch2_positive_susceptible_3hpi)



# 2. Extract genotype and data

genotype <- as.character(unlist(batch2_positive_susceptible_3hpi[1, ])) # first row = genotype
data <- batch2_positive_susceptible_3hpi[-1, ]                          # remaining rows = metabolites
data <- apply(data, 2, as.numeric)            # convert to numeric
rownames(data) <- rownames(batch2_positive_susceptible_3hpi)[-1]

# 3. Compute log2 fold change

mean_R <- rowMeans(data[, genotype == "E2-AP"])
mean_S <- rowMeans(data[, genotype == "E2-C"])

log2FC <- log2(mean_R / mean_S)   


# 4. Compute p-values (t-test)

pvalues <- apply(data, 1, function(x) {
  t.test(x[genotype == "E2-AP"], x[genotype == "E2-C"])$p.value
})


# 5. Adjust p-values (FDR)

padj <- p.adjust(pvalues, method = "fdr")


# 6. Create results_batch2_positive_susceptible_3hpi table

results_batch2_positive_susceptible_3hpi <- data.frame(
  metabolite = rownames(data),
  log2FC = log2FC,
  pvalue = pvalues,
  padj = padj
)

results_batch2_positive_susceptible_3hpi$negLog10P <- -log10(results_batch2_positive_susceptible_3hpi$padj)


# 7. Volcano plot

# Define significance
results_batch2_positive_susceptible_3hpi$significant <- "Not Significant"
results_batch2_positive_susceptible_3hpi$significant[results_batch2_positive_susceptible_3hpi$log2FC > 1 & results_batch2_positive_susceptible_3hpi$padj < 0.05] <- "Up"
results_batch2_positive_susceptible_3hpi$significant[results_batch2_positive_susceptible_3hpi$log2FC < -1 & results_batch2_positive_susceptible_3hpi$padj < 0.05] <- "Down"

Volcano_batch2_positive_susceptible_3hpi <- ggplot(results_batch2_positive_susceptible_3hpi, aes(x = log2FC, y = negLog10P, color = significant)) +
  geom_point(size = 3, alpha = 0.8) +
  scale_color_manual(values = c("blue", "grey", "red")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  theme_minimal() +
  labs(
    x = "Log2 Fold Change (E2-AP / E2-C)",
    y = "-Log10(FDR)",
    title = "Volcano Plot"
  ) +
  geom_text_repel(
    data = subset(results_batch2_positive_susceptible_3hpi, significant != "Not Significant"),
    aes(label = metabolite)
  )+
  # Set fixed x and y axis limits
  scale_x_continuous(limits = c(-5, 5)) +  # Adjust min/max for log2FC
  scale_y_continuous(limits = c(0, 10))    # Adjust min/max for -log10(FDR)

results_batch2_positive_susceptible_3hpi$significant<-as.factor(results_batch2_positive_susceptible_3hpi$significant)
summary(results_batch2_positive_susceptible_3hpi)




#24hpi----
# Keep only rows where genotype is NOT "E2-R"
batch2_positive_susceptible <- batch2_positive[batch2_positive$genotype != "E2-R", ]
batch2_positive_susceptible_24hpi <- batch2_positive_susceptible[batch2_positive_susceptible$hpi == "24", ]
batch2_positive_susceptible_24hpi <- batch2_positive_susceptible_24hpi[, !(names(batch2_positive_susceptible_24hpi) %in% c("FEATURE_ID", "genotype", "hpi", "Sample_ID", "Injection_order", "Group", "QC"))]
batch2_positive_susceptible_24hpi <- t(batch2_positive_susceptible_24hpi)



# 2. Extract genotype and data

genotype <- as.character(unlist(batch2_positive_susceptible_24hpi[1, ])) # first row = genotype
data <- batch2_positive_susceptible_24hpi[-1, ]                          # remaining rows = metabolites
data <- apply(data, 2, as.numeric)            # convert to numeric
rownames(data) <- rownames(batch2_positive_susceptible_24hpi)[-1]

# 3. Compute log2 fold change

mean_R <- rowMeans(data[, genotype == "E2-AP"])
mean_S <- rowMeans(data[, genotype == "E2-C"])

log2FC <- log2(mean_R / mean_S)   


# 4. Compute p-values (t-test)

pvalues <- apply(data, 1, function(x) {
  t.test(x[genotype == "E2-AP"], x[genotype == "E2-C"])$p.value
})


# 5. Adjust p-values (FDR)

padj <- p.adjust(pvalues, method = "fdr")


# 6. Create results_batch2_positive_susceptible_24hpi table

results_batch2_positive_susceptible_24hpi <- data.frame(
  metabolite = rownames(data),
  log2FC = log2FC,
  pvalue = pvalues,
  padj = padj
)

results_batch2_positive_susceptible_24hpi$negLog10P <- -log10(results_batch2_positive_susceptible_24hpi$padj)


# 7. Volcano plot

# Define significance
results_batch2_positive_susceptible_24hpi$significant <- "Not Significant"
results_batch2_positive_susceptible_24hpi$significant[results_batch2_positive_susceptible_24hpi$log2FC > 1 & results_batch2_positive_susceptible_24hpi$padj < 0.05] <- "Up"
results_batch2_positive_susceptible_24hpi$significant[results_batch2_positive_susceptible_24hpi$log2FC < -1 & results_batch2_positive_susceptible_24hpi$padj < 0.05] <- "Down"

Volcano_batch2_positive_susceptible_24hpi <- ggplot(results_batch2_positive_susceptible_24hpi, aes(x = log2FC, y = negLog10P, color = significant)) +
  geom_point(size = 3, alpha = 0.8) +
  scale_color_manual(values = c("blue", "grey", "red")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  theme_minimal() +
  labs(
    x = "Log2 Fold Change (E2-AP / E2-C)",
    y = "-Log10(FDR)",
    title = "Volcano Plot"
  ) +
  geom_text_repel(
    data = subset(results_batch2_positive_susceptible_24hpi, significant != "Not Significant"),
    aes(label = metabolite)
  )+
  # Set fixed x and y axis limits
  scale_x_continuous(limits = c(-5, 5)) +  # Adjust min/max for log2FC
  scale_y_continuous(limits = c(0, 10))    # Adjust min/max for -log10(FDR)

results_batch2_positive_susceptible_24hpi$significant<-as.factor(results_batch2_positive_susceptible_24hpi$significant)
summary(results_batch2_positive_susceptible_24hpi)




#72hpi----
# Keep only rows where genotype is NOT "E2-R"
batch2_positive_susceptible <- batch2_positive[batch2_positive$genotype != "E2-R", ]
batch2_positive_susceptible_72hpi <- batch2_positive_susceptible[batch2_positive_susceptible$hpi == "72", ]
batch2_positive_susceptible_72hpi <- batch2_positive_susceptible_72hpi[, !(names(batch2_positive_susceptible_72hpi) %in% c("FEATURE_ID", "genotype", "hpi", "Sample_ID", "Injection_order", "Group", "QC"))]
batch2_positive_susceptible_72hpi <- t(batch2_positive_susceptible_72hpi)



# 2. Extract genotype and data

genotype <- as.character(unlist(batch2_positive_susceptible_72hpi[1, ])) # first row = genotype
data <- batch2_positive_susceptible_72hpi[-1, ]                          # remaining rows = metabolites
data <- apply(data, 2, as.numeric)            # convert to numeric
rownames(data) <- rownames(batch2_positive_susceptible_72hpi)[-1]

# 3. Compute log2 fold change

mean_R <- rowMeans(data[, genotype == "E2-AP"])
mean_S <- rowMeans(data[, genotype == "E2-C"])

log2FC <- log2(mean_R / mean_S)   


# 4. Compute p-values (t-test)

pvalues <- apply(data, 1, function(x) {
  t.test(x[genotype == "E2-AP"], x[genotype == "E2-C"])$p.value
})


# 5. Adjust p-values (FDR)

padj <- p.adjust(pvalues, method = "fdr")


# 6. Create results_batch2_positive_susceptible_72hpi table

results_batch2_positive_susceptible_72hpi <- data.frame(
  metabolite = rownames(data),
  log2FC = log2FC,
  pvalue = pvalues,
  padj = padj
)

results_batch2_positive_susceptible_72hpi$negLog10P <- -log10(results_batch2_positive_susceptible_72hpi$padj)


# 7. Volcano plot

# Define significance
results_batch2_positive_susceptible_72hpi$significant <- "Not Significant"
results_batch2_positive_susceptible_72hpi$significant[results_batch2_positive_susceptible_72hpi$log2FC > 1 & results_batch2_positive_susceptible_72hpi$padj < 0.05] <- "Up"
results_batch2_positive_susceptible_72hpi$significant[results_batch2_positive_susceptible_72hpi$log2FC < -1 & results_batch2_positive_susceptible_72hpi$padj < 0.05] <- "Down"

Volcano_batch2_positive_susceptible_72hpi <- ggplot(results_batch2_positive_susceptible_72hpi, aes(x = log2FC, y = negLog10P, color = significant)) +
  geom_point(size = 3, alpha = 0.8) +
  scale_color_manual(values = c("blue", "grey", "red")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  theme_minimal() +
  labs(
    x = "Log2 Fold Change (E2-AP / E2-C)",
    y = "-Log10(FDR)",
    title = "Volcano Plot"
  ) +
  geom_text_repel(
    data = subset(results_batch2_positive_susceptible_72hpi, significant != "Not Significant"),
    aes(label = metabolite)
  )+
  # Set fixed x and y axis limits
  scale_x_continuous(limits = c(-5, 5)) +  # Adjust min/max for log2FC
  scale_y_continuous(limits = c(0, 10))    # Adjust min/max for -log10(FDR)

results_batch2_positive_susceptible_72hpi$significant<-as.factor(results_batch2_positive_susceptible_72hpi$significant)
summary(results_batch2_positive_susceptible_72hpi)

#120hpi----
# Keep only rows where genotype is NOT "E2-R"
batch2_positive_susceptible <- batch2_positive[batch2_positive$genotype != "E2-R", ]
batch2_positive_susceptible_120hpi <- batch2_positive_susceptible[batch2_positive_susceptible$hpi == "120", ]
batch2_positive_susceptible_120hpi <- batch2_positive_susceptible_120hpi[, !(names(batch2_positive_susceptible_120hpi) %in% c("FEATURE_ID", "genotype", "hpi", "Sample_ID", "Injection_order", "Group", "QC"))]
batch2_positive_susceptible_120hpi <- t(batch2_positive_susceptible_120hpi)



# 2. Extract genotype and data

genotype <- as.character(unlist(batch2_positive_susceptible_120hpi[1, ])) # first row = genotype
data <- batch2_positive_susceptible_120hpi[-1, ]                          # remaining rows = metabolites
data <- apply(data, 2, as.numeric)            # convert to numeric
rownames(data) <- rownames(batch2_positive_susceptible_120hpi)[-1]

# 3. Compute log2 fold change

mean_R <- rowMeans(data[, genotype == "E2-AP"])
mean_S <- rowMeans(data[, genotype == "E2-C"])

log2FC <- log2(mean_R / mean_S)   


# 4. Compute p-values (t-test)

pvalues <- apply(data, 1, function(x) {
  t.test(x[genotype == "E2-AP"], x[genotype == "E2-C"])$p.value
})


# 5. Adjust p-values (FDR)

padj <- p.adjust(pvalues, method = "fdr")


# 6. Create results_batch2_positive_susceptible_120hpi table

results_batch2_positive_susceptible_120hpi <- data.frame(
  metabolite = rownames(data),
  log2FC = log2FC,
  pvalue = pvalues,
  padj = padj
)

results_batch2_positive_susceptible_120hpi$negLog10P <- -log10(results_batch2_positive_susceptible_120hpi$padj)


# 7. Volcano plot

# Define significance
results_batch2_positive_susceptible_120hpi$significant <- "Not Significant"
results_batch2_positive_susceptible_120hpi$significant[results_batch2_positive_susceptible_120hpi$log2FC > 1 & results_batch2_positive_susceptible_120hpi$padj < 0.05] <- "Up"
results_batch2_positive_susceptible_120hpi$significant[results_batch2_positive_susceptible_120hpi$log2FC < -1 & results_batch2_positive_susceptible_120hpi$padj < 0.05] <- "Down"

Volcano_batch2_positive_susceptible_120hpi <- ggplot(results_batch2_positive_susceptible_120hpi, aes(x = log2FC, y = negLog10P, color = significant)) +
  geom_point(size = 3, alpha = 0.8) +
  scale_color_manual(values = c("blue", "grey", "red")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  theme_minimal() +
  labs(
    x = "Log2 Fold Change (E2-AP / E2-C)",
    y = "-Log10(FDR)",
    title = "Volcano Plot"
  ) +
  geom_text_repel(
    data = subset(results_batch2_positive_susceptible_120hpi, significant != "Not Significant"),
    aes(label = metabolite)
  )+
  # Set fixed x and y axis limits
  scale_x_continuous(limits = c(-5, 5)) +  # Adjust min/max for log2FC
  scale_y_continuous(limits = c(0, 10))    # Adjust min/max for -log10(FDR)

results_batch2_positive_susceptible_120hpi$significant<-as.factor(results_batch2_positive_susceptible_120hpi$significant)
summary(results_batch2_positive_susceptible_120hpi)







##batch2_positive_resistant----
#batch2_positive_resistant----
# Keep only rows where genotype is NOT "E2-S"
batch2_positive_resistant <- batch2_positive[batch2_positive$genotype != "E2-S", ]
batch2_positive_resistant <- batch2_positive_resistant[, !(names(batch2_positive_resistant) %in% c("FEATURE_ID", "genotype", "hpi", "Sample_ID", "Injection_order", "Group", "QC"))]
batch2_positive_resistant <- t(batch2_positive_resistant)



# 2. Extract genotype and data

genotype <- as.character(unlist(batch2_positive_resistant[1, ])) # first row = genotype
data <- batch2_positive_resistant[-1, ]                          # remaining rows = metabolites
data <- apply(data, 2, as.numeric)            # convert to numeric
rownames(data) <- rownames(batch2_positive_resistant)[-1]

# 3. Compute log2 fold change

mean_R <- rowMeans(data[, genotype == "E2-AP"])
mean_S <- rowMeans(data[, genotype == "E2-C"])

log2FC <- log2(mean_R / mean_S)   


# 4. Compute p-values (t-test)

pvalues <- apply(data, 1, function(x) {
  t.test(x[genotype == "E2-AP"], x[genotype == "E2-C"])$p.value
})


# 5. Adjust p-values (FDR)

padj <- p.adjust(pvalues, method = "fdr")


# 6. Create results_batch2_positive_resistant table

results_batch2_positive_resistant <- data.frame(
  metabolite = rownames(data),
  log2FC = log2FC,
  pvalue = pvalues,
  padj = padj
)

results_batch2_positive_resistant$negLog10P <- -log10(results_batch2_positive_resistant$padj)


# 7. Volcano plot

# Define significance
results_batch2_positive_resistant$significant <- "Not Significant"
results_batch2_positive_resistant$significant[results_batch2_positive_resistant$log2FC > 1 & results_batch2_positive_resistant$padj < 0.05] <- "Up"
results_batch2_positive_resistant$significant[results_batch2_positive_resistant$log2FC < -1 & results_batch2_positive_resistant$padj < 0.05] <- "Down"

Volcano_batch2_positive_resistant <- ggplot(results_batch2_positive_resistant, aes(x = log2FC, y = negLog10P, color = significant)) +
  geom_point(size = 3, alpha = 0.8) +
  scale_color_manual(values = c("blue", "grey", "red")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  theme_minimal() +
  labs(
    x = "Log2 Fold Change (E2-AP / E2-C)",
    y = "-Log10(FDR)",
    title = "Volcano Plot"
  ) +
  geom_text_repel(
    data = subset(results_batch2_positive_resistant, significant != "Not Significant"),
    aes(label = metabolite)
  )+
  # Set fixed x and y axis limits
  scale_x_continuous(limits = c(-5, 5)) +  # Adjust min/max for log2FC
  scale_y_continuous(limits = c(0, 10))    # Adjust min/max for -log10(FDR)

results_batch2_positive_resistant$significant<-as.factor(results_batch2_positive_resistant$significant)
summary(results_batch2_positive_resistant)




#early----
# Keep only rows where genotype is NOT "E2-S"
batch2_positive_resistant <- batch2_positive[batch2_positive$genotype != "E2-S", ]
batch2_positive_resistant_early <- batch2_positive_resistant[
  batch2_positive_resistant$hpi %in% c("3", "24"),]
batch2_positive_resistant_early <- batch2_positive_resistant_early[, !(names(batch2_positive_resistant_early) %in% c("FEATURE_ID", "genotype", "hpi", "Sample_ID", "Injection_order", "Group", "QC"))]
batch2_positive_resistant_early <- t(batch2_positive_resistant_early)



# 2. Extract genotype and data

genotype <- as.character(unlist(batch2_positive_resistant_early[1, ])) # first row = genotype
data <- batch2_positive_resistant_early[-1, ]                          # remaining rows = metabolites
data <- apply(data, 2, as.numeric)            # convert to numeric
rownames(data) <- rownames(batch2_positive_resistant_early)[-1]

# 3. Compute log2 fold change

mean_R <- rowMeans(data[, genotype == "E2-AP"])
mean_S <- rowMeans(data[, genotype == "E2-C"])

log2FC <- log2(mean_R / mean_S)   


# 4. Compute p-values (t-test)

pvalues <- apply(data, 1, function(x) {
  t.test(x[genotype == "E2-AP"], x[genotype == "E2-C"])$p.value
})


# 5. Adjust p-values (FDR)

padj <- p.adjust(pvalues, method = "fdr")


# 6. Create results_batch2_positive_resistant_early table

results_batch2_positive_resistant_early <- data.frame(
  metabolite = rownames(data),
  log2FC = log2FC,
  pvalue = pvalues,
  padj = padj
)

results_batch2_positive_resistant_early$negLog10P <- -log10(results_batch2_positive_resistant_early$padj)


# 7. Volcano plot

# Define significance
results_batch2_positive_resistant_early$significant <- "Not Significant"
results_batch2_positive_resistant_early$significant[results_batch2_positive_resistant_early$log2FC > 1 & results_batch2_positive_resistant_early$padj < 0.05] <- "Up"
results_batch2_positive_resistant_early$significant[results_batch2_positive_resistant_early$log2FC < -1 & results_batch2_positive_resistant_early$padj < 0.05] <- "Down"

Volcano_batch2_positive_resistant_early <- ggplot(results_batch2_positive_resistant_early, aes(x = log2FC, y = negLog10P, color = significant)) +
  geom_point(size = 3, alpha = 0.8) +
  scale_color_manual(values = c("blue", "grey", "red")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  theme_minimal() +
  labs(
    x = "Log2 Fold Change (E2-AP / E2-C)",
    y = "-Log10(FDR)",
    title = "Volcano Plot"
  ) +
  geom_text_repel(
    data = subset(results_batch2_positive_resistant_early, significant != "Not Significant"),
    aes(label = metabolite)
  )+
  # Set fixed x and y axis limits
  scale_x_continuous(limits = c(-5, 5)) +  # Adjust min/max for log2FC
  scale_y_continuous(limits = c(0, 10))    # Adjust min/max for -log10(FDR)

results_batch2_positive_resistant_early$significant<-as.factor(results_batch2_positive_resistant_early$significant)
summary(results_batch2_positive_resistant_early)





#late----
# Keep only rows where genotype is NOT "E2-S"
batch2_positive_resistant <- batch2_positive[batch2_positive$genotype != "E2-S", ]
batch2_positive_resistant_late <- batch2_positive_resistant[
  batch2_positive_resistant$hpi %in% c("72", "120"),
]
batch2_positive_resistant_late <- batch2_positive_resistant_late[, !(names(batch2_positive_resistant_late) %in% c("FEATURE_ID", "genotype", "hpi", "Sample_ID", "Injection_order", "Group", "QC"))]
batch2_positive_resistant_late <- t(batch2_positive_resistant_late)



# 2. Extract genotype and data

genotype <- as.character(unlist(batch2_positive_resistant_late[1, ])) # first row = genotype
data <- batch2_positive_resistant_late[-1, ]                          # remaining rows = metabolites
data <- apply(data, 2, as.numeric)            # convert to numeric
rownames(data) <- rownames(batch2_positive_resistant_late)[-1]

# 3. Compute log2 fold change

mean_R <- rowMeans(data[, genotype == "E2-AP"])
mean_S <- rowMeans(data[, genotype == "E2-C"])

log2FC <- log2(mean_R / mean_S)   


# 4. Compute p-values (t-test)

pvalues <- apply(data, 1, function(x) {
  t.test(x[genotype == "E2-AP"], x[genotype == "E2-C"])$p.value
})


# 5. Adjust p-values (FDR)

padj <- p.adjust(pvalues, method = "fdr")


# 6. Create results_batch2_positive_resistant_late table

results_batch2_positive_resistant_late <- data.frame(
  metabolite = rownames(data),
  log2FC = log2FC,
  pvalue = pvalues,
  padj = padj
)

results_batch2_positive_resistant_late$negLog10P <- -log10(results_batch2_positive_resistant_late$padj)


# 7. Volcano plot

# Define significance
results_batch2_positive_resistant_late$significant <- "Not Significant"
results_batch2_positive_resistant_late$significant[results_batch2_positive_resistant_late$log2FC > 1 & results_batch2_positive_resistant_late$padj < 0.05] <- "Up"
results_batch2_positive_resistant_late$significant[results_batch2_positive_resistant_late$log2FC < -1 & results_batch2_positive_resistant_late$padj < 0.05] <- "Down"

Volcano_batch2_positive_resistant_late <- ggplot(results_batch2_positive_resistant_late, aes(x = log2FC, y = negLog10P, color = significant)) +
  geom_point(size = 3, alpha = 0.8) +
  scale_color_manual(values = c("blue", "grey", "red")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  theme_minimal() +
  labs(
    x = "Log2 Fold Change (E2-AP / E2-C)",
    y = "-Log10(FDR)",
    title = "Volcano Plot"
  ) +
  geom_text_repel(
    data = subset(results_batch2_positive_resistant_late, significant != "Not Significant"),
    aes(label = metabolite)
  )+
  # Set fixed x and y axis limits
  scale_x_continuous(limits = c(-5, 5)) +  # Adjust min/max for log2FC
  scale_y_continuous(limits = c(0, 10))    # Adjust min/max for -log10(FDR)

results_batch2_positive_resistant_late$significant<-as.factor(results_batch2_positive_resistant_late$significant)
summary(results_batch2_positive_resistant_late)




#3hpi----
# Keep only rows where genotype is NOT "E2-S"
batch2_positive_resistant <- batch2_positive[batch2_positive$genotype != "E2-S", ]
batch2_positive_resistant_3hpi <- batch2_positive_resistant[batch2_positive_resistant$hpi == "3", ]
batch2_positive_resistant_3hpi <- batch2_positive_resistant_3hpi[, !(names(batch2_positive_resistant_3hpi) %in% c("FEATURE_ID", "genotype", "hpi", "Sample_ID", "Injection_order", "Group", "QC"))]
batch2_positive_resistant_3hpi <- t(batch2_positive_resistant_3hpi)



# 2. Extract genotype and data

genotype <- as.character(unlist(batch2_positive_resistant_3hpi[1, ])) # first row = genotype
data <- batch2_positive_resistant_3hpi[-1, ]                          # remaining rows = metabolites
data <- apply(data, 2, as.numeric)            # convert to numeric
rownames(data) <- rownames(batch2_positive_resistant_3hpi)[-1]

# 3. Compute log2 fold change

mean_R <- rowMeans(data[, genotype == "E2-AP"])
mean_S <- rowMeans(data[, genotype == "E2-C"])

log2FC <- log2(mean_R / mean_S)   


# 4. Compute p-values (t-test)

pvalues <- apply(data, 1, function(x) {
  t.test(x[genotype == "E2-AP"], x[genotype == "E2-C"])$p.value
})


# 5. Adjust p-values (FDR)

padj <- p.adjust(pvalues, method = "fdr")


# 6. Create results_batch2_positive_resistant_3hpi table

results_batch2_positive_resistant_3hpi <- data.frame(
  metabolite = rownames(data),
  log2FC = log2FC,
  pvalue = pvalues,
  padj = padj
)

results_batch2_positive_resistant_3hpi$negLog10P <- -log10(results_batch2_positive_resistant_3hpi$padj)


# 7. Volcano plot

# Define significance
results_batch2_positive_resistant_3hpi$significant <- "Not Significant"
results_batch2_positive_resistant_3hpi$significant[results_batch2_positive_resistant_3hpi$log2FC > 1 & results_batch2_positive_resistant_3hpi$padj < 0.05] <- "Up"
results_batch2_positive_resistant_3hpi$significant[results_batch2_positive_resistant_3hpi$log2FC < -1 & results_batch2_positive_resistant_3hpi$padj < 0.05] <- "Down"

Volcano_batch2_positive_resistant_3hpi <- ggplot(results_batch2_positive_resistant_3hpi, aes(x = log2FC, y = negLog10P, color = significant)) +
  geom_point(size = 3, alpha = 0.8) +
  scale_color_manual(values = c("blue", "grey", "red")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  theme_minimal() +
  labs(
    x = "Log2 Fold Change (E2-AP / E2-C)",
    y = "-Log10(FDR)",
    title = "Volcano Plot"
  ) +
  geom_text_repel(
    data = subset(results_batch2_positive_resistant_3hpi, significant != "Not Significant"),
    aes(label = metabolite)
  )+
  # Set fixed x and y axis limits
  scale_x_continuous(limits = c(-5, 5)) +  # Adjust min/max for log2FC
  scale_y_continuous(limits = c(0, 10))    # Adjust min/max for -log10(FDR)

results_batch2_positive_resistant_3hpi$significant<-as.factor(results_batch2_positive_resistant_3hpi$significant)
summary(results_batch2_positive_resistant_3hpi)




#24hpi----
# Keep only rows where genotype is NOT "E2-S"
batch2_positive_resistant <- batch2_positive[batch2_positive$genotype != "E2-S", ]
batch2_positive_resistant_24hpi <- batch2_positive_resistant[batch2_positive_resistant$hpi == "24", ]
batch2_positive_resistant_24hpi <- batch2_positive_resistant_24hpi[, !(names(batch2_positive_resistant_24hpi) %in% c("FEATURE_ID", "genotype", "hpi", "Sample_ID", "Injection_order", "Group", "QC"))]
batch2_positive_resistant_24hpi <- t(batch2_positive_resistant_24hpi)



# 2. Extract genotype and data

genotype <- as.character(unlist(batch2_positive_resistant_24hpi[1, ])) # first row = genotype
data <- batch2_positive_resistant_24hpi[-1, ]                          # remaining rows = metabolites
data <- apply(data, 2, as.numeric)            # convert to numeric
rownames(data) <- rownames(batch2_positive_resistant_24hpi)[-1]

# 3. Compute log2 fold change

mean_R <- rowMeans(data[, genotype == "E2-AP"])
mean_S <- rowMeans(data[, genotype == "E2-C"])

log2FC <- log2(mean_R / mean_S)   


# 4. Compute p-values (t-test)

pvalues <- apply(data, 1, function(x) {
  t.test(x[genotype == "E2-AP"], x[genotype == "E2-C"])$p.value
})


# 5. Adjust p-values (FDR)

padj <- p.adjust(pvalues, method = "fdr")


# 6. Create results_batch2_positive_resistant_24hpi table

results_batch2_positive_resistant_24hpi <- data.frame(
  metabolite = rownames(data),
  log2FC = log2FC,
  pvalue = pvalues,
  padj = padj
)

results_batch2_positive_resistant_24hpi$negLog10P <- -log10(results_batch2_positive_resistant_24hpi$padj)


# 7. Volcano plot

# Define significance
results_batch2_positive_resistant_24hpi$significant <- "Not Significant"
results_batch2_positive_resistant_24hpi$significant[results_batch2_positive_resistant_24hpi$log2FC > 1 & results_batch2_positive_resistant_24hpi$padj < 0.05] <- "Up"
results_batch2_positive_resistant_24hpi$significant[results_batch2_positive_resistant_24hpi$log2FC < -1 & results_batch2_positive_resistant_24hpi$padj < 0.05] <- "Down"

Volcano_batch2_positive_resistant_24hpi <- ggplot(results_batch2_positive_resistant_24hpi, aes(x = log2FC, y = negLog10P, color = significant)) +
  geom_point(size = 3, alpha = 0.8) +
  scale_color_manual(values = c("blue", "grey", "red")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  theme_minimal() +
  labs(
    x = "Log2 Fold Change (E2-AP / E2-C)",
    y = "-Log10(FDR)",
    title = "Volcano Plot"
  ) +
  geom_text_repel(
    data = subset(results_batch2_positive_resistant_24hpi, significant != "Not Significant"),
    aes(label = metabolite)
  )+
  # Set fixed x and y axis limits
  scale_x_continuous(limits = c(-5, 5)) +  # Adjust min/max for log2FC
  scale_y_continuous(limits = c(0, 10))    # Adjust min/max for -log10(FDR)

results_batch2_positive_resistant_24hpi$significant<-as.factor(results_batch2_positive_resistant_24hpi$significant)
summary(results_batch2_positive_resistant_24hpi)




#72hpi----
# Keep only rows where genotype is NOT "E2-S"
batch2_positive_resistant <- batch2_positive[batch2_positive$genotype != "E2-S", ]
batch2_positive_resistant_72hpi <- batch2_positive_resistant[batch2_positive_resistant$hpi == "72", ]
batch2_positive_resistant_72hpi <- batch2_positive_resistant_72hpi[, !(names(batch2_positive_resistant_72hpi) %in% c("FEATURE_ID", "genotype", "hpi", "Sample_ID", "Injection_order", "Group", "QC"))]
batch2_positive_resistant_72hpi <- t(batch2_positive_resistant_72hpi)



# 2. Extract genotype and data

genotype <- as.character(unlist(batch2_positive_resistant_72hpi[1, ])) # first row = genotype
data <- batch2_positive_resistant_72hpi[-1, ]                          # remaining rows = metabolites
data <- apply(data, 2, as.numeric)            # convert to numeric
rownames(data) <- rownames(batch2_positive_resistant_72hpi)[-1]

# 3. Compute log2 fold change

mean_R <- rowMeans(data[, genotype == "E2-AP"])
mean_S <- rowMeans(data[, genotype == "E2-C"])

log2FC <- log2(mean_R / mean_S)   


# 4. Compute p-values (t-test)

pvalues <- apply(data, 1, function(x) {
  t.test(x[genotype == "E2-AP"], x[genotype == "E2-C"])$p.value
})


# 5. Adjust p-values (FDR)

padj <- p.adjust(pvalues, method = "fdr")


# 6. Create results_batch2_positive_resistant_72hpi table

results_batch2_positive_resistant_72hpi <- data.frame(
  metabolite = rownames(data),
  log2FC = log2FC,
  pvalue = pvalues,
  padj = padj
)

results_batch2_positive_resistant_72hpi$negLog10P <- -log10(results_batch2_positive_resistant_72hpi$padj)


# 7. Volcano plot

# Define significance
results_batch2_positive_resistant_72hpi$significant <- "Not Significant"
results_batch2_positive_resistant_72hpi$significant[results_batch2_positive_resistant_72hpi$log2FC > 1 & results_batch2_positive_resistant_72hpi$padj < 0.05] <- "Up"
results_batch2_positive_resistant_72hpi$significant[results_batch2_positive_resistant_72hpi$log2FC < -1 & results_batch2_positive_resistant_72hpi$padj < 0.05] <- "Down"

Volcano_batch2_positive_resistant_72hpi <- ggplot(results_batch2_positive_resistant_72hpi, aes(x = log2FC, y = negLog10P, color = significant)) +
  geom_point(size = 3, alpha = 0.8) +
  scale_color_manual(values = c("blue", "grey", "red")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  theme_minimal() +
  labs(
    x = "Log2 Fold Change (E2-AP / E2-C)",
    y = "-Log10(FDR)",
    title = "Volcano Plot"
  ) +
  geom_text_repel(
    data = subset(results_batch2_positive_resistant_72hpi, significant != "Not Significant"),
    aes(label = metabolite)
  )+
  # Set fixed x and y axis limits
  scale_x_continuous(limits = c(-5, 5)) +  # Adjust min/max for log2FC
  scale_y_continuous(limits = c(0, 10))    # Adjust min/max for -log10(FDR)

results_batch2_positive_resistant_72hpi$significant<-as.factor(results_batch2_positive_resistant_72hpi$significant)
summary(results_batch2_positive_resistant_72hpi)

#120hpi----
# Keep only rows where genotype is NOT "E2-S"
batch2_positive_resistant <- batch2_positive[batch2_positive$genotype != "E2-S", ]
batch2_positive_resistant_120hpi <- batch2_positive_resistant[batch2_positive_resistant$hpi == "120", ]
batch2_positive_resistant_120hpi <- batch2_positive_resistant_120hpi[, !(names(batch2_positive_resistant_120hpi) %in% c("FEATURE_ID", "genotype", "hpi", "Sample_ID", "Injection_order", "Group", "QC"))]
batch2_positive_resistant_120hpi <- t(batch2_positive_resistant_120hpi)



# 2. Extract genotype and data

genotype <- as.character(unlist(batch2_positive_resistant_120hpi[1, ])) # first row = genotype
data <- batch2_positive_resistant_120hpi[-1, ]                          # remaining rows = metabolites
data <- apply(data, 2, as.numeric)            # convert to numeric
rownames(data) <- rownames(batch2_positive_resistant_120hpi)[-1]

# 3. Compute log2 fold change

mean_R <- rowMeans(data[, genotype == "E2-AP"])
mean_S <- rowMeans(data[, genotype == "E2-C"])

log2FC <- log2(mean_R / mean_S)   


# 4. Compute p-values (t-test)

pvalues <- apply(data, 1, function(x) {
  t.test(x[genotype == "E2-AP"], x[genotype == "E2-C"])$p.value
})


# 5. Adjust p-values (FDR)

padj <- p.adjust(pvalues, method = "fdr")


# 6. Create results_batch2_positive_resistant_120hpi table

results_batch2_positive_resistant_120hpi <- data.frame(
  metabolite = rownames(data),
  log2FC = log2FC,
  pvalue = pvalues,
  padj = padj
)

results_batch2_positive_resistant_120hpi$negLog10P <- -log10(results_batch2_positive_resistant_120hpi$padj)


# 7. Volcano plot

# Define significance
results_batch2_positive_resistant_120hpi$significant <- "Not Significant"
results_batch2_positive_resistant_120hpi$significant[results_batch2_positive_resistant_120hpi$log2FC > 1 & results_batch2_positive_resistant_120hpi$padj < 0.05] <- "Up"
results_batch2_positive_resistant_120hpi$significant[results_batch2_positive_resistant_120hpi$log2FC < -1 & results_batch2_positive_resistant_120hpi$padj < 0.05] <- "Down"

Volcano_batch2_positive_resistant_120hpi <- ggplot(results_batch2_positive_resistant_120hpi, aes(x = log2FC, y = negLog10P, color = significant)) +
  geom_point(size = 3, alpha = 0.8) +
  scale_color_manual(values = c("blue", "grey", "red")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  theme_minimal() +
  labs(
    x = "Log2 Fold Change (E2-AP / E2-C)",
    y = "-Log10(FDR)",
    title = "Volcano Plot"
  ) +
  geom_text_repel(
    data = subset(results_batch2_positive_resistant_120hpi, significant != "Not Significant"),
    aes(label = metabolite)
  )+
  # Set fixed x and y axis limits
  scale_x_continuous(limits = c(-5, 5)) +  # Adjust min/max for log2FC
  scale_y_continuous(limits = c(0, 10))    # Adjust min/max for -log10(FDR)

results_batch2_positive_resistant_120hpi$significant<-as.factor(results_batch2_positive_resistant_120hpi$significant)
summary(results_batch2_positive_resistant_120hpi)





#####

results_batch2_positive_susceptible<-results_batch2_positive_susceptible[,c(1,6)]
colnames(results_batch2_positive_susceptible)[colnames(results_batch2_positive_susceptible) == "significant"] <- "batch2_positive_susceptible"
results_batch2_positive_resistant<-results_batch2_positive_resistant[,c(1,6)]
colnames(results_batch2_positive_resistant)[colnames(results_batch2_positive_resistant) == "significant"] <- "batch2_positive_resistant"

results_batch2_positive_susceptible_early<-results_batch2_positive_susceptible_early[,c(1,6)]
colnames(results_batch2_positive_susceptible_early)[colnames(results_batch2_positive_susceptible_early) == "significant"] <- "batch2_positive_susceptible_early"

results_batch2_positive_resistant_early<-results_batch2_positive_resistant_early[,c(1,6)]
colnames(results_batch2_positive_resistant_early)[colnames(results_batch2_positive_resistant_early) == "significant"] <- "batch2_positive_resistant_early"


results_batch2_positive_susceptible_late<-results_batch2_positive_susceptible_late[,c(1,6)]
colnames(results_batch2_positive_susceptible_late)[colnames(results_batch2_positive_susceptible_late) == "significant"] <- "batch2_positive_susceptible_late"
results_batch2_positive_resistant_late<-results_batch2_positive_resistant_late[,c(1,6)]
colnames(results_batch2_positive_resistant_late)[colnames(results_batch2_positive_resistant_late) == "significant"] <- "batch2_positive_resistant_late"
results_batch2_positive_susceptible_3hpi<-results_batch2_positive_susceptible_3hpi[,c(1,6)]
colnames(results_batch2_positive_susceptible_3hpi)[colnames(results_batch2_positive_susceptible_3hpi) == "significant"] <- "batch2_positive_susceptible_3hpi"
results_batch2_positive_resistant_3hpi<-results_batch2_positive_resistant_3hpi[,c(1,6)]
colnames(results_batch2_positive_resistant_3hpi)[colnames(results_batch2_positive_resistant_3hpi) == "significant"] <- "batch2_positive_resistant_3hpi"
results_batch2_positive_susceptible_24hpi<-results_batch2_positive_susceptible_24hpi[,c(1,6)]
colnames(results_batch2_positive_susceptible_24hpi)[colnames(results_batch2_positive_susceptible_24hpi) == "significant"] <- "batch2_positive_susceptible_24hpi"
results_batch2_positive_resistant_24hpi<-results_batch2_positive_resistant_24hpi[,c(1,6)]
colnames(results_batch2_positive_resistant_24hpi)[colnames(results_batch2_positive_resistant_24hpi) == "significant"] <- "batch2_positive_resistant_24hpi"
results_batch2_positive_susceptible_72hpi<-results_batch2_positive_susceptible_72hpi[,c(1,6)]
colnames(results_batch2_positive_susceptible_72hpi)[colnames(results_batch2_positive_susceptible_72hpi) == "significant"] <- "batch2_positive_susceptible_72hpi"
results_batch2_positive_resistant_72hpi<-results_batch2_positive_resistant_72hpi[,c(1,6)]
colnames(results_batch2_positive_resistant_72hpi)[colnames(results_batch2_positive_resistant_72hpi) == "significant"] <- "batch2_positive_resistant_72hpi"
results_batch2_positive_susceptible_120hpi<-results_batch2_positive_susceptible_120hpi[,c(1,6)]
colnames(results_batch2_positive_susceptible_120hpi)[colnames(results_batch2_positive_susceptible_120hpi) == "significant"] <- "batch2_positive_susceptible_120hpi"
results_batch2_positive_resistant_120hpi<-results_batch2_positive_resistant_120hpi[,c(1,6)]
colnames(results_batch2_positive_resistant_120hpi)[colnames(results_batch2_positive_resistant_120hpi) == "significant"] <- "batch2_positive_resistant_120hpi"




all_results <- list(
  results_batch2_positive_susceptible,
  results_batch2_positive_resistant,
  results_batch2_positive_susceptible_early,
  results_batch2_positive_resistant_early,
  results_batch2_positive_susceptible_late,
  results_batch2_positive_resistant_late,
  results_batch2_positive_susceptible_3hpi,
  results_batch2_positive_resistant_3hpi,
  results_batch2_positive_susceptible_24hpi,
  results_batch2_positive_resistant_24hpi,
  results_batch2_positive_susceptible_72hpi,
  results_batch2_positive_resistant_72hpi,
  results_batch2_positive_susceptible_120hpi,
  results_batch2_positive_resistant_120hpi
)

merged_results_results_batch2_positive <- reduce(all_results, full_join, by = "metabolite")
merged_results_results_batch2_positive <- merged_results_results_batch2_positive[
  rowSums(merged_results_results_batch2_positive[, -1] != "Not Significant", na.rm = TRUE) > 0,
]

merged_results_results_batch2_positive$metabolite <- paste0(
  merged_results_results_batch2_positive$metabolite,
  "pos")

heatmap_batch2_positive <- merged_results_results_batch2_positive


heatmap_batch2_positive <- as.data.frame(lapply(heatmap_batch2_positive, as.character))
rownames(heatmap_batch2_positive) <- heatmap_batch2_positive$metabolite
heatmap_batch2_positive <- heatmap_batch2_positive[, -1]

heatmap_batch2_positive[] <- lapply(heatmap_batch2_positive, function(x) {
  x[x == "Up"] <- 1
  x[x == "Down"] <- -1
  x[x == "Not Significant"] <- 0
  return(as.numeric(x))
})

# Convert to matrix
heatmap_matrix <- as.matrix(heatmap_batch2_positive)

# Convert to numeric
mode(heatmap_matrix) <- "numeric"

# compute row sums (total signal per metabolite)
row_score <- rowSums(heatmap_matrix, na.rm = TRUE)

# order from highest to lowest
heatmap_matrix_batch2_positive <- heatmap_matrix[order(row_score, decreasing = TRUE), ]


# Plot heatmap
pheatmap(
  heatmap_matrix_batch2_positive,
  color = c("blue", "white", "red"),
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  border_color = "grey80"
)


#join positive and negative#####

fix_names <- function(x) {
  parts <- strsplit(x, "_")
  
  sapply(parts, function(p) {
    keep <- p[c(1, 3:min(length(p), 4))]  # keep 1st, 3rd, 4th if exists
    paste(keep, collapse = "_")
  })
}
colnames(heatmap_matrix_batch1_negative) <- fix_names(colnames(heatmap_matrix_batch1_negative))
colnames(heatmap_matrix_batch1_positive) <- fix_names(colnames(heatmap_matrix_batch1_positive))

colnames(heatmap_matrix_batch2_negative) <- fix_names(colnames(heatmap_matrix_batch2_negative))
colnames(heatmap_matrix_batch2_positive) <- fix_names(colnames(heatmap_matrix_batch2_positive))

batch1_combined_matrix <- rbind(heatmap_matrix_batch1_negative, heatmap_matrix_batch1_positive)
batch2_combined_matrix <- rbind(heatmap_matrix_batch2_negative, heatmap_matrix_batch2_positive)



# compute row sums (total signal per metabolite)
row_score <- rowSums(batch1_combined_matrix, na.rm = TRUE)

# order from highest to lowest
batch1_combined_matrix <- batch1_combined_matrix[order(row_score, decreasing = TRUE), ]


# Plot heatmap
pheatmap(
  batch1_combined_matrix,
  color = c("blue", "white", "red"),
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  border_color = "grey80"
)

# compute row sums (total signal per metabolite)
row_score <- rowSums(batch2_combined_matrix, na.rm = TRUE)

# order from highest to lowest
batch2_combined_matrix <- batch1_combined_matrix[order(row_score, decreasing = TRUE), ]


# Plot heatmap
pheatmap(
  batch2_combined_matrix,
  color = c("blue", "white", "red"),
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  border_color = "grey80"
)


#####

colnames(merged_results_results_batch1_negative) <- fix_names(colnames(merged_results_results_batch1_negative))
colnames(merged_results_results_batch2_negative) <- fix_names(colnames(merged_results_results_batch2_negative))
colnames(merged_results_results_batch1_positive) <- fix_names(colnames(merged_results_results_batch1_positive))
colnames(merged_results_results_batch2_positive) <- fix_names(colnames(merged_results_results_batch2_positive))

merged_results_results_batch1 <- rbind(merged_results_results_batch1_negative, merged_results_results_batch1_positive)
colnames(merged_results_results_batch1)[colnames(merged_results_results_batch1) == "metabolite_NA_NA_metabolite"] <- "metabolite"

merged_results_results_batch2 <- rbind(merged_results_results_batch2_positive, merged_results_results_batch2_negative)
colnames(merged_results_results_batch2)[colnames(merged_results_results_batch2) == "metabolite_NA_NA_metabolite"] <- "metabolite"






ID_positive1 <- read_tsv('../Positive mode/metabolites_identification/SIRIUS_positive_mode/canopus_formula_summary.tsv')
colnames(ID_positive1)[colnames(ID_positive1) == "mappingFeatureId"] <- "metabolite"
ID_positive1<-ID_positive1[,c(5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,28)]
ID_positive2 <- read_tsv('../Positive mode/metabolites_identification/SIRIUS_positive_mode/denovo_structure_identifications.tsv')
colnames(ID_positive2)[colnames(ID_positive2) == "mappingFeatureId"] <- "metabolite"
ID_positive2<-ID_positive2[,c(13,21)]
ID_positive3 <- read_tsv('../Positive mode/metabolites_identification/GNPS_Positive_mode/DB_result/GNPS_positive.tsv')
colnames(ID_positive3)[colnames(ID_positive3) == "#Scan#"] <- "metabolite"
ID_positive3<-ID_positive3[,c(2,21)]

merged_ID_positive_pre <- merge(ID_positive1,ID_positive2, by ="metabolite",all=TRUE)
merged_ID_positive <- merge(merged_ID_positive_pre,ID_positive3, by ="metabolite",all=TRUE)

merged_ID_positive$metabolite <- str_c(
  merged_ID_positive$metabolite,
  "pos"
)

ID_negative1 <- read_tsv('../Negative mode/metabolites identification/SIRIUS_Negative_mode/canopus_formula_summary.tsv')
colnames(ID_negative1)[colnames(ID_negative1) == "mappingFeatureId"] <- "metabolite"
ID_negative1<-ID_negative1[,c(5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,28)]
ID_negative2 <- read_tsv('../Negative mode/metabolites identification/SIRIUS_Negative_mode/denovo_structure_identifications.tsv')
colnames(ID_negative2)[colnames(ID_negative2) == "mappingFeatureId"] <- "metabolite"
ID_negative2<-ID_negative2[,c(13,21)]
ID_negative3 <- read_tsv('../Negative mode/metabolites identification/GNPS_Negative mode/DB_result/GNPS_negative.tsv')
colnames(ID_negative3)[colnames(ID_negative3) == "#Scan#"] <- "metabolite"
ID_negative3<-ID_negative3[,c(2,21)]

merged_ID_negative_pre <- merge(ID_negative1,ID_negative2, by ="metabolite",all=TRUE)
merged_ID_negative <- merge(merged_ID_negative_pre,ID_negative3, by ="metabolite",all=TRUE)

merged_ID_negative$metabolite <- str_c(
  merged_ID_negative$metabolite,
  "neg"
)

merged_ID <- rbind(merged_ID_negative, merged_ID_positive)


merged_results_batch1_ID <- merge(merged_results_results_batch1, merged_ID, by ="metabolite",all=FALSE)
merged_results_batch2_ID <- merge(merged_results_results_batch2, merged_ID, by ="metabolite",all=FALSE)






