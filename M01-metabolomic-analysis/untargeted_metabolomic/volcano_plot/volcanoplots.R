###packages----
library(dplyr)
library(pheatmap)
library(ggplot2)
library(ggrepel)
library(readr)
library(readxl)
library(tidyr)
library(ggvenn)
#----
setwd("~/Documents/metabolomic_too_heavy_for_git/untargeted/Negative mode/")
##batch_1_negative----

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

#batch1_negative_constitutive----
# Keep only rows where treatment is NOT "E1-AF"
batch1_negative_constitutive <- batch1_negative[batch1_negative$treatment != "E1-AF", ]
batch1_negative_constitutive <- batch1_negative_constitutive[, !(names(batch1_negative_constitutive) %in% c("FEATURE_ID", "treatment", "hpi", "Sample_ID", "Injection_order", "Group", "QC"))]
batch1_negative_constitutive <- t(batch1_negative_constitutive)



# 2. Extract genotype and data

genotype <- as.character(unlist(batch1_negative_constitutive[1, ])) # first row = genotype
data <- batch1_negative_constitutive[-1, ]                          # remaining rows = metabolites
data <- apply(data, 2, as.numeric)            # convert to numeric
rownames(data) <- rownames(batch1_negative_constitutive)[-1]

# 3. Compute log2 fold change

mean_R <- rowMeans(data[, genotype == "E1-R"])
mean_S <- rowMeans(data[, genotype == "E1-S"])

log2FC <- log2(mean_R / mean_S)   # E1-R / E1-S


# 4. Compute p-values (t-test)

pvalues <- apply(data, 1, function(x) {
  t.test(x[genotype == "E1-R"], x[genotype == "E1-S"])$p.value
})

# 5. Adjust p-values (FDR)

padj <- p.adjust(pvalues, method = "fdr")


# 6. Create results_batch1_negative_constitutive table

results_batch1_negative_constitutive <- data.frame(
  metabolite = rownames(data),
  log2FC = log2FC,
  pvalue = pvalues,
  padj = padj
)

results_batch1_negative_constitutive$negLog10P <- -log10(results_batch1_negative_constitutive$padj)


# 7. Volcano plot

# Define significance
results_batch1_negative_constitutive$significant <- "Not Significant"
results_batch1_negative_constitutive$significant[results_batch1_negative_constitutive$log2FC > 1 & results_batch1_negative_constitutive$padj < 0.05] <- "Up"
results_batch1_negative_constitutive$significant[results_batch1_negative_constitutive$log2FC < -1 & results_batch1_negative_constitutive$padj < 0.05] <- "Down"

Volcano_batch1_negative_constitutive <- ggplot(results_batch1_negative_constitutive, aes(x = log2FC, y = negLog10P, color = significant)) +
  geom_point(size = 3, alpha = 0.8) +
  scale_color_manual(values = c("blue", "grey", "red")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  theme_minimal() +
  labs(
    x = "Log2 Fold Change (E1-R / E1-S)",
    y = "-Log10(FDR)",
    title = "Volcano Plot"
  ) +
  geom_text_repel(
    data = subset(results_batch1_negative_constitutive, significant != "Not Significant"),
    aes(label = metabolite)
  )
results_batch1_negative_constitutive$significant<-as.factor(results_batch1_negative_constitutive$significant)
summary(results_batch1_negative_constitutive)



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


##batch_2_negative----

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


#batch2_negative_constitutive----
# Keep only rows where treatment is NOT "E2-AP"
batch2_negative_constitutive <- batch2_negative[batch2_negative$treatment != "E2-AP", ]
batch2_negative_constitutive <- batch2_negative_constitutive[, !(names(batch2_negative_constitutive) %in% c("FEATURE_ID", "treatment", "hpi", "Sample_ID", "Injection_order", "Group", "QC"))]
batch2_negative_constitutive <- t(batch2_negative_constitutive)




# 2. Extract genotype and data

genotype <- as.character(unlist(batch2_negative_constitutive[1, ])) # first row = genotype
data <- batch2_negative_constitutive[-1, ]                          # remaining rows = metabolites
data <- apply(data, 2, as.numeric)            # convert to numeric
rownames(data) <- rownames(batch2_negative_constitutive)[-1]


# 3. Compute log2 fold change

mean_R <- rowMeans(data[, genotype == "E2-R"])
mean_S <- rowMeans(data[, genotype == "E2-S"])

log2FC <- log2(mean_R / mean_S)   # E2-R / E2-S


# 4. Compute p-values (t-test)

pvalues <- apply(data, 1, function(x) {
  t.test(x[genotype == "E2-R"], x[genotype == "E2-S"])$p.value
})

# 5. Adjust p-values (FDR)

padj <- p.adjust(pvalues, method = "fdr")


# 6. Create results_batch2_negative_constitutive table

results_batch2_negative_constitutive <- data.frame(
  metabolite = rownames(data),
  log2FC = log2FC,
  pvalue = pvalues,
  padj = padj
)

results_batch2_negative_constitutive$negLog10P <- -log10(results_batch2_negative_constitutive$padj)


# 7. Volcano plot

# Define significance
results_batch2_negative_constitutive$significant <- "Not Significant"
results_batch2_negative_constitutive$significant[results_batch2_negative_constitutive$log2FC > 1 & results_batch2_negative_constitutive$padj < 0.05] <- "Up"
results_batch2_negative_constitutive$significant[results_batch2_negative_constitutive$log2FC < -1 & results_batch2_negative_constitutive$padj < 0.05] <- "Down"

Volcano_batch2_negative_constitutive <- ggplot(results_batch2_negative_constitutive, aes(x = log2FC, y = negLog10P, color = significant)) +
  geom_point(size = 3, alpha = 0.8) +
  scale_color_manual(values = c("blue", "grey", "red")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  theme_minimal() +
  labs(
    x = "Log2 Fold Change (E2-R / E2-S)",
    y = "-Log10(FDR)",
    title = "Volcano Plot"
  ) +
  geom_text_repel(
    data = subset(results_batch2_negative_constitutive, significant != "Not Significant"),
    aes(label = metabolite)
  )
results_batch2_negative_constitutive$significant<-as.factor(results_batch2_negative_constitutive$significant)
summary(results_batch2_negative_constitutive)



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
    data = subset(results_batch2_negative_susceptible, significant != "Not Significant"),
    aes(label = metabolite)
  )+
  # Set fixed x and y axis limits
  scale_x_continuous(limits = c(-2, 4)) +  # Adjust min/max for log2FC
  scale_y_continuous(limits = c(0, 4))    # Adjust min/max for -log10(FDR)

results_batch2_negative_susceptible$significant<-as.factor(results_batch2_negative_susceptible$significant)
summary(results_batch2_negative_susceptible)



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
  scale_x_continuous(limits = c(-2, 4)) +  # Adjust min/max for log2FC
  scale_y_continuous(limits = c(0, 4))    # Adjust min/max for -log10(FDR)

results_batch2_negative_resistant$significant<-as.factor(results_batch2_negative_resistant$significant)
summary(results_batch2_negative_resistant)



#----   
setwd("~/Documents/metabolomic_too_heavy_for_git/untargeted/Positive mode")
##batch_1_positive----

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




#batch1_positive_constitutive----
# Keep only rows where treatment is NOT "E1-AF"
batch1_positive_constitutive <- batch1_positive[batch1_positive$treatment != "E1-AF", ]
batch1_positive_constitutive <- batch1_positive_constitutive[, !(names(batch1_positive_constitutive) %in% c("FEATURE_ID", "treatment", "hpi", "Sample_ID", "Injection_order", "Group", "QC"))]
batch1_positive_constitutive <- t(batch1_positive_constitutive)











# 2. Extract genotype and data

genotype <- as.character(unlist(batch1_positive_constitutive[1, ])) # first row = genotype
data <- batch1_positive_constitutive[-1, ]                          # remaining rows = metabolites
data <- apply(data, 2, as.numeric)            # convert to numeric
rownames(data) <- rownames(batch1_positive_constitutive)[-1]


# 3. Compute log2 fold change

mean_R <- rowMeans(data[, genotype == "E1-R"])
mean_S <- rowMeans(data[, genotype == "E1-S"])

log2FC <- log2(mean_R / mean_S)   # E1-R / E1-S


# 4. Compute p-values (t-test)

pvalues <- apply(data, 1, function(x) {
  t.test(x[genotype == "E1-R"], x[genotype == "E1-S"])$p.value
})

# 5. Adjust p-values (FDR)

padj <- p.adjust(pvalues, method = "fdr")


# 6. Create results_batch1_positive_constitutive table

results_batch1_positive_constitutive <- data.frame(
  metabolite = rownames(data),
  log2FC = log2FC,
  pvalue = pvalues,
  padj = padj
)

results_batch1_positive_constitutive$negLog10P <- -log10(results_batch1_positive_constitutive$padj)


# 7. Volcano plot

# Define significance
results_batch1_positive_constitutive$significant <- "Not Significant"
results_batch1_positive_constitutive$significant[results_batch1_positive_constitutive$log2FC > 1 & results_batch1_positive_constitutive$padj < 0.05] <- "Up"
results_batch1_positive_constitutive$significant[results_batch1_positive_constitutive$log2FC < -1 & results_batch1_positive_constitutive$padj < 0.05] <- "Down"

Volcano_batch1_positive_constitutive <- ggplot(results_batch1_positive_constitutive, aes(x = log2FC, y = negLog10P, color = significant)) +
  geom_point(size = 3, alpha = 0.8) +
  scale_color_manual(values = c("blue", "grey", "red")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  theme_minimal() +
  labs(
    x = "Log2 Fold Change (E1-R / E1-S)",
    y = "-Log10(FDR)",
    title = "Volcano Plot"
  ) +
  geom_text_repel(
    data = subset(results_batch1_positive_constitutive, significant != "Not Significant"),
    aes(label = metabolite)
  )
results_batch1_positive_constitutive$significant<-as.factor(results_batch1_positive_constitutive$significant)
summary(results_batch1_positive_constitutive)



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
  )
results_batch1_positive_susceptible$significant<-as.factor(results_batch1_positive_susceptible$significant)
summary(results_batch1_positive_susceptible)




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
  )
results_batch1_positive_resistant$significant<-as.factor(results_batch1_positive_resistant$significant)
summary(results_batch1_positive_resistant)



##batch_2_positive----

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

#batch2_positive_constitutive----
# Keep only rows where treatment is NOT "E2-AP"
batch2_positive_constitutive <- batch2_positive[batch2_positive$treatment != "E2-AP", ]
batch2_positive_constitutive <- batch2_positive_constitutive[, !(names(batch2_positive_constitutive) %in% c("FEATURE_ID", "treatment", "hpi", "Sample_ID", "Injection_order", "Group", "QC"))]
batch2_positive_constitutive <- t(batch2_positive_constitutive)




# 2. Extract genotype and data

genotype <- as.character(unlist(batch2_positive_constitutive[1, ])) # first row = genotype
data <- batch2_positive_constitutive[-1, ]                          # remaining rows = metabolites
data <- apply(data, 2, as.numeric)            # convert to numeric
rownames(data) <- rownames(batch2_positive_constitutive)[-1]


# 3. Compute log2 fold change

mean_R <- rowMeans(data[, genotype == "E2-R"])
mean_S <- rowMeans(data[, genotype == "E2-S"])

log2FC <- log2(mean_R / mean_S)   # E2-R / E2-S


# 4. Compute p-values (t-test)

pvalues <- apply(data, 1, function(x) {
  t.test(x[genotype == "E2-R"], x[genotype == "E2-S"])$p.value
})

# 5. Adjust p-values (FDR)

padj <- p.adjust(pvalues, method = "fdr")


# 6. Create results_batch2_positive_constitutive table

results_batch2_positive_constitutive <- data.frame(
  metabolite = rownames(data),
  log2FC = log2FC,
  pvalue = pvalues,
  padj = padj
)

results_batch2_positive_constitutive$negLog10P <- -log10(results_batch2_positive_constitutive$padj)


# 7. Volcano plot

# Define significance
results_batch2_positive_constitutive$significant <- "Not Significant"
results_batch2_positive_constitutive$significant[results_batch2_positive_constitutive$log2FC > 1 & results_batch2_positive_constitutive$padj < 0.05] <- "Up"
results_batch2_positive_constitutive$significant[results_batch2_positive_constitutive$log2FC < -1 & results_batch2_positive_constitutive$padj < 0.05] <- "Down"

Volcano_batch2_positive_constitutive <- ggplot(results_batch2_positive_constitutive, aes(x = log2FC, y = negLog10P, color = significant)) +
  geom_point(size = 3, alpha = 0.8) +
  scale_color_manual(values = c("blue", "grey", "red")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  theme_minimal() +
  labs(
    x = "Log2 Fold Change (E2-R / E2-S)",
    y = "-Log10(FDR)",
    title = "Volcano Plot"
  ) +
  geom_text_repel(
    data = subset(results_batch2_positive_constitutive, significant != "Not Significant"),
    aes(label = metabolite)
  )
results_batch2_positive_constitutive$significant<-as.factor(results_batch2_positive_constitutive$significant)
summary(results_batch2_positive_constitutive)




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
    data = subset(results_batch2_positive_susceptible, significant != "Not Significant"),
    aes(label = metabolite)
  )
results_batch2_positive_susceptible$significant<-as.factor(results_batch2_positive_susceptible$significant)
summary(results_batch2_positive_susceptible)



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
    data = subset(results_batch2_positive_resistant, significant != "Not Significant"),
    aes(label = metabolite)
  )
results_batch2_positive_resistant$significant<-as.factor(results_batch2_positive_resistant$significant)
summary(results_batch2_positive_resistant)






#bind ----

ID_positive <- read_tsv('metabolites_identification/SIRIUS_Positive_mode/canopus_formula_summary.tsv')
colnames(ID_positive)[colnames(ID_positive) == "mappingFeatureId"] <- "metabolite"

ID_negative <- read_tsv('../Negative mode/metabolites identification/SIRIUS_Negative_mode/canopus_formula_summary.tsv')
colnames(ID_negative)[colnames(ID_negative) == "mappingFeatureId"] <- "metabolite"


results_batch1_negative_resistant_sub <- results_batch1_negative_resistant[results_batch1_negative_resistant$significant != "Not Significant", ]
colnames(results_batch1_negative_resistant_sub)[colnames(results_batch1_negative_resistant_sub) == "significant"] <- "significant_resistant"
summary(results_batch1_negative_resistant_sub)
results_batch1_negative_susceptible_sub <- results_batch1_negative_susceptible[results_batch1_negative_susceptible$significant != "Not Significant", ]
colnames(results_batch1_negative_susceptible_sub)[colnames(results_batch1_negative_susceptible_sub) == "significant"] <- "significant_susceptible"

results_batch1_negative_susceptible_resistant_sub <- merge(results_batch1_negative_resistant_sub,results_batch1_negative_susceptible_sub, by ="metabolite",all=TRUE)
results_batch1_negative_susceptible_resistant_sub_ID <- merge(results_batch1_negative_susceptible_resistant_sub, ID_negative, by ="metabolite",all=FALSE)


results_batch2_negative_resistant_sub <- results_batch2_negative_resistant[results_batch2_negative_resistant$significant != "Not Significant", ]
colnames(results_batch2_negative_resistant_sub)[colnames(results_batch2_negative_resistant_sub) == "significant"] <- "significant_resistant"
summary(results_batch2_negative_resistant_sub)
results_batch2_negative_susceptible_sub <- results_batch2_negative_susceptible[results_batch2_negative_susceptible$significant != "Not Significant", ]
colnames(results_batch2_negative_susceptible_sub)[colnames(results_batch2_negative_susceptible_sub) == "significant"] <- "significant_susceptible"

results_batch2_negative_susceptible_resistant_sub <- merge(results_batch2_negative_resistant_sub,results_batch2_negative_susceptible_sub, by ="metabolite",all=TRUE)
results_batch2_negative_susceptible_resistant_sub_ID <- merge(results_batch2_negative_susceptible_resistant_sub, ID_negative, by ="metabolite",all.x = TRUE)



results_batch1_positive_resistant_sub <- results_batch1_positive_resistant[results_batch1_positive_resistant$significant != "Not Significant", ]
colnames(results_batch1_positive_resistant_sub)[colnames(results_batch1_positive_resistant_sub) == "significant"] <- "significant_resistant"
summary(results_batch1_positive_resistant_sub)
results_batch1_positive_susceptible_sub <- results_batch1_positive_susceptible[results_batch1_positive_susceptible$significant != "Not Significant", ]
colnames(results_batch1_positive_susceptible_sub)[colnames(results_batch1_positive_susceptible_sub) == "significant"] <- "significant_susceptible"

results_batch1_positive_susceptible_resistant_sub <- merge(results_batch1_positive_resistant_sub,results_batch1_positive_susceptible_sub, by ="metabolite",all=TRUE)
results_batch1_positive_susceptible_resistant_sub_ID <- merge(results_batch1_positive_susceptible_resistant_sub, ID_positive, by ="metabolite",all.x = TRUE)


results_batch2_positive_resistant_sub <- results_batch2_positive_resistant[results_batch2_positive_resistant$significant != "Not Significant", ]
colnames(results_batch2_positive_resistant_sub)[colnames(results_batch2_positive_resistant_sub) == "significant"] <- "significant_resistant"
summary(results_batch2_positive_resistant_sub)
results_batch2_positive_susceptible_sub <- results_batch2_positive_susceptible[results_batch2_positive_susceptible$significant != "Not Significant", ]
colnames(results_batch2_positive_susceptible_sub)[colnames(results_batch2_positive_susceptible_sub) == "significant"] <- "significant_susceptible"

results_batch2_positive_susceptible_resistant_sub <- merge(results_batch2_positive_resistant_sub,results_batch2_positive_susceptible_sub, by ="metabolite",all=TRUE)
results_batch2_positive_susceptible_resistant_sub_ID <- merge(results_batch2_positive_susceptible_resistant_sub, ID_positive, by ="metabolite",all.x = TRUE)

####################################################################################


# Create a list of metabolite names for each group
up_susceptible <- results_batch1_negative_susceptible_resistant_sub_ID %>%
  filter(significant_susceptible == "Up") %>%
  pull(metabolite)

down_susceptible <- results_batch1_negative_susceptible_resistant_sub_ID %>%
  filter(significant_susceptible == "Down") %>%
  pull(metabolite)

up_resistant <- results_batch1_negative_susceptible_resistant_sub_ID %>%
  filter(`significant_resistant` == "Up") %>%
  pull(metabolite)

down_resistant <- results_batch1_negative_susceptible_resistant_sub_ID %>%
  filter(significant_resistant == "Down") %>%
  pull(metabolite)

# Name the lists
list_sets <- list(
  "Up Resistant" = up_resistant,
  "Down Susceptible" = down_susceptible,
  "Down Resistant" = down_resistant,
  "Up Susceptible" = up_susceptible
)


batch1_negative_venn<-ggvenn(list_sets, fill_color = c("#1b9e77", "#d95f02", "#7570b3", "#e7298a"),
                             stroke_size = 0.5, set_name_size = 4.5,show_percentage = FALSE)





# Create a list of metabolite names for each group
up_susceptible <- results_batch1_positive_susceptible_resistant_sub_ID %>%
  filter(significant_susceptible == "Up") %>%
  pull(metabolite)

down_susceptible <- results_batch1_positive_susceptible_resistant_sub_ID %>%
  filter(significant_susceptible == "Down") %>%
  pull(metabolite)

up_resistant <- results_batch1_positive_susceptible_resistant_sub_ID %>%
  filter(significant_resistant == "Up") %>%
  pull(metabolite)

down_resistant <- results_batch1_positive_susceptible_resistant_sub_ID %>%
  filter(significant_resistant == "Down") %>%
  pull(metabolite)

# Name the lists
list_sets <- list(
  "Up Resistant" = up_resistant,
  "Down Susceptible" = down_susceptible,
  "Down Resistant" = down_resistant,
  "Up Susceptible" = up_susceptible
)


batch1_positive_venn<-ggvenn(list_sets, fill_color = c("#1b9e77", "#d95f02", "#7570b3", "#e7298a"),
                             stroke_size = 0.5, set_name_size = 4.5,show_percentage = FALSE)




# Create a list of metabolite names for each group
up_susceptible <- results_batch2_negative_susceptible_resistant_sub_ID %>%
  filter(significant_susceptible == "Up") %>%
  pull(metabolite)

down_susceptible <- results_batch2_negative_susceptible_resistant_sub_ID %>%
  filter(significant_susceptible == "Down") %>%
  pull(metabolite)

up_resistant <- results_batch2_negative_susceptible_resistant_sub_ID %>%
  filter(`significant_resistant` == "Up") %>%
  pull(metabolite)

down_resistant <- results_batch2_negative_susceptible_resistant_sub_ID %>%
  filter(significant_resistant == "Down") %>%
  pull(metabolite)

# Name the lists
list_sets <- list(
  "Up Resistant" = up_resistant,
  "Down Susceptible" = down_susceptible,
  "Down Resistant" = down_resistant,
  "Up Susceptible" = up_susceptible
)


batch2_negative_venn<-ggvenn(list_sets, fill_color = c("#1b9e77", "#d95f02", "#7570b3", "#e7298a"),
                             stroke_size = 0.5, set_name_size = 4.5,show_percentage = FALSE)





# Create a list of metabolite names for each group
up_susceptible <- results_batch2_positive_susceptible_resistant_sub_ID %>%
  filter(significant_susceptible == "Up") %>%
  pull(metabolite)

down_susceptible <- results_batch2_positive_susceptible_resistant_sub_ID %>%
  filter(significant_susceptible == "Down") %>%
  pull(metabolite)

up_resistant <- results_batch2_positive_susceptible_resistant_sub_ID %>%
  filter(`significant_resistant` == "Up") %>%
  pull(metabolite)

down_resistant <- results_batch2_positive_susceptible_resistant_sub_ID %>%
  filter(significant_resistant == "Down") %>%
  pull(metabolite)

# Name the lists
list_sets <- list(
  "Up Resistant" = up_resistant,
  "Down Susceptible" = down_susceptible,
  "Down Resistant" = down_resistant,
  "Up Susceptible" = up_susceptible
)


batch2_positive_venn<-ggvenn(list_sets, fill_color = c("#1b9e77", "#d95f02", "#7570b3", "#e7298a"),
                             stroke_size = 0.5, set_name_size = 4.5,show_percentage = FALSE)


batch1_negative_venn
batch1_positive_venn
batch2_negative_venn
batch2_positive_venn










combined <- bind_rows(results_batch2_positive_susceptible_resistant_sub_ID, results_batch2_negative_susceptible_resistant_sub_ID)



# =============================================
# 1. Define the grouping logic for Up and Down
# =============================================

# Function to create groups for a given regulation type ("Up" or "Down")
create_groups <- function(data, regulation) {
  data %>%
    mutate(
      # Trim whitespace from the regulation columns to avoid mismatches
      significant_susceptible = trimws(significant_susceptible),
      significant_resistant = trimws(significant_resistant),
      group = case_when(
        significant_susceptible == regulation & (significant_resistant != regulation | is.na(significant_resistant)) ~
          paste(regulation, "in significant_susceptible only"),
        significant_resistant == regulation & (significant_susceptible != regulation | is.na(significant_susceptible)) ~
          paste(regulation, "in significant_resistant only"),
        significant_susceptible == regulation & significant_resistant == regulation ~
          paste(regulation, "in both"),
        TRUE ~ NA_character_
      )
    ) %>%
    filter(!is.na(group)) %>%
    group_by(group, `NPC#pathway`) %>%
    summarise(count = n(), .groups = "drop")
}

# =============================================
# 2. Create groups for Up and Down
# =============================================

# Up groups
up_groups <- create_groups(combined, "Up")

# Down groups
down_groups <- create_groups(combined, "Down")

# Debug: Print the first few rows of down_groups
print("Down groups:")
print(head(down_groups))

# =============================================
# 3. Prepare data for plotting
# =============================================

# Get all unique NPC#pathway levels
all_pathways <- unique(combined$`NPC#pathway`)

# Create a data frame with all possible groups and pathways for Up
all_up_combinations <- expand.grid(
  group = c(
    "Up in significant_susceptible only",
    "Up in significant_resistant only",
    "Up in both"
  ),
  `NPC#pathway` = all_pathways
)

# Merge with actual counts (filling missing combinations with 0)
up_groups <- all_up_combinations %>%
  left_join(up_groups, by = c("group", "NPC#pathway")) %>%
  replace_na(list(count = 0))

# Create a data frame with all possible groups and pathways for Down
all_down_combinations <- expand.grid(
  group = c(
    "Down in significant_susceptible only",
    "Down in significant_resistant only",
    "Down in both"
  ),
  `NPC#pathway` = all_pathways
)

# Merge with actual counts (filling missing combinations with 0)
down_groups <- all_down_combinations %>%
  left_join(down_groups, by = c("group", "NPC#pathway")) %>%
  replace_na(list(count = 0))

# Debug: Print the first few rows of down_groups after merging
print("Down groups after merging:")
print(head(down_groups))

# =============================================
# 4. Plot stacked bars for Up and Down
# =============================================

# Plot for Up
up_plot <- ggplot(up_groups, aes(x = group, y = count, fill = `NPC#pathway`)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(
    title = "Up Regulation Groups by NPC Pathway",
    x = "",
    y = "Count"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Plot for Down
down_plot <- ggplot(down_groups, aes(x = group, y = count, fill = `NPC#pathway`)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(
    title = "Down Regulation Groups by NPC Pathway",
    x = "",
    y = "Count"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Display plots
print(up_plot)
print(down_plot)
