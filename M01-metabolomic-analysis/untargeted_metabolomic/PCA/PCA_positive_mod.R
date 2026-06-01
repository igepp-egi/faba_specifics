#load packages ----
library(dplyr)
library(tidyr)
library(ggplot2)
library(factoextra)
library(patchwork)
library(factoextra)
library(viridis)
library(readxl)

setwd("~/Downloads/Paul Fremaux/untargeted/Positive mode/final_filtered_files")




# load datasets ----
batch1 <- read_excel("batch_1_positive_Intensity_filtered_notame.xlsx", col_names = FALSE)
batch2 <- read_excel("batch_2_positive_Intensity_filtered_notame.xlsx", col_names = FALSE)
batch3 <- read_excel("batch_3_positive_Intensity_filtered_notame.xlsx", col_names = FALSE)

#prepare the datasets ----
batch1 <- batch1 %>%
  mutate(across(34, ~ ifelse(!is.na(batch1[[3]]) & batch1[[3]] != "", batch1[[3]], .)))
batch1 <- batch1[,-c(1:33)]
batch1 <- t(batch1)
batch1 <- as.data.frame(batch1)
colnames(batch1) <- batch1[1, ]

batch1 <- batch1[-1, ]
batch1_row_name <- batch1 %>%
  filter(Group != "Blank") %>%
  separate(Group, into = c("Batch","pheno","treatmentA","hpi"), sep = "-") %>%
  mutate(genotype = paste(Batch, pheno, sep = "-")) %>%
  mutate(treatment = paste(Batch, treatmentA, sep = "-")) %>%
  select(FEATURE_ID,genotype,treatment, hpi)
batch1 <- merge(batch1_row_name,batch1, by ="FEATURE_ID")


batch2 <- batch2 %>%
  mutate(across(35, ~ ifelse(!is.na(batch2[[3]]) & batch2[[3]] != "", batch2[[3]], .)))
batch2 <- batch2[,-c(1:34)]
batch2 <- t(batch2)
batch2 <- as.data.frame(batch2)
colnames(batch2) <- batch2[1, ]
batch2 <- batch2[-1, ]
batch2_row_name <- batch2 %>%
  filter(Group != "Blank") %>%
  separate(Group, into = c("Batch","pheno","treatmentA","hpi"), sep = "-") %>%
  mutate(genotype = paste(Batch, pheno, sep = "-")) %>%
  mutate(treatment = paste(Batch, treatmentA, sep = "-")) %>%
  select(FEATURE_ID,genotype,treatment, hpi)
batch2 <- merge(batch2_row_name,batch2, by ="FEATURE_ID")
batch2 <- batch2[batch2$treatment != "QC_-NA", ]

batch3 <- batch3 %>%
  mutate(across(36, ~ ifelse(!is.na(batch3[[3]]) & batch3[[3]] != "", batch3[[3]], .)))
batch3 <- batch3[,-c(1:35)]
batch3 <- t(batch3)
batch3 <- as.data.frame(batch3)
colnames(batch3) <- batch3[1, ]
batch3 <- batch3[-1, ]
batch3_row_name <- batch3 %>%
  filter(Group != "Blank") %>%
  separate(Group, into = c("Batch","treatment","hpi"), sep = "-") %>%
  select(FEATURE_ID,Batch,treatment, hpi)
batch3 <- merge(batch3_row_name,batch3, by ="FEATURE_ID")



## PCA ----
#PCA batch1----
# Step 1: Select quantitative data (columns 9 onward)
quant_data <- batch1[, 9:ncol(batch1)]
quant_data <- apply(quant_data, 2, as.numeric)  # ensure numeric
quant_data[is.na(quant_data)] <- 0



# 3. Power transformation: square root
quant_data_log10 <- log10(quant_data + 1)

# 4. Centering: subtract metabolite mean (column-wise)
col_means <- colMeans(quant_data_log10, na.rm = TRUE)
quant_data_centered <- sweep(quant_data_log10, 2, col_means, "-")

# 5. Autoscaling: divide by metabolite SD (column-wise)
col_sds <- apply(quant_data_centered, 2, sd, na.rm = TRUE)
col_sds[col_sds == 0] <- 1  # prevent division by zero
quant_data <- sweep(quant_data_centered, 2, col_sds, "/")


rownames(quant_data) <- batch1$sample  # optional row names

# Step 2: Run PCA
pca_res <- prcomp(quant_data, scale. = FALSE)

pca_batch1 <- fviz_pca_ind(
  pca_res,
  geom.ind = "point",
  pointshape = 16,
  addEllipses = FALSE,
  mean.point = FALSE,
  axes = c(1, 2),
  legend.title = "hpi"
) +
  theme_minimal()

pca_hpi_batch1 <- fviz_pca_ind(
  pca_res,
  geom.ind = "point",
  habillage = batch1$hpi,
  palette = c("green", "red","yellow", "blue"),
  pointshape = 16,
  addEllipses = FALSE,
  mean.point = FALSE,
  axes = c(1, 2),
  legend.title = "hpi"
) +
  theme_minimal()


pca_genotype_batch1 <- fviz_pca_ind(
  pca_res,
  geom.ind = "point",
  habillage = batch1$genotype,
  palette = c("green", "red"),
  pointshape = 16,
  addEllipses = FALSE,
  mean.point = FALSE,
  axes = c(1, 2),
  legend.title = "genotype"
) +
  theme_minimal()

pca_treatment_batch1 <- fviz_pca_ind(
  pca_res,
  geom.ind = "point",
  habillage = batch1$treatment,
  palette = c("green", "red"),
  pointshape = 16,
  addEllipses = FALSE,
  mean.point = FALSE,
  axes = c(1, 2),
  legend.title = "treatment"
) +
  theme_minimal()



# Step 1: Run PCA
pca_res <- prcomp(quant_data, scale. = FALSE)

# Step 2: Extract variance explained
var_explained <- summary(pca_res)$importance[2, 1:2] * 100
pc1_var <- round(var_explained[1], 1)
pc2_var <- round(var_explained[2], 1)

# Step 3: Extract PCA coordinates
pca_coords <- as.data.frame(pca_res$x)
pca_coords$genotype <- batch1$genotype
pca_coords$treatment <- batch1$treatment

# Step 4: Plot with ggplot2
library(ggplot2)

ggplot(pca_coords, aes(x = PC1, y = PC2)) +
  geom_point(
    aes(shape = genotype, color = treatment),
    size = 3
  ) +
  scale_shape_manual(values = c(16, 17)) +
  scale_color_manual(values = c("green", "red")) +
  labs(
    x = paste0("PC1 (", pc1_var, "%)"),
    y = paste0("PC2 (", pc2_var, "%)"),
    shape = "Genotype",
    color = "Treatment"
  ) +
  theme_minimal()

PCA_axes_batch1<-fviz_eig(pca_res, addlabels = TRUE, ylim = c(0, 100))


#PCA batch2----
# Step 1: Select quantitative data (columns 9 onward)
quant_data <- batch2[, 9:ncol(batch2)]
quant_data <- apply(quant_data, 2, as.numeric)  # ensure numeric
quant_data[is.na(quant_data)] <- 0

# 3. Power transformation: square root
quant_data_log10 <- log10(quant_data + 1)

# 4. Centering: subtract metabolite mean (column-wise)
col_means <- colMeans(quant_data_log10, na.rm = TRUE)
quant_data_centered <- sweep(quant_data_log10, 2, col_means, "-")

# 5. Autoscaling: divide by metabolite SD (column-wise)
col_sds <- apply(quant_data_centered, 2, sd, na.rm = TRUE)
col_sds[col_sds == 0] <- 1  # prevent division by zero
quant_data <- sweep(quant_data_centered, 2, col_sds, "/")

# Result: quant_data_scaled is ready for PCA / PLS-DA
rownames(quant_data) <- batch2$sample  # optional row names

# Step 2: Run PCA
pca_res <- prcomp(quant_data, scale. = FALSE)

pca_batch2 <- fviz_pca_ind(
  pca_res,
  geom.ind = "point",
  pointshape = 16,
  addEllipses = FALSE,
  mean.point = FALSE,
  axes = c(1, 2),
  legend.title = "hpi"
) +
  theme_minimal()

pca_hpi_batch2 <- fviz_pca_ind(
  pca_res,
  geom.ind = "point",
  habillage = batch2$hpi,
  palette = c("green", "red","yellow", "blue"),
  pointshape = 16,
  addEllipses = FALSE,
  mean.point = FALSE,
  axes = c(1, 2),
  legend.title = "hpi"
) +
  theme_minimal()


pca_genotype_batch2 <- fviz_pca_ind(
  pca_res,
  geom.ind = "point",
  habillage = batch2$genotype,
  palette = c("green", "red"),
  pointshape = 16,
  addEllipses = FALSE,
  mean.point = FALSE,
  axes = c(1, 12),
  legend.title = "genotype"
) +
  theme_minimal()

pca_treatment_batch2 <- fviz_pca_ind(
  pca_res,
  geom.ind = "point",
  habillage = batch2$treatment,
  palette = c("green", "red"),
  pointshape = 16,
  addEllipses = FALSE,
  mean.point = FALSE,
  axes = c(1, 12),
  legend.title = "treatment"
) +
  theme_minimal()





# Step 1: Run PCA
pca_res <- prcomp(quant_data, scale. = FALSE)

# Step 2: Extract variance explained for PC1 and PC12
var_explained <- summary(pca_res)$importance[2, ]
pc1_var <- round(var_explained[1], 1)
pc12_var <- round(var_explained[12], 1)

# Step 3: Extract PCA coordinates for PC1 and PC12
pca_coords <- as.data.frame(pca_res$x)
pca_coords$genotype <- batch2$genotype
pca_coords$treatment <- batch2$treatment

# Step 4: Plot with ggplot2
library(ggplot2)

ggplot(pca_coords, aes(x = PC1, y = PC12)) +
  geom_point(
    aes(shape = genotype, color = treatment),
    size = 3
  ) +
  scale_shape_manual(values = c(16, 17)) +
  scale_color_manual(values = c("green", "red")) +
  labs(
    x = paste0("PC1 (", pc1_var, "%)"),
    y = paste0("PC12 (", pc12_var, "%)"),
    shape = "Genotype",
    color = "Treatment"
  ) +
  theme_minimal()

PCA_axes_batch1<-fviz_eig(pca_res, addlabels = TRUE, ylim = c(0, 100))



#PCA batch3----
# Step 1: Select quantitative data (columns 9 onward)
quant_data <- batch3[, 9:ncol(batch3)]
quant_data <- apply(quant_data, 2, as.numeric)  # ensure numeric
quant_data[is.na(quant_data)] <- 0


# 3. Power transformation: square root
quant_data_log10 <- log10(quant_data + 1)

# 4. Centering: subtract metabolite mean (column-wise)
col_means <- colMeans(quant_data_log10, na.rm = TRUE)
quant_data_centered <- sweep(quant_data_log10, 2, col_means, "-")

# 5. Autoscaling: divide by metabolite SD (column-wise)
col_sds <- apply(quant_data_centered, 2, sd, na.rm = TRUE)
col_sds[col_sds == 0] <- 1  # prevent division by zero
quant_data <- sweep(quant_data_centered, 2, col_sds, "/")


rownames(quant_data) <- batch3$sample  # optional row names

# Step 2: Run PCA
pca_res <- prcomp(quant_data, scale. = FALSE)


num_hpi <- length(unique(batch3$hpi))
palette <- viridis(num_hpi, option = "turbo")  # colorblind-friendly

pca_hpi_batch3 <- fviz_pca_ind(pca_res,
                                      geom.ind = "point",       # only points
                                      col.ind = batch3$hpi, # color by hpi
                                      palette = palette,
                                      pointshape = 16,          # all points same shape (solid circle)
                                      addEllipses = FALSE,      # remove ellipses if you want only points
                                      legend.title = "HPI") +
  theme_minimal()


num_treatment <- length(unique(batch3$treatment))
palette <- viridis(num_treatment, option = "turbo")  # colorblind-friendly

pca_treatment_batch3 <- fviz_pca_ind(pca_res,
                                            geom.ind = "point",       # only points
                                            col.ind = batch3$treatment, # color by treatment
                                            palette = palette,
                                            pointshape = 16,          # all points same shape (solid circle)
                                            addEllipses = FALSE,      # remove ellipses if you want only points
                                            legend.title = "treatment") +
  theme_minimal()

PCA_axes_batch3<-fviz_eig(pca_res, addlabels = TRUE, ylim = c(0, 100))

