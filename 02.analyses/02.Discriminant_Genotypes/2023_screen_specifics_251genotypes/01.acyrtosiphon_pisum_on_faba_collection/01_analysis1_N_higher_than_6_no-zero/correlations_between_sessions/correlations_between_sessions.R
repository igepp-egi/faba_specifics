# Package----
suppressPackageStartupMessages(library(corrplot))
suppressPackageStartupMessages(library(Hmisc))
suppressPackageStartupMessages(library(tidyverse))


# Data preparation ----
data250pisum_00 <- read.csv("01.data/2023_screen_specifics_251genotypes/01.acyrthosiphon_pisum_P123amp_on_faba_collection/acyrthospipon_pisum_fecundity_on_faba_collection.csv") %>% 
  mutate(genotype = as.factor(genotype)) %>%
  filter(!is.na(larvae_total_number)) %>%
  filter(adult_alive != "0") #remove if no adult

# Step 0: Filter the data if n replicates < 6 ----
data250pisum_01 <- data250pisum_00 %>%
  group_by(short_code) %>%           # Group by genotype
  filter(n() >= 6) %>%             # Keep only groups with 6 or more repetitions
  ungroup()  

# Step 1: Summarize duplicates by calculating the mean----
data250pisum_02 <- data250pisum_01 %>% 
  group_by(short_code, session) %>% # Group by the relevant columns
  summarise(
    larvae_total_number = mean(larvae_total_number, na.rm = TRUE), # Calculate the mean
    .groups = "drop" # Remove grouping after summarizing
  )

# Step 2: Pivot the summarized dataset----
data250pisum_03 <- data250pisum_02 %>%
  pivot_wider(
    names_from = session,              # Use 'session' to create new columns
    values_from = larvae_total_number  # Fill new columns with the mean values
  ) %>%
  dplyr::select(- short_code)  # Remove the 'short_code' column

# Step 3: computation of the correlation between sessions ----

# ++++++++++++++++++++++++++++
# flattenCorrMatrix
# ++++++++++++++++++++++++++++
# cormat : matrix of the correlation coefficients
# pmat : matrix of the correlation p-values
flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
  )
}

rest <- rcorr(as.matrix(data250pisum_03),type="pearson") # « pearson » : données quantitatives, « spearman » : données ordinales
flattenCorrMatrix(rest$r, rest$P)
M2 <- rest$r
diag(M2) = NA
corrplot(M2,
         title='Rep',
         na.label = "NA",
         type="upper",
         method = c("square"),
         p.mat = rest$P,
         sig.level = 0.05,
         tl.cex=1.2,
         tl.srt=90,
         tl.col="black",
         addCoef.col = "black",
         mar=c(10,0,1,0))
plot.new()
dev.off()
corrplot(M2, title = "Rep", na.label = "NA", type = "upper", tl.cex = 0.8)
