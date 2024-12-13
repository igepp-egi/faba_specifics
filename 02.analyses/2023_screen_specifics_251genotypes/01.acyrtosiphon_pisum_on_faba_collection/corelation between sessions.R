
#data preparation ----
data250pisum<-read.csv("acyrthospipon_pisum_fecundity_on_faba_collection.csv", sep =",")
data250pisum$genotype<-as.factor(data250pisum$genotype)
data250pisum <- data250pisum %>% filter(!is.na(larvae_total_number))
data_filtered_250_pisum <- data250pisum %>% filter(adult_alive != "0") #remove if no adult


data_filtered_250_pisum <- data_filtered_250_pisum %>%
  group_by(short_code) %>%           # Group by genotype
  filter(n() >= 6) %>%             # Keep only groups with 6 or more repetitions
  ungroup()  
#Package----
library(corrplot)
library(Hmisc)
# Step 1: Summarize duplicates by calculating the mean----
datacorpisum250<-data_filtered_250_pisum[,c(2,5,12)]
datacorpisum250_summarized <- datacorpisum250 %>%
  group_by(genotype, session) %>% # Group by the relevant columns
  summarise(
    larvae_total_number = mean(larvae_total_number, na.rm = TRUE), # Calculate the mean
    .groups = "drop" # Remove grouping after summarizing
  )


# Step 2: Pivot the summarized dataset----
wide_datapisum250 <- datacorpisum250_summarized %>%
  pivot_wider(
    names_from = session,        # Use 'replicate' to create new columns
    values_from = larvae_total_number # Fill new columns with the mean values
  )

#step 3:calcul of the correlation between session ----

corpisum250<-wide_datapisum250[,-c(1)]

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

rest <- rcorr(as.matrix(corpisum250[,1:7]),type="pearson") # « pearson » : données quantitatives, « spearman » : données ordinales
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
