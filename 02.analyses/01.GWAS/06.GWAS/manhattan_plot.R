library(qqman)
df<-read.csv("../../../Downloads/GAPIT.Association.GWAS_Results.BLINK.emmean.csv"
             , sep = ",",header = TRUE, stringsAsFactors = FALSE)

df<-read.csv("../../../Documents/results_with_hedin/GWAS outpout/fabae/GAPIT.Association.GWAS_Results.BLINK.emmean(NYC).csv"
             , sep = ",",header = TRUE, stringsAsFactors = FALSE)

df$P_Bonferroni <- p.adjust(df$P.value, method = "bonferroni")
df$P_BH <- p.adjust(df$P.value, method = "BH")        # Benjamini-Hochberg

threshold_bonf <- max(df$P.value[df$P_Bonferroni < 0.05], na.rm = TRUE)
threshold_bh   <- max(df$P.value[df$P_BH < 0.05], na.rm = TRUE)

# Convert to -log10(p)
thr_bonf_log <- -log10(threshold_bonf)
thr_bh_log   <- -log10(threshold_bh)


df$P <- df$P.value
summary(df)




df$Chr_num <- NA

df$Chr_num[df$Chr == "chr1S"] <- 1
df$Chr_num[df$Chr == "chr1L"] <- 2
df$Chr_num[df$Chr == "chr2"]  <- 3
df$Chr_num[df$Chr == "chr3"]  <- 4
df$Chr_num[df$Chr == "chr4"]  <- 5
df$Chr_num[df$Chr == "chr5"]  <- 6
df$Chr_num[df$Chr == "chr6"]  <- 7

df$Chr_num <- as.numeric(df$Chr_num)

qtls <- c("chr1L_1703281021", "chr3_641834492", "chr5_1299427225", "chr6_369308865")


# Plot base
manhattan(df,
          chr = "Chr_num", bp = "Pos", P = "P", snp = "SNP", highlight = qtls,
          genomewideline = FALSE,
          suggestiveline = FALSE,
          main = "Manhattan Plot with Multiple Corrections")

# Add custom thresholds
abline(h = thr_bonf_log, col = "red", lty = 2)  # Bonferroni
abline(h = thr_bh_log, col = "blue", lty = 2)   # BH
abline(h = thr_fdr_log, col = "green", lty = 2) # FDR
legend("topright", legend = c("Bonferroni", "BH", "FDR"),
       col = c("red", "blue", "green"), lty = 2, cex = 0.8)

manhattan(df,
          chr = "Chr_num",
          bp  = "Pos",
          P   = "P",
          snp = "SNP",
          highlight = qtls,
          col = c("grey70", "grey40"),
          genomewideline = FALSE,
          suggestiveline = FALSE,
          main = "Manhattan Plot with Highlighted QTLs")
dev.off()  # reset graphics device

