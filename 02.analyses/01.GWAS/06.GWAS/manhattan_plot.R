library(qqman)
df<-read.csv("../../../Downloads/GAPIT.Association.GWAS_Results.BLINK.emmean.csv"
             , sep = ",",header = TRUE, stringsAsFactors = FALSE)

df$P_Bonferroni <- p.adjust(df$P.value, method = "bonferroni")
df$P_BH <- p.adjust(df$P.value, method = "BH")        # Benjamini-Hochberg

threshold_bonf <- max(df$P.value[df$P_Bonferroni < 0.05], na.rm = TRUE)
threshold_bh   <- max(df$P.value[df$P_BH < 0.05], na.rm = TRUE)

# Convert to -log10(p)
thr_bonf_log <- -log10(threshold_bonf)
thr_bh_log   <- log10(threshold_bh)
test   <- -log10(0.05/(nrow(df)-1))



# Plot base
manhattan(df,
          chr = "Chr", bp = "Pos", P = "P.value", snp = "SNP",
          genomewideline = FALSE,
          suggestiveline = FALSE,
          main = "Manhattan Plot with Multiple Corrections")
df$chr_numeric <- as.numeric(factor(df$Chr))

manhattan(df,
          chr = "chr_numeric", bp = "Pos", p = "P.value", snp = "SNP",
          genomewideline = FALSE,
          suggestiveline = FALSE,
          main = "Manhattan Plot with Multiple Corrections")
summary(df)

# Add custom thresholds
abline(h = thr_bonf_log, col = "red", lty = 2)  # Bonferroni
abline(h = thr_bh_log, col = "blue", lty = 2)   # BH
abline(h = thr_fdr_log, col = "green", lty = 2) # FDR
legend("topright", legend = c("Bonferroni", "BH", "FDR"),
       col = c("red", "blue", "green"), lty = 2, cex = 0.8)
