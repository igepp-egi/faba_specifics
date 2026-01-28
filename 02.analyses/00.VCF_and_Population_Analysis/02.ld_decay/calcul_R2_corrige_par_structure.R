library(vcfR)
library(LDcorSV)
library(ggplot2)
library(tibble)

data <- read.table("../../../crosshap/metadata.csv",
                   header = FALSE, sep = ",")

subpops <- unique(data$V2)
P <- length(subpops)

N <- nrow(data)
S <- matrix(0, nrow = N, ncol = P - 1)
colnames(S) <- subpops[-P]

for (i in seq_len(N)) {
  for (j in seq_len(P - 1)) {
    if (data$V2[i] == subpops[j]) {
      S[i, j] <- 1
    }
  }
}

rownames(S) <- data$V1


vcf <- read.vcfR("region_subset.recode.vcf")

gt <- extract.gt(vcf, element = "GT")

Geno <- matrix(NA, nrow = nrow(gt), ncol = ncol(gt))

for (i in seq_len(nrow(gt))) {
  for (j in seq_len(ncol(gt))) {
    alleles <- unlist(strsplit(as.character(gt[i, j]), "[/|]"))
    Geno[i, j] <- sum(as.numeric(alleles))
  }
}

rownames(Geno) <- rownames(gt)
colnames(Geno) <- colnames(gt)

Geno <- t(Geno)


common <- intersect(rownames(S), rownames(Geno))

S    <- S[common, , drop = FALSE]
Geno <- Geno[common, , drop = FALSE]


best_snp <- "chr1S_1240520268"

best_vec <- Geno[, best_snp]

other_snps <- setdiff(colnames(Geno), best_snp)


r2_corrected <- sapply(other_snps, function(snp) {
  
  biloci <- data.frame(
    best  = best_vec,
    other = Geno[, snp]
  )
  
  keep <- complete.cases(biloci, S)
  
  Measure.R2S(
    biloci[keep, , drop = FALSE],
    S[keep, , drop = FALSE],
    na.presence = FALSE
  )
})


r2_corrected <- as.data.frame(r2_corrected)
r2_corrected <- rownames_to_column(r2_corrected, var = "SNP")

write.csv(r2_corrected, "r2_corrected_chr1S_1240520268.csv")

ggplot(r2_corrected,
       aes(x = SNP, y = r2_corrected, fill = r2_corrected)) +
  geom_col() +
  scale_fill_gradient(low = "red", high = "yellow") +
  theme_minimal() +
  labs(
    title = "LD structure-corrected (R²S)",
    x = "SNP",
    y = "R² corrigé",
    fill = "R²S"
  )

selected_snp <- r2_corrected[r2_corrected$r2_corrected > 0.75, ]

lines <- c(
  selected_snp$SNP,
  "chr1S_1240520268"
)

writeLines(lines, "test.txt")



