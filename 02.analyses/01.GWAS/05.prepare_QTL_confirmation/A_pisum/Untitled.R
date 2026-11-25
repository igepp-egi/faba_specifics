library(vcfR)
library(dplyr)
library(readxl)
library(writexl)

# ---------------- INPUTS (edit file names / paths) ----------------
snps <- read_excel("haplotype_snps1.xlsx")        # columns: QTL, CHROM, POS (ID optional)
geno <- read_excel("genotypes.xlsx")[[1]]        # first column = genotypes of interest
refs <- read_excel("reference_genotypes.xlsx")   # columns: QTL, Reference

# One truncated VCF per QTL. Names MUST match the QTL names in your SNP table exactly.
vcf_files <- c(
  "Ap-Vf1S.1" = "../../../../../../nextcloud/2024-2027_Paul_PhD/01.data/01.SNPs_VCF_faba_collection/chr1S-1240510268-1240530268.vcf",
  "Ap-Vf2.1" = "../../../../../../nextcloud/2024-2027_Paul_PhD/01.data/01.SNPs_VCF_faba_collection/chr2-35673372-35693372.vcf",
  "Ap-Vf2.2" = "../../../../../../nextcloud/2024-2027_Paul_PhD/01.data/01.SNPs_VCF_faba_collection/chr2-920108494-920308494.vcf",
  "Ap-Vf2.3" = "../../../../../../nextcloud/2024-2027_Paul_PhD/01.data/01.SNPs_VCF_faba_collection/chr2-1460736473-1462736473.vcf"
)

# QTLs where a heterozygous call (0/1, i.e. "1" in 0/1/2 coding) counts as
# compatible with BOTH homozygous states (0 and 2)
het_wildcard <- c("Ap-Vf2.3")

# ---------------- CLEAN NAMES / KEYS ----------------
# remove ANY number of leading "chr" (chr, chrchr, Chr...) and spaces
clean_chr <- function(x) sub("^(chr)+", "", trimws(as.character(x)), ignore.case = TRUE)

geno <- trimws(as.character(geno))
refs$QTL       <- trimws(as.character(refs$QTL))
refs$Reference <- trimws(as.character(refs$Reference))
snps$QTL       <- trimws(as.character(snps$QTL))


snps$key <- paste0(clean_chr(snps$CHROM), "_",
                   format(snps$POS, scientific = FALSE, trim = TRUE))

print(setdiff(unique(snps$QTL), names(vcf_files)))   # should be character(0)
print(setdiff(unique(snps$QTL), refs$QTL))           # should be character(0)
stopifnot(all(unique(snps$QTL) %in% names(vcf_files)))
stopifnot(all(file.exists(vcf_files)))

# heterozygous = the two alleles of the call differ (0/1, 1/0, 1/2 ...); NA -> FALSE
is_het <- function(x) {
  a <- sub("/.*", "", x)
  b <- sub(".*/", "", x)
  !is.na(x) & a != b
}

# ---------------- COMPARE ONE QTL (reads only its own small VCF) ----------------
compare_qtl <- function(q) {
  cat("\n=== ", q, " ===\n", sep = "")
  vcf <- read.vcfR(vcf_files[[q]], verbose = FALSE)
  
  vcf_key <- paste0(clean_chr(vcf@fix[, "CHROM"]), "_", trimws(vcf@fix[, "POS"]))
  gt <- extract.gt(vcf, element = "GT")
  gt <- matrix(gsub("\\|", "/", gt), nrow = nrow(gt), ncol = ncol(gt),
               dimnames = list(vcf_key, colnames(gt)))
  gt[gt %in% c("./.", ".")] <- NA
  
  keys  <- unique(snps$key[snps$QTL == q])
  found <- keys %in% rownames(gt)
  cat("SNPs found in VCF:", sum(found), "of", length(keys), "\n")
  cat("Genotypes found in VCF:", sum(geno %in% colnames(gt)), "of", length(geno), "\n")
  cat("Genotypes NOT in VCF:", paste(setdiff(geno, colnames(gt)), collapse = ", "), "\n")
  
  ref <- refs$Reference[refs$QTL == q]
  if (sum(found) == 0 || length(ref) != 1 || !ref %in% colnames(gt)) {
    warning("Cannot compute ", q, " (no SNPs found or reference missing in VCF)")
    return(rep(NA_real_, length(geno)))
  }
  
  use_wild <- q %in% het_wildcard
  if (use_wild) cat("Heterozygous calls treated as matching 0 and 2 for this QTL\n")
  
  idx    <- match(keys[found], rownames(gt))
  ref_gt <- gt[idx, ref]
  
  vapply(geno, function(g) {
    if (!g %in% colnames(gt)) return(NA_real_)
    g_gt  <- gt[idx, g]
    valid <- !is.na(ref_gt) & !is.na(g_gt)       # ignore missing calls
    if (!any(valid)) return(NA_real_)
    same <- g_gt == ref_gt
    if (use_wild) same <- same | is_het(g_gt) | is_het(ref_gt)
    mean(same[valid]) * 100                      # % identity with the reference
  }, numeric(1))
}

# ---------------- RUN ALL QTLs ----------------
qtls <- unique(snps$QTL)
res <- data.frame(Genotype = geno, stringsAsFactors = FALSE)
for (q in qtls) res[[q]] <- compare_qtl(q)

# Yes/No version: Yes = 100% identity with the reference haplotype
res_yesno <- res
res_yesno[qtls] <- lapply(res[qtls], function(x)
  ifelse(is.na(x), "NA", ifelse(x == 100, "Yes", "No")))

# ---------------- EXPORT ----------------
write_xlsx(list(Percent_identity = res, Carries_haplotype = res_yesno),
           "haplotype_results.xlsx")
print(head(res))







