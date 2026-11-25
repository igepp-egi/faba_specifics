library(vcfR)
library(dplyr)
library(readxl)
library(writexl)

# ---------------- INPUTS (edit file names / paths) ----------------
snps <- read_excel("haplotype_snps.xlsx")        # columns: QTL, CHROM, POS (ID optional)
geno <- read_excel("genotypes.xlsx")[[1]]         # first column = genotypes of interest
refs <- read_excel("reference_genotypes.xlsx")   # columns: QTL, Reference

# folder containing the sub-VCFs (chr1L/, chr3/, chr5/, chr6/ created with bcftools)
vcf_dir <- "../../../../../../../nextcloud/2024-2027_Paul_PhD/01.data/01.SNPs_VCF_faba_collection"

# One truncated VCF per QTL. Names MUST match the QTL names in your SNP and reference tables exactly.
vcf_files <- c(
  "Af-Vf1L.1" = file.path(vcf_dir, "chr1L-1703271021-1703291123.vcf"),
  "Af-Vf3.1"  = file.path(vcf_dir, "chr3-641824491-641844492.vcf"),
  "Af-Vf5.1"  = file.path(vcf_dir, "chr5-1299417225-1299437278.vcf"),
  "Af-Vf6.1"  = file.path(vcf_dir, "chr6-369298865-369318865.vcf")
)

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
  
  idx    <- match(keys[found], rownames(gt))
  ref_gt <- gt[idx, ref]
  
  vapply(geno, function(g) {
    if (!g %in% colnames(gt)) return(NA_real_)
    g_gt  <- gt[idx, g]
    valid <- !is.na(ref_gt) & !is.na(g_gt)       # ignore missing calls
    if (!any(valid)) return(NA_real_)
    mean(g_gt[valid] == ref_gt[valid]) * 100     # % identity with the reference
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
           "haplotype_results2.xlsx")
print(head(res))