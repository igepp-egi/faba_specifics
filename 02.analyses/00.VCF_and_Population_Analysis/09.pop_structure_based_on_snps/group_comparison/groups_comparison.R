groups<-read.csv("../../../../00.info/11.infos_genotypes/metadata.csv", sep= ",",header = FALSE)
colnames(groups) <- c("ivip","groups")

pheno_fabae <- read.csv("../../../../02.analyses/01.GWAS/01_phenotipic_analysis/01.dataset_filtering/01_A.fabae/final filters/aphis_fabae_data_filtered_4.csv", sep =",")
pheno_pisum <- read.csv("../../../../02.analyses/01.GWAS/01_phenotipic_analysis/01.dataset_filtering/02_A.pisum/final filters/acyrthosiphon_pisum_filtered_6.csv", sep =",")

merge_pisum<-merge(groups,pheno_pisum, by ="ivip", all.x = TRUE)
merge_fabae<-merge(groups,pheno_fabae, by ="ivip", all.x = TRUE)
merge_pisum$groups <- as.factor(merge_pisum$groups)
merge_fabae$groups <- as.factor(merge_fabae$groups)

# Pisum
postscript(
  "boxplot_groups_comparison_pisum.eps",
  onefile = FALSE,
  paper = "special",
  width = 7,
  height = 6
)

boxplot(
  larvae_total_number ~ groups,
  data = merge_pisum,
  main = "Boxplot of aphid count by Group",
  xlab = "Group",
  ylab = "aphid count",
  col = "grey"
)

stripchart(
  larvae_total_number ~ groups,
  data = merge_pisum,
  vertical = TRUE,
  method = "jitter",
  pch = 14,
  col = "black",
  add = TRUE
)

dev.off()


# Fabae
postscript(
  "boxplot_groups_comparison_fabae.eps",
  onefile = FALSE,
  paper = "special",
  width = 7,
  height = 6
)

boxplot(
  larvae_total_number ~ groups,
  data = merge_fabae,
  main = "Boxplot of aphid count by Group",
  xlab = "Group",
  ylab = "aphid count",
  col = "grey"
)

stripchart(
  larvae_total_number ~ groups,
  data = merge_fabae,
  vertical = TRUE,
  method = "jitter",
  pch = 14,
  col = "black",
  add = TRUE
)

dev.off()

# For fabae
mod_fabae <- lm(larvae_total_number ~ groups + session/block, data = merge_fabae)
anova_fabae <- anova(mod_fabae)
print(anova_fabae)

# Run Tukey HSD
tukey_fabae <- TukeyHSD(aov(larvae_total_number ~ groups, data = merge_fabae))
print(tukey_fabae)

# For pisum
mod_pisum <- lm(larvae_total_number ~ groups + session/block, data = merge_pisum)
anova_pisum <- anova(mod_pisum)
print(anova_pisum)

# Run Tukey HSD
tukey_pisum <- TukeyHSD(aov(larvae_total_number ~ groups, data = merge_pisum))
print(tukey_pisum)
