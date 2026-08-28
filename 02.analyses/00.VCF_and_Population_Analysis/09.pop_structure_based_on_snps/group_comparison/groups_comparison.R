library(lme4)
library(multcomp)
library(emmeans)

groups<-read.csv("../../../../00.info/11.infos_genotypes/metadata.csv", sep= ",",header = FALSE)
colnames(groups) <- c("ivip","groups")

emmean_fabae <- read.csv("../../../../02.analyses/01.GWAS/01_phenotipic_analysis/02.ajusted_means/fabae/emmeans_fabae_4.csv", sep =",")
emmean_pisum <- read.csv("../../../../02.analyses/01.GWAS/01_phenotipic_analysis/02.ajusted_means/pisum/emmeans_pisum_6.csv", sep =",")
  
merge_pisum_emmean<-merge(groups,emmean_pisum, by ="ivip", all.x = TRUE)
merge_fabae_emmean<-merge(groups,emmean_fabae, by ="ivip", all.x = TRUE)
merge_pisum_emmean$groups <- as.factor(merge_pisum_emmean$groups)
merge_fabae_emmean$groups <- as.factor(merge_fabae_emmean$groups)

# Pisum
postscript(
  "boxplot_groups_comparison_pisum.eps",
  onefile = FALSE,
  paper = "special",
  width = 7,
  height = 6
)

boxplot(
  emmean ~ groups,
  data = merge_pisum_emmean,
  outline = FALSE,
  main = "Boxplot of aphid count by Group",
  xlab = "Group",
  ylab = "Ajusted mean of aphid count for each genotypes",
  col = "grey"
)

stripchart(
  emmean ~ groups,
  data = merge_pisum_emmean,
  vertical = TRUE,
  method = "jitter",
  pch = 1,
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
  emmean ~ groups,
  data = merge_fabae_emmean,
  outline = FALSE,
  main = "Boxplot of aphid count by Group",
  xlab = "Group",
  ylab = "Ajusted mean of aphid count for each genotypes",
  col = "grey"
)

stripchart(
  emmean ~ groups,
  data = merge_fabae_emmean,
  vertical = TRUE,
  method = "jitter",
  pch = 1,
  col = "black",
  add = TRUE
)

dev.off()




pheno_fabae <- read.csv("../../../../02.analyses/01.GWAS/01_phenotipic_analysis/01.dataset_filtering/01_A.fabae/final filters/aphis_fabae_data_filtered_4.csv", sep =",")
pheno_pisum <- read.csv("../../../../02.analyses/01.GWAS/01_phenotipic_analysis/01.dataset_filtering/02_A.pisum/final filters/acyrthosiphon_pisum_filtered_6.csv", sep =",")

merge_pisum<-merge(groups,pheno_pisum, by ="ivip", all.x = TRUE)
merge_fabae<-merge(groups,pheno_fabae, by ="ivip", all.x = TRUE)
merge_pisum$groups <- as.factor(merge_pisum$groups)
merge_fabae$groups <- as.factor(merge_fabae$groups)



# For fabae
mod_fabae <- lmer(larvae_total_number ~ groups + (1 |session/block) + (1 | groups:ivip), data = merge_fabae)

anova_fabae <- anova(mod_fabae)
print(anova_fabae)

# Estimated marginal means for groups
emm_fabae <- emmeans(mod_fabae, ~ groups)

# Tukey pairwise comparisons
pairs(emm_fabae, adjust = "tukey")

letters_fabae <- cld(
  emm_fabae,
  adjust = "tukey",
  Letters = letters
)

print(letters_fabae)




# For pisum
mod_pisum <- lmer(larvae_total_number ~ groups + (1 |session/block) + (1 | groups:ivip), data = merge_pisum)

anova_pisum <- anova(mod_pisum)
print(anova_pisum)

# Estimated marginal means for groups
emm_pisum <- emmeans(mod_pisum, ~ groups)

# Tukey pairwise comparisons
pairs(emm_pisum, adjust = "tukey")

letters_pisum <- cld(
  emm_pisum,
  adjust = "tukey",
  Letters = letters
)

print(letters_pisum)
