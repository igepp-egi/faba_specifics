library(lme4)
library(multcomp)
library(emmeans)

groups<-read.csv("00.info/metadata.csv", sep= ",",header = FALSE)
colnames(groups) <- c("ivip","groups")

emmean_fabae <- read.csv("02.analyses/01.GWAS/01_phenotipic_analysis/02.ajusted_means/fabae/emmeans_fabae_4.csv", sep =",")
emmean_pisum <- read.csv("02.analyses/01.GWAS/01_phenotipic_analysis/02.ajusted_means/pisum/emmeans_pisum_6.csv", sep =",")
  
merge_pisum_emmean<-merge(groups,emmean_pisum, by ="ivip", all.x = TRUE)
merge_fabae_emmean<-merge(groups,emmean_fabae, by ="ivip", all.x = TRUE)
merge_pisum_emmean$groups <- as.factor(merge_pisum_emmean$groups)
merge_fabae_emmean$groups <- as.factor(merge_fabae_emmean$groups)

library(ggplot2)
library(patchwork)

outdir <- "02.analyses/00.VCF_and_Population_Analysis/03.pop_structure_based_on_snps/03_phenotype_group_comparison"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

make_group_plot <- function(data, plot_title) {
  
  ggplot(data, aes(x = factor(groups), y = emmean)) +
    geom_boxplot(fill = "grey85", colour = "black", outlier.shape = NA) +
    geom_jitter(width = 0.15, shape = 1, colour = "black", size = 2) +
    scale_y_continuous(limits = c(0, 45), breaks = seq(0, 45, by = 5)) +
    labs(title = plot_title,
         x = "Group",
         y = "Adjusted mean of aphid count for each genotype") +
    theme_classic(base_size = 18) +
    theme(
      axis.text         = element_text(size = 18, colour = "black"),
      axis.title        = element_text(size = 20),
      axis.ticks        = element_line(linewidth = 0.8),
      axis.ticks.length = unit(0.25, "cm"),
      axis.line         = element_line(linewidth = 0.8),
      plot.title        = element_text(size = 22, face = "bold.italic")
    )
}

p_pisum <- make_group_plot(merge_pisum_emmean, "A. pisum")
p_fabae <- make_group_plot(merge_fabae_emmean, "A. fabae")

# pisum first (left), fabae second (right)
p_groups <- p_pisum | p_fabae

ggsave(file.path(outdir, "boxplot_groups_comparison_pisum_fabae.svg"),
       plot = p_groups, device = "svg", width = 14, height = 6)


pheno_fabae <- read.csv("02.analyses/01.GWAS/01_phenotipic_analysis/01.dataset_filtering/01_A.fabae/final filters/aphis_fabae_data_filtered_4.csv", sep =",")
pheno_pisum <- read.csv("02.analyses/01.GWAS/01_phenotipic_analysis/01.dataset_filtering/02_A.pisum/final filters/acyrthosiphon_pisum_filtered_6.csv", sep =",")

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
