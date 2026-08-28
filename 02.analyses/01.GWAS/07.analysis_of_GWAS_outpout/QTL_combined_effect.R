library(pheatmap)
library(ggplot2)
library(openxlsx)
library(lme4)
library(emmeans)
library(multcomp)
library(multcompView)

cluster<-read.xlsx("../../../00.info/11.infos_genotypes/metadata_passportdata_emmeans.xlsx")
chr1<-read.csv("../../../02.analyses/01.GWAS/07.analysis_of_GWAS_outpout/analysis of GWAS outpout fabae /chr1L_1703281021/chr1L_final_window_genotype_groups.csv")
chr3<-read.csv("../../../02.analyses/01.GWAS/07.analysis_of_GWAS_outpout/analysis of GWAS outpout fabae /chr3_641834492/haplo_chr3_641834492.csv")
chr5<-read.csv("../../../02.analyses/01.GWAS/07.analysis_of_GWAS_outpout/analysis of GWAS outpout fabae /chr5_1299427225/haplo_chr5_1299427225.csv")
chr6<-read.csv("../../../02.analyses/01.GWAS/07.analysis_of_GWAS_outpout/analysis of GWAS outpout fabae /chr6_369308865/haplo_chr6_369308865.csv")
emmean_fabae<-read.csv("../../../02.analyses/01.GWAS/01_phenotipic_analysis/02.ajusted_means/fabae/emmeans_fabae_4.csv")
emmean_fabae<-emmean_fabae[,c(1,2)]
emmean_fabae$ivip[emmean_fabae$ivip == "hedin"] <- "Hedin"
merge1<-merge(chr1,chr3, by="ivip",all = TRUE)
merge2<-merge(merge1,chr5, by="ivip",all = TRUE)
merge3<-merge(merge2,chr6, by="ivip",all = TRUE)
merge_haplotype_fabae<-merge(merge3,emmean_fabae, by="ivip",all = TRUE)
merge_haplotype_fabae[is.na(merge_haplotype_fabae)] <- "Unfavorable"
write.xlsx(merge_haplotype_fabae,"../../../02.analyses/01.GWAS/07.analysis_of_GWAS_outpout/analysis of GWAS outpout fabae")
haplotype_cluster_fabae<-merge(merge_haplotype_fabae,cluster, by="ivip",all = TRUE)


# Count the number of Favorable alleles in columns 2-5
haplotype_cluster_fabae$Favorable_Count <- rowSums(haplotype_cluster_fabae[, 2:5] == "Favorable", na.rm = TRUE)

# Rename column 6 (optional)
names(haplotype_cluster_fabae)[6] <- "emmeans"

p <- ggplot(haplotype_cluster_fabae, aes(x = factor(Favorable_Count), y = emmeans)) +
  geom_jitter(
    width = 0.15,
    size = 2,
    colour = "black"
  ) +
  labs(
    x = "Number of favorable alleles",
    y = "Estimated marginal mean"
  ) +
  theme_classic(base_size = 11, base_family = "Helvetica")

ggsave(
  "../../../02.analyses/01.GWAS/07.analysis_of_GWAS_outpout/Figure_cumulated_effect_fabae.eps",
  plot = p,
  device = "eps",
  width = 3.5,
  height = 3.0,
  units = "in"
)



raw_data_fabae<-read.csv("../../../02.analyses/01.GWAS/01_phenotipic_analysis/01.dataset_filtering/01_A.fabae/final filters/aphis_fabae_data_filtered_4.csv")

merge5<-merge(merge3,raw_data_fabae, by="ivip",all = TRUE)
merge5$Favorable_Count <- rowSums(merge5[, 2:5] == "Favorable", na.rm = TRUE)
merge5$Favorable_Count<-as.factor(merge5$Favorable_Count)
mod <- lmer(larvae_total_number ~ Favorable_Count + (1 |ivip) + (1 |session/block), data = merge5)

emm <- emmeans(mod, "Favorable_Count")

letters_fabae <- cld(emm, adjust = "tukey", Letters = letters)
letters_fabae

chr1<-read.csv("../../../02.analyses/01.GWAS/07.analysis_of_GWAS_outpout/analysis of GWAS outpout pisum/chr1S_1240520268/haplo_chr1S_1240520268.csv")
chr2A<-read.csv("../../../02.analyses/01.GWAS/07.analysis_of_GWAS_outpout/analysis of GWAS outpout pisum/chr2_1461736473/haplo_chr2_1461736473.csv")
chr2B<-read.csv("../../../02.analyses/01.GWAS/07.analysis_of_GWAS_outpout/analysis of GWAS outpout pisum/chr2_35683372/haplo_chr2_35683372.csv")
chr2C<-read.csv("../../../02.analyses/01.GWAS/07.analysis_of_GWAS_outpout/analysis of GWAS outpout pisum/chr2_920208494/haplo_chr2_920208494.csv")

emmean_pisum<-read.csv("../../../02.analyses/01.GWAS/01_phenotipic_analysis/02.ajusted_means/pisum/emmeans_pisum_6.csv")
emmean_pisum<-emmean_pisum[,c(1,2)]
emmean_pisum$ivip[emmean_pisum$ivip == "hedin"] <- "Hedin"
merge1<-merge(chr1,chr2A, by="ivip",all = TRUE)
merge2<-merge(merge1,chr2B, by="ivip",all = TRUE)
merge3<-merge(merge2,chr2C, by="ivip",all = TRUE)
merge_haplotype_pisum<-merge(merge3,emmean_pisum, by="ivip",all = TRUE)
merge_haplotype_pisum[is.na(merge_haplotype_pisum)] <- "Unfavorable"
write.xlsx(merge_haplotype_pisum,"../../../02.analyses/01.GWAS/07.analysis_of_GWAS_outpout/analysis of GWAS outpout pisum")
haplotype_cluster_pisum<-merge(merge_haplotype_pisum,cluster, by="ivip",all = TRUE)


# Count the number of Favorable alleles in columns 2-5
haplotype_cluster_pisum$Favorable_Count <- rowSums(haplotype_cluster_pisum[, 2:5] == "Favorable", na.rm = TRUE)

# Rename column 6 
names(haplotype_cluster_pisum)[6] <- "emmeans"



p <- ggplot(haplotype_cluster_pisum, aes(x = factor(Favorable_Count), y = emmeans)) +
  geom_jitter(
    width = 0.15,
    size = 2,
    colour = "black"
  ) +
  labs(
    x = "Number of favorable alleles ",
    y = "Ajusted mean of aphid count"
  ) +
  theme_classic(base_size = 11, base_family = "Helvetica")

ggsave(
  "../../../02.analyses/01.GWAS/07.analysis_of_GWAS_outpout/Figure_cumulated_effect_pisum.eps",
  plot = p,
  device = "eps",
  width = 3.5,
  height = 3.0,
  units = "in"
)


raw_data_pisum<-read.csv("../../../02.analyses/01.GWAS/01_phenotipic_analysis/01.dataset_filtering/02_A.pisum/final filters/acyrthosiphon_pisum_filtered_6.csv")

merge5<-merge(merge3,raw_data_pisum, by="ivip",all = TRUE)
merge5$Favorable_Count <- rowSums(merge5[, 2:5] == "Favorable", na.rm = TRUE)
merge5$Favorable_Count<-as.factor(merge5$Favorable_Count)
mod <- lmer(larvae_total_number ~ Favorable_Count + (1 |ivip) + (1 |session/block), data = merge5)


emm <- emmeans(mod, "Favorable_Count")

letters_pisum <- cld(emm, adjust = "tukey", Letters = letters)
letters_pisum

all_haplotype<- merge(merge_haplotype_pisum, merge_haplotype_fabae, by = "ivip", all = TRUE)
all_haplotype_group <- merge(all_haplotype, cluster, by = "ivip", all = TRUE)
all_haplotype_group <- all_haplotype_group[!is.na(all_haplotype_group$group), ]
write.xlsx(all_haplotype_group,"../../../02.analyses/01.GWAS/07.analysis_of_GWAS_outpout/all_haplotype_group_emmean.xls")


