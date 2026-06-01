setwd("~/Documents/gitrepo/faba_specifics/00.info/11.infos_genotypes")
passport<-read.csv("../00.passeport_data/2024-11_passeport_data_and_genetic_cluster_viciafaba_peamust.csv")
passport[nrow(passport) + 1, "ivip"] <- "Hedin"
group_info<-read.csv("../11.infos_genotypes/metadata.csv", header = FALSE)
colnames(group_info)[1] <- "ivip"
colnames(group_info)[2] <- "group"
emmeans_fabae<-read.csv("../../../faba_specifics/02.analyses/01.GWAS/01_phenotipic_analysis/02.ajusted_means/fabae/emmeans_fabae_4.csv")
emmeans_fabae<-emmeans_fabae[,-c(3:6)]
colnames(emmeans_fabae)[2] <- "emmeans_fabae"
emmeans_pisum<-read.csv("../../../faba_specifics/02.analyses/01.GWAS/01_phenotipic_analysis/02.ajusted_means/pisum/emmeans_pisum_6.csv")
colnames(emmeans_pisum)[2] <- "emmeans_pisum"
passport_fabae <-merge(passport, emmeans_fabae, by = "ivip", all = TRUE)
passport_fabae_pisum<-merge(passport_fabae, emmeans_pisum, by = "ivip", all = TRUE)
passport_fabae_pisum_group<-merge(passport_fabae_pisum, group_info, by = "ivip", all = TRUE)

write_xlsx(passport_fabae_pisum_group, "metadata_passportdata_emmeans.xlsx")
