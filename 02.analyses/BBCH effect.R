










#emmeans even when we remove genetic effect it seem like we have an effect of the BBCH 25 on the larvae count ...----
model2<- lmer(larvae_total_number~bbch_10days + (1 |genotype), data = combined_data)
model2b<- lmer(larvae_total_number~bbch_25days + (1 |genotype), data = combined_data)

emm_genotypeBBCH10 <- emmeans(model2, ~ bbch_10days)
emm_dfBBCH10 <- as.data.frame(emm_genotypeBBCH10)

emm_genotypebbch25 <- emmeans(model2b, ~ bbch_25days)
emm_dfBBCH25 <- as.data.frame(emm_genotypebbch25)


