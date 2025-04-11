library(dplyr)





dataset1<-read.csv("01.data/2023_screen_specifics_251genotypes/01.acyrthosiphon_pisum_P123amp_on_faba_collection/acyrthospipon_pisum_fecundity_on_faba_collection.csv", sep=",")
dataset1<-read.csv("01.data/2023_screen_specifics_251genotypes/02.aphis_fabae_on_faba_collection/aphis_fabae_fecundity_on_faba_collection_2024-04-02.csv", sep=",")

dataset2 <-read.csv("00.info/00.seed_stocks_and_multiplication/morphological_data.csv",sep = ",")

dataset2 <- dataset2 %>%
  rename(ivip = CODE.IVIP)

dataset1 <- dataset1 %>%
  left_join(dataset2 %>% select(ivip, size), by = "ivip")

dataset1$size<-as.numeric(dataset1$size)

summary(dataset1)

median <- dataset1 %>%
  group_by(ivip) %>%
  summarise(Median_fecondity = median(larvae_total_number))
median <- median %>%
  left_join(dataset2 %>% select(ivip, size), by = "ivip")


mod<-lm(larvae_total_number ~ ivip, data=dataset1)
mod<-lm(size ~ ivip , data=dataset1)

anova(mod)

# Basic scatterplot
plot(x = dataset1$larvae_total_number, y = dataset1$size,
     main = "Scatterplot Title",
     xlab = "larvae_total_number", ylab = "plant size",
     pch = 19, col = "blue")

plot(x = median$Median_fecondity, y = median$size,
     main = "Scatterplot Title",
     xlab = "larvae_total_number", ylab = "plant size",
     pch = 10, col = "blue")
