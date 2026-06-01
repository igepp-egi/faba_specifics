
library(data.table)
library(vcfR)   
library(tidyverse)
library(dplyr)


###### GWAS A.fabae ######

pheno<-read.csv("Documents/gitrepo/faba_specifics/02.analyses/01.GWAS/01_phenotipic_analysis/02.ajusted_means/fabae/emmeans_fabae_4.csv", sep = ",")        

myGD1 <-fread(file = "Documents/results_with_hedin/vfaba_afabae_n4.filtered_imputed_maf.vcf.hedin.vcf.012",header = FALSE, stringsAsFactors =  FALSE) #%>% select(-1)
myGD1 <- myGD1[, -1, with = FALSE]  # Removes the first column

myGD2 <-fread(file = "Documents/results_with_hedin/vfaba_afabae_n4.filtered_imputed_maf.vcf.hedin.vcf.012.indv",header = FALSE, stringsAsFactors =  FALSE) #%>% select(-1)
myGM <-fread(file = "Documents/results_with_hedin/vfaba_afabae_n4.filtered_imputed_maf.vcf.hedin.vcf.012.pos") %>% rename(
  "chromosome" = "V1",
  "position" = "V2"
)%>% mutate(snp_id=paste(chromosome,position, sep = "_")) %>% relocate(snp_id,.before = chromosome) 

colnames(myGD2)[1]<-"ivip"
myY <- merge(myGD2, pheno, by = "ivip", all.x = TRUE)

colnames(myGD1)<- myGM$snp_id
myGD <-bind_cols(myGD2,myGD1)


myKI <- fread(file = "Documents/results_with_hedin/output/vfaba_GRM_fabae.cXX.txt")

# Add snp_ID column from dataset2 as the first column of dataset1
myKI <- cbind(snp_id = myY$ivip, myKI)

source("http://zzlab.net/GAPIT/GAPIT.library.R")
source("http://zzlab.net/GAPIT/gapit_functions.txt")


myGAPIT_MLM <- GAPIT(
  Y=myY[,c(1,2)], 
  GD=myGD,
  GM=myGM,
  KI = myKI,
  PCA.total=3,
  model=c("MLMM", "Blink"),
  Multiple_analysis=TRUE)


###### GWAS A.pisum ######

pheno<-read.csv("Documents/gitrepo/faba_specifics/02.analyses/01.GWAS/01_phenotipic_analysis/02.ajusted_means/pisum/emmeans_pisum_4.csv", sep = ",")        

myGD1 <-fread(file = "Documents/results_with_hedin/vfaba_afabae_n4.filtered_imputed_maf.vcf.hedin.vcf.012",header = FALSE, stringsAsFactors =  FALSE) #%>% select(-1)
myGD1 <- myGD1[, -1, with = FALSE]  # Removes the first column

myGD2 <-fread(file = "Documents/results_with_hedin/vfaba_afabae_n4.filtered_imputed_maf.vcf.hedin.vcf.012.indv",header = FALSE, stringsAsFactors =  FALSE) #%>% select(-1)
myGM <-fread(file = "Documents/results_with_hedin/vfaba_afabae_n4.filtered_imputed_maf.vcf.hedin.vcf.012.pos") %>% rename(
  "chromosome" = "V1",
  "position" = "V2"
)%>% mutate(snp_id=paste(chromosome,position, sep = "_")) %>% relocate(snp_id,.before = chromosome) 

colnames(myGD2)[1]<-"ivip"
myY <- merge(myGD2, pheno, by = "ivip", all.x = TRUE)

colnames(myGD1)<- myGM$snp_id
myGD <-bind_cols(myGD2,myGD1)


myKI <- fread(file = "Documents/results_with_hedin/output/vfaba_GRM_fabae.cXX.txt")

# Add snp_ID column from dataset2 as the first column of dataset1
myKI <- cbind(snp_id = myY$ivip, myKI)

source("http://zzlab.net/GAPIT/GAPIT.library.R")
source("http://zzlab.net/GAPIT/gapit_functions.txt")


myGAPIT_MLM <- GAPIT(
  Y=myY[,c(1,2)], 
  GD=myGD,
  GM=myGM,
  KI = myKI,
  PCA.total=3,
  model=c("MLMM", "Blink"),
  Multiple_analysis=TRUE)

