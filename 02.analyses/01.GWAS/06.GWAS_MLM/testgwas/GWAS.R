testmyY <- read.table("http://zzlab.net/GAPIT/data/mdp_traits.txt", head = TRUE)
testmyGD=read.table(file="http://zzlab.net/GAPIT/data/mdp_numeric.txt",head=T)
testmyGM=read.table(file="http://zzlab.net/GAPIT/data/mdp_SNP_information.txt",head=T)
# install.packages("BiocManager")
# BiocManager::install("Biobase")
# BiocManager::install("multtest")
# install.packages("GAPIT3")

# install.packages("devtools")
# devtools::install_github("jiabowang/GAPIT3",force=TRUE)
library(data.table)
library(vcfR)   
library(tidyverse)



pheno<-read.csv("02.analyses/01.GWAS/01_phenotipic_analysis/02.ajusted_means/emmeans_fabae.csv", sep = ",")        
myGD1 <-fread(file = "02.analyses/01.GWAS/06.GWAS_MLM/testgwas/PEAMUST_vfaba_F1_SNP_imputed_MAF5.012",header = FALSE, stringsAsFactors =  FALSE)%>% select(-1)



myGD2 <-fread(file = "02.analyses/01.GWAS/06.GWAS_MLM/testgwas/PEAMUST_vfaba_F1_SNP_imputed_MAF5.012.indv",header = FALSE,stringsAsFactors =  FALSE)

myKI <- fread(file = "02.analyses/01.GWAS/06.GWAS_MLM/testgwas/PEAMUST_vfaba_F1_SNP_imputed_MAF5_kinship.cXX.txt")
# Load datasets



myGM <-fread(file = "02.analyses/01.GWAS/06.GWAS_MLM/testgwas/PEAMUST_vfaba_F1_SNP_imputed_MAF5.012.pos") %>% rename(
    "chromosome" = "V1",
    "position" = "V2"
  )%>% mutate(snp_id=paste(chromosome,position, sep = "_")) %>% relocate(snp_id,.before = chromosome) 

colnames(myGD2)[1]<-"ivip"
myY <- merge(myGD2, pheno, by = "ivip", all.x = TRUE)

colnames(myGD1)<- myGM$snp_id
myGD <-bind_cols(myGD2,myGD1)

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
  model=c("GLM", "MLM", "CMLM", "FarmCPU", "Blink"),
  Multiple_analysis=TRUE)

